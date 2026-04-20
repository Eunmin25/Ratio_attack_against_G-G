#include "GG_scheme.h"
#include "GG_sampler.h"
#include "GG_klein.h"
#include "GG_poly.h"
#include <openssl/evp.h>
#include <cmath>
#include <complex>
#include <ctime>
#include <cstdlib>
#include <fstream>

namespace GG {

#ifndef GG_DEBUG
#define GG_DEBUG 0
#endif

#if GG_DEBUG
#define GG_DEBUG_DO(stmt) do { stmt; } while (0)
#else
#define GG_DEBUG_DO(stmt) do { } while (0)
#endif

static bool env_flag_enabled_local(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

static int env_int_or_default_local(const char* name, int default_value) {
    const char* value = std::getenv(name);
    if (value == nullptr || value[0] == '\0') return default_value;
    return std::atoi(value);
}

static double env_double_or_default_local(const char* name, double default_value) {
    const char* value = std::getenv(name);
    if (value == nullptr || value[0] == '\0') return default_value;
    return std::atof(value);
}

static std::complex<double> eval_poly_at_negacyclic_root(const PolyInt& p, int freq_idx) {
    const double pi = std::acos(-1.0);
    const double theta = pi * (2.0 * (double)freq_idx + 1.0) / (double)GG_N;
    const std::complex<double> xi(std::cos(theta), std::sin(theta));

    std::complex<double> val(0.0, 0.0);
    std::complex<double> pow_xi(1.0, 0.0);
    for (int t = 0; t < GG_N; ++t) {
        val += (double)p.a[t] * pow_xi;
        pow_xi *= xi;
    }
    return val;
}

// Sigma_1 for S built from stacked skew-circulant blocks rot(s_i), i=0..k-1.
// We intentionally use s (no zeta multiplier) as requested.
static double sigma1_skew_circulant_stack(const Poly2QVec& s_vec) {
    if (s_vec.empty()) return 0.0;

    std::vector<PolyInt> lifted;
    lifted.reserve(s_vec.size());
    for (const auto& p2q : s_vec) {
        lifted.push_back(polyint_from_poly2q_centered(p2q));
    }

    double sigma1 = 0.0;
    for (int j = 0; j < GG_N; ++j) {
        double energy = 0.0;
        for (const auto& p : lifted) {
            const std::complex<double> v = eval_poly_at_negacyclic_root(p, j);
            energy += std::norm(v);
        }
        const double sj = std::sqrt(energy);
        if (sj > sigma1) sigma1 = sj;
    }
    return sigma1;
}

void keygen(PublicKey& pk, SecretKey& sk) {
    const double sigma_base_for_rejection = env_double_or_default_local(
        "GG_SIGMA_BASE_OVERRIDE", GG_SIGMA_BASE);
    const double sigma_u_for_rejection = GG_SIGMA_U;
    const double sqrt2 = std::sqrt(2.0);
    const double derived_threshold =
        (sigma_u_for_rejection > 0.0)
            ? (sigma_base_for_rejection / (sqrt2 * sigma_u_for_rejection))
            : 0.0;
    const double threshold_override = env_double_or_default_local(
        "GG_KEYGEN_S_THRESHOLD_OVERRIDE", 0.0);
    const double s_threshold = (threshold_override > 0.0) ? threshold_override : derived_threshold;
    const int max_retries_env = env_int_or_default_local(
        "GG_KEYGEN_MAX_RETRIES", GG_KEYGEN_MAX_RETRIES);
    const int max_retries = (max_retries_env > 0) ? max_retries_env : GG_KEYGEN_MAX_RETRIES;
    const bool enforce_rejection = (s_threshold > 0.0);
    const bool log_sigma1 = env_flag_enabled_local("GG_KEYGEN_LOG_SIGMA1");

    const int k1 = GG_k - GG_m - 1;
    int attempt = 0;
    while (true) {
        ++attempt;
        // Generate s1, s2 as PolyVecs (논문: s1 ∈ χ_η^{k-m-1}, s2 ∈ χ_η^m)
        PolyIntVec s1_int = uniform_vectorint_bounded(k1, GG_ETA); // without mod q
        PolyIntVec s2_int = uniform_vectorint_bounded(GG_m, GG_ETA);

        // For all mod-q matrix operations below, we need Poly (mod q).
        PolyVec s1 = polyvec_from_intvec_modq(s1_int);
        PolyVec s2 = polyvec_from_intvec_modq(s2_int);

        // Paper (Alg.1): choose A0 <- U(R_q^{m x (k-m-1)})
        MatrixR A0 = MatrixR::random_uniform(GG_m, k1); // 내부적으로 Poly::random_uniform 호출, mod_q 연산과 관련

        // b^T := A0 s1^T + s2^T (mod q)  => b is length-m PolyVec
        PolyVec b(GG_m); // mod_q 연산과 관련
        for (int row = 0; row < GG_m; ++row) {
            Poly acc = s2[row];
            for (int col = 0; col < k1; ++col) {
                acc = acc + (A0.at(row, col) * s1[col]);
            }
            b[row] = acc;
        }

        // A := (-2b^T + q j^T | 2A0 | 2I_m)
        // -2b^T + q j^T
        pk.A = Matrix2Q(GG_m, GG_k);
        // Paper-style mod 2q version: (-2b + qj) mod 2q
        for (int row = 0; row < GG_m; ++row) {
            Poly2Q col0 = poly2q_from_poly_mod2q(b[row]);// Poly 타입이었던 b[row]를 Poly2Q 타입으로 변환
            col0.mul_scalar(-2);
            if (row == 0) {
                Poly2Q qj = Poly2Q::zero();
                qj.a[0] = mod_2q((int64_t)GG_q);
                col0 = col0 + qj;
            }
            pk.A.at(row, 0) = col0;
        }

        // middle block: 2A0
        for (int row = 0; row < GG_m; ++row) {
            for (int col = 0; col < k1; ++col) {
                Poly2Q v = poly2q_from_poly_mod2q(A0.at(row, col));
                v.mul_scalar(2);
                pk.A.at(row, 1 + col) = v;
            }
        }
        // last block: 2I_m
        for (int row = 0; row < GG_m; ++row) {
            for (int col = 0; col < GG_m; ++col) {
                Poly2Q v = Poly2Q::zero();
                if (row == col) v.a[0] = mod_2q(2);
                pk.A.at(row, 1 + k1 + col) = v;
            }
        }

        // s := (1 | s1 | s2) mod q
        sk.s.clear();
        const Poly2QVec s1_2q = poly2qvec_from_polyvec_mod2q(s1);
        const Poly2QVec s2_2q = poly2qvec_from_polyvec_mod2q(s2);

        sk.s.resize(1 + s1_2q.size() + s2_2q.size());
        sk.s[0] = Poly2Q::zero();
        sk.s[0].a[0] = mod_2q(1);
        for (int i = 0; i < (int)s1_2q.size(); ++i) sk.s[1 + i] = s1_2q[i];
        for (int i = 0; i < (int)s2_2q.size(); ++i) sk.s[1 + (int)s1_2q.size() + i] = s2_2q[i];

        const double sigma1_s = sigma1_skew_circulant_stack(sk.s);
        if (log_sigma1) {
            std::cerr << "[keygen] attempt=" << attempt
                      << ", sigma1(rot(s))=" << sigma1_s;
            if (enforce_rejection) {
                std::cerr << ", threshold=" << s_threshold
                          << ", sigma_base=" << sigma_base_for_rejection
                          << ", sigma_u=" << sigma_u_for_rejection
                          << ", rule=sigma >= sqrt(2)*sigma_u*sigma1";
            }
            std::cerr << "\n";
        }

        if (!enforce_rejection || sigma1_s < s_threshold) {
            break;
        }
        if (attempt >= max_retries) {
            std::cerr << "[keygen] max retries reached (" << max_retries
                      << "), accepting last key with sigma1=" << sigma1_s << "\n";
            break;
        }
    }

    // ====== DEBUG: Print As and qj^T for verification ======
    // Compute As = pk.A * sk.s
    Poly2QVec As = mat_vec_mul_mod2q(pk.A, sk.s);
    // Compute qj^T (first basis vector times q)
    Poly2QVec qj(pk.A.rows, Poly2Q::zero());
    // Only the constant term of the first row is q, all others are 0
    if (pk.A.rows > 0) {
        qj[0].a[0] = mod_2q((int64_t)GG_q);
    }
    /*
    std::cout << "\n[DEBUG] As (A * s):\n";
    for (size_t row = 0; row < As.size(); ++row) {
        std::cout << "  As[" << row << "]: [";
        for (int i = 0; i < GG_N; ++i) {
            std::cout << As[row].a[i];
            if (i < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "[DEBUG] qj^T :\n";
    for (size_t row = 0; row < qj.size(); ++row) {
        std::cout << "  qj[" << row << "]: [";
        for (int i = 0; i < GG_N; ++i) {
            std::cout << qj[row].a[i];
            if (i < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "====================================================\n";*/
}



Signature sign(const PublicKey& pk, const SecretKey& sk, const std::string& message) {
    // Compute sigma for covariance matrix Σ(S) = σ² I_m - s² SS^T
    double sigma_base = env_double_or_default_local("GG_SIGMA_BASE_OVERRIDE", GG_SIGMA_BASE);
    
    // 1: Sample masking vector y ~ D_{Z^m, Σ(S)} using Klein-style sampling
    //    Σ(S) = σ² I_m - s² SS^T where S = rot(ζs)
    //    Klein sampling ensures y has the correct covariance structure
    double sigma_u = GG_SIGMA_U; // (혹은 적절한 값/상수)
    PolyIntVec y = sample_y_klein(sk, sigma_base, sigma_u);
    // Always print y for debugging
    //print_polyintvec_full(y, "[DEBUG] Masking vector y (Klein sampled)");


    // 2: Compute w = A * y (mod q)
    Poly2QVec w = mat_vec_mul_mod2q(pk.A, y);

    /*// Debug: print hash input (w, message) in sign
    std::cout << "[sign] hash input (w, message):\n";
    for (size_t row = 0; row < w.size(); ++row) {
        std::cout << "  w[" << row << "]: [";
        for (int i = 0; i < GG_N; ++i) {
            std::cout << w[row].a[i];
            if (i < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "  message: " << message << "\n";*/

    // 3: Compute challenge c = H(w, message)
    PolyInt c = hash_to_challenge(w, message);
    //GG_DEBUG_DO(print_polyint(c, "[DEBUG] Challenge c"));

    // 4: Sample u <- D_{R, sigma_u^2 I_n, -c/2}
    //    Coefficient-wise discrete Gaussian centered at -c_i/2.
    PolyInt u;
    //std::cerr << "[DEBUG] u (sampled in sign) = [";
    for (int i = 0; i < GG_N; ++i) {
        const double center = -0.5 * (double)(c.a[i]);
        const int32_t ui = sample_gaussian_centered(center, sigma_u);
        u.a[i] = ui;
        //std::cerr << ui;
    }
    //std::cerr << "]\n";
    //GG_DEBUG_DO(print_polyint(u, "[DEBUG] Sampled u (integer coeffs)"));

    // 5: Print secret key s (full, but only once per process)
    /*static bool printed_secret_s = false;
    if (!printed_secret_s) {
        GG_DEBUG_DO(print_polyvec_full(polyvec_from_poly2qvec_modq(sk.s), "[DEBUG] Secret key s (mod q view)"));
        printed_secret_s = true;
    }*/

    // 6: Compute signature z := y + (2u + c)s
    //    Here z has the same dimension as s and y (length k).
    PolyInt two_u_plus_c = u;
    two_u_plus_c.mul_scalar(2);
    two_u_plus_c += c;
    //GG_DEBUG_DO(print_polyint(two_u_plus_c, "[DEBUG] (2u + c) mod q")); //Poly::mul_scalar 또는 Poly::+=에서 mod_q 연산 함

    PolyIntVec z = y;
    const int k = (int)sk.s.size();
    if ((int)z.size() != k) z.resize(k);

    const bool decomp_enabled = env_flag_enabled_local("GG_LOG_Z_DECOMP");
    const int decomp_target_poly = env_int_or_default_local("GG_DECOMP_TARGET_POLY", 1);
    const int decomp_target_coeff = env_int_or_default_local("GG_DECOMP_TARGET_COEFF", 0);
    const int decomp_max_samples = env_int_or_default_local("GG_DECOMP_MAX_SAMPLES", 20000);
    static long long decomp_sample_idx = 0;
    const long long sid = decomp_sample_idx++;
    const bool do_decomp_log = decomp_enabled && sid < decomp_max_samples &&
                               decomp_target_poly >= 0 && decomp_target_poly < k &&
                               decomp_target_coeff >= 0 && decomp_target_coeff < GG_N;
    static std::ofstream decomp_log;
    static bool decomp_log_initialized = false;
    if (decomp_enabled && !decomp_log_initialized) {
        decomp_log.open("z_decomp.csv", std::ios::out | std::ios::trunc);
        if (decomp_log.is_open()) {
            decomp_log << "sample,target_poly,target_coeff,"
                      << "s00,sij,"
                      << "y00,w00,z00,"
                      << "yij,wij,zij,"
                      << "ratio" << "\n";
        }
        decomp_log_initialized = true;
    }

    int64_t y00 = 0, w00 = 0, z00 = 0;
    int64_t yij = 0, wij = 0, zij = 0;
    int64_t s00 = 0, sij = 0;
    bool have00 = false, haveij = false;

    for (int i = 0; i < k; ++i) {
        const PolyInt si = polyint_from_poly2q_centered(sk.s[i]);
        const PolyInt w_part = (two_u_plus_c * si);
        z[i] = z[i] + w_part;
        if (do_decomp_log) {
            // Denominator is fixed to z_{0,0} (paper ratio target).
            if (i == 0) {
                y00 = y[i].a[0];
                w00 = w_part.a[0];
                z00 = z[i].a[0];
                s00 = si.a[0];
                have00 = true;
            }
            if (i == decomp_target_poly) {
                yij = y[i].a[(size_t)decomp_target_coeff];
                wij = w_part.a[(size_t)decomp_target_coeff];
                zij = z[i].a[(size_t)decomp_target_coeff];
                sij = si.a[(size_t)decomp_target_coeff];
                haveij = true;
            }
        }
    }

    if (do_decomp_log && have00 && haveij && decomp_log.is_open()) {
        const double ratio = (z00 != 0) ? ((double)zij / (double)z00) : 0.0;
        decomp_log << sid << "," << decomp_target_poly << "," << decomp_target_coeff << ","
                   << s00 << "," << sij << ","
                   << y00 << "," << w00 << "," << z00 << ","
                   << yij << "," << wij << "," << zij << ","
                   << ratio << "\n";
        decomp_log.flush();
    }
    //GG_DEBUG_DO(print_polyintvec(z, "[DEBUG] Signature z"));
    
    Signature sig;
    sig.z = z;
    sig.c = c;
    return sig;
}

bool verify(const PublicKey& pk, const std::string& message, const Signature& sig) {
    // Algorithm 3: Verification of the Generic G+G CGS
    // Input: message, pk = A, signature sig = (z, c)
    // Output: validity of the signature

    // 1. Compute v := Az^T - qcj^T mod 2q
    // Az^T
    Poly2QVec Az = mat_vec_mul_mod2q(pk.A, sig.z);
    // -qcj^T (j^T: first basis vector, only a[0]=1)
    Poly2QVec qcj(pk.A.rows, Poly2Q::zero());
    // Only the first row: constant term is c.a[i]*q, others are 0
    if (pk.A.rows > 0) {
        qcj[0] = Poly2Q::zero();
        for (int i = 0; i < GG_N; ++i) {
            qcj[0].a[i] = mod_2q(sig.c.a[i] * GG_q);
        }
    }
    // v = Az^T - qcj^T mod 2q
    Poly2QVec v(pk.A.rows);
    for (int row = 0; row < pk.A.rows; ++row) {
        v[row] = Az[row] - qcj[row];
        // mod 2q
        for (int i = 0; i < GG_N; ++i) {
            v[row].a[i] = mod_2q(v[row].a[i]);
        }
    }


    // Debug: print hash input (v, message) in verify
    /*std::cout << "[verify] hash input (v, message):\n";
    for (size_t row = 0; row < v.size(); ++row) {
        std::cout << "  v[" << row << "]: [";
        for (int i = 0; i < GG_N; ++i) {
            std::cout << v[row].a[i];
            if (i < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "  message: " << message << "\n";*/

    // Print challenge comparison
    PolyInt challenge = hash_to_challenge(v, message);
    /*std::cout << "[verify] challenge (recomputed): [";
    for (int i = 0; i < std::min(8, GG_N); ++i) {
        std::cout << challenge.a[i];
        if (i < std::min(8, GG_N) - 1) std::cout << ", ";
    }
    if (GG_N > 8) std::cout << ", ...";
    std::cout << "]\n";

    std::cout << "[verify] sig.c (from signature): [";
    for (int i = 0; i < std::min(8, GG_N); ++i) {
        std::cout << sig.c.a[i];
        if (i < std::min(8, GG_N) - 1) std::cout << ", ";
    }
    if (GG_N > 8) std::cout << ", ...";
    std::cout << "]\n";*/

    if (!polyint_equal(challenge, sig.c)) return false;

    // ||z|| <= gamma, gamma is defined as GG_BZ_BOUND.
    double norm = 0.0;
    for (const auto& poly : sig.z) {
        for (int i = 0; i < GG_N; ++i) {
            norm += (double)(poly.a[i]) * (double)(poly.a[i]);
        }
    }
    norm = std::sqrt(norm);
    if (norm > GG_BZ_BOUND) return false;

    return true;
}

PolyInt hash_to_challenge(const Poly2QVec& w, const std::string& message) {
    // Build a 0/1 challenge polynomial c from SHA-256 bits (paper-style).
    // We hash (w || message) to a 32-byte seed, then expand with SHA-256(seed || counter)
    // to get enough bits for GG_N coefficients.

    auto sha256 = [](const unsigned char* data, size_t len, unsigned char out[32]) {
        EVP_MD_CTX* ctx = EVP_MD_CTX_new();
        EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr);
        EVP_DigestUpdate(ctx, data, len);
        unsigned int out_len = 0;
        EVP_DigestFinal_ex(ctx, out, &out_len);
        EVP_MD_CTX_free(ctx);
    };

    // First digest seed = SHA256(serialize(w) || message)
    EVP_MD_CTX* ctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr);

    for (const auto& poly : w) {
        for (int i = 0; i < GG_N; ++i) {
            // Canonical serialization: map coefficient from (-q, q] to [0, 2q).
            // Prevents equivalent values (e.g., -1 and 2q-1) from hashing differently.
            int32_t t = poly.a[i];
            if (t < 0) t += (int32_t)((int64_t)GG_q * 2);
            uint16_t coeff = (uint16_t)t; // safe here because 2*GG_q = 24578 < 65536
            unsigned char bytes[2] = {
                (unsigned char)(coeff & 0xff),
                (unsigned char)(coeff >> 8)
            };
            EVP_DigestUpdate(ctx, bytes, 2);
        }
    }

    EVP_DigestUpdate(ctx, message.data(), message.size());

    unsigned char seed[32];
    unsigned int seed_len = 0;
    EVP_DigestFinal_ex(ctx, seed, &seed_len);
    EVP_MD_CTX_free(ctx);

    PolyInt challenge;

    unsigned char block[32];
    uint32_t counter = 0;
    int bits_left = 0;

    auto refill_block = [&]() {
        unsigned char buf[36];
        for (int i = 0; i < 32; ++i) buf[i] = seed[i];
        buf[32] = (unsigned char)(counter & 0xff);
        buf[33] = (unsigned char)((counter >> 8) & 0xff);
        buf[34] = (unsigned char)((counter >> 16) & 0xff);
        buf[35] = (unsigned char)((counter >> 24) & 0xff);
        sha256(buf, sizeof(buf), block);
        ++counter;
        bits_left = 256;
    };

    refill_block();
    for (int i = 0; i < GG_N; ++i) {
        if (bits_left == 0) {
            refill_block();
        }
        int bit_index = 256 - bits_left;
        unsigned char byte = block[bit_index / 8];
        int bit = (byte >> (bit_index % 8)) & 1;
        challenge.a[i] = bit; // exactly 0/1
        --bits_left;
    }

    return challenge;
}

//============================================================
// Utility Functions for Printing
//============================================================

void print_poly(const Poly& p, const std::string& name) {
    std::cout << name << " (first 8 coeffs): [";
    for (int i = 0; i < std::min(8, GG_N); ++i) {
        std::cout << p.a[i];
        if (i < std::min(8, GG_N) - 1) std::cout << ", ";
    }
    if (GG_N > 8) std::cout << ", ...";
    std::cout << "]\n";
}

void print_polyint(const PolyInt& p, const std::string& name) {
    std::cout << name << " (first 8 coeffs): [";
    for (int i = 0; i < std::min(8, GG_N); ++i) {
        std::cout << p.a[i];
        if (i < std::min(8, GG_N) - 1) std::cout << ", ";
    }
    if (GG_N > 8) std::cout << ", ...";
    std::cout << "]\n";
}

void print_polyvec(const PolyVec& v, const std::string& name) {
    std::cout << name << " (length " << v.size() << "):\n";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "  [" << i << "]: [";
        for (int j = 0; j < std::min(4, GG_N); ++j) {
            std::cout << v[i].a[j];
            if (j < std::min(4, GG_N) - 1) std::cout << ", ";
        }
        if (GG_N > 4) std::cout << ", ...";
        std::cout << "]\n";
    }
}

void print_polyvec_full(const PolyVec& v, const std::string& name) {
    std::cout << name << " (length " << v.size() << "):\n";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "  [" << i << "]: [";
        for (int j = 0; j < GG_N; ++j) {
            std::cout << v[i].a[j];
            if (j < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

void print_polyintvec(const PolyIntVec& v, const std::string& name) {
    std::cout << name << " (length " << v.size() << "):\n";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "  [" << i << "]: [";
        for (int j = 0; j < std::min(4, GG_N); ++j) {
            std::cout << v[i].a[j];
            if (j < std::min(4, GG_N) - 1) std::cout << ", ";
        }
        if (GG_N > 4) std::cout << ", ...";
        std::cout << "]\n";
    }
}

void print_polyintvec_full(const PolyIntVec& v, const std::string& name) {
    std::cout << name << " (length " << v.size() << "):\n";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "  [" << i << "]: [";
        for (int j = 0; j < GG_N; ++j) {
            std::cout << v[i].a[j];
            if (j < GG_N - 1) std::cout << ", ";
        }
        std::cout << "]\n";
    }
}

void print_matrix(const MatrixR& M, const std::string& name) {
    std::cout << name << " (" << M.rows << " x " << M.cols << "):\n";
    for (int i = 0; i < M.rows; ++i) {
        for (int j = 0; j < M.cols; ++j) {
            std::cout << "  A[" << i << "][" << j << "][0:3] = [";
            for (int k = 0; k < std::min(4, GG_N); ++k) {
                std::cout << M.at(i,j).a[k];
                if (k < std::min(4, GG_N) - 1) std::cout << ", ";
            }
            if (GG_N > 4) std::cout << ", ...";
            std::cout << "]\n";
        }
    }
}

void print_matrix(const Matrix2Q& M, const std::string& name) {
    std::cout << name << " (" << M.rows << " x " << M.cols << "):\n";
    for (int i = 0; i < M.rows; ++i) {
        for (int j = 0; j < M.cols; ++j) {
            std::cout << "  A[" << i << "][" << j << "][0:3] = [";
            for (int k = 0; k < std::min(4, GG_N); ++k) {
                std::cout << M.at(i, j).a[k];
                if (k < std::min(4, GG_N) - 1) std::cout << ", ";
            }
            if (GG_N > 4) std::cout << ", ...";
            std::cout << "]\n";
        }
    }
}

void print_public_key(const PublicKey& pk) {
    std::cout << "\n========== PUBLIC KEY ==========\n";
        // Print all coefficients for each entry in A
        std::cout << "Public Matrix A (full):\n";
        for (int i = 0; i < pk.A.rows; ++i) {
            for (int j = 0; j < pk.A.cols; ++j) {
                std::cout << "  A[" << i << "][" << j << "] = [";
                for (int k = 0; k < GG_N; ++k) {
                    std::cout << pk.A.at(i, j).a[k];
                    if (k < GG_N - 1) std::cout << ", ";
                }
                std::cout << "]\n";
            }
        }
    std::cout << "================================\n";
}

void print_secret_key(const SecretKey& sk) {
    std::cout << "\n========== SECRET KEY ==========\n";
        PolyVec s_modq = polyvec_from_poly2qvec_modq(sk.s);
        std::cout << "Secret Vector s (mod q view) (length " << s_modq.size() << "):\n";
        for (size_t i = 0; i < s_modq.size(); ++i) {
            std::cout << "  [" << i << "]: [";
            for (int j = 0; j < GG_N; ++j) {
                std::cout << s_modq[i].a[j];
                if (j < GG_N - 1) std::cout << ", ";
            }
            std::cout << "]\n";
        }
    std::cout << "================================\n";
}

void print_signature(const Signature& sig) {
    std::cout << "\n========== SIGNATURE ==========\n";
    print_polyintvec(sig.z, "Response z");
    std::cout << "\n";
    print_polyint(sig.c, "Challenge c");
    std::cout << "================================\n";
}

} // namespace GG

