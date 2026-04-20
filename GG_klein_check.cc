#include "GG_scheme.h"
#include "GG_klein.h"
#include "GG_params.h"
#include "GG_poly.h"
#include "GG_sampler.h"

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

using namespace GG;

static int env_int_or_default(const char* name, int default_value) {
    const char* value = std::getenv(name);
    if (value == nullptr || value[0] == '\0') return default_value;
    return std::atoi(value);
}

// Mirror of `circ_poly` in GG_klein.cc (used only for expected covariance computation here).
static std::vector<std::vector<double>> circ_poly_expected(const Poly2Q& s) {
    const int N = GG_N;
    std::vector<std::vector<double>> mat((size_t)N, std::vector<double>((size_t)N, 0.0));
    for (int row = 0; row < N; ++row) {
        for (int col = 0; col < N; ++col) {
            const int idx = (row - col + N) % N;
            const int32_t c = s.a[(size_t)idx];
            double val = (double)c;
            if (row < col) val = -val;
            mat[(size_t)row][(size_t)col] = val;
        }
    }
    return mat;
}

int main() {
    const int num_samples = env_int_or_default("GG_CHECK_SAMPLES", 200);
    const int coeff_index = env_int_or_default("GG_CHECK_COEFF", 0);
    if (coeff_index < 0 || coeff_index >= GG_N) {
        std::cerr << "GG_CHECK_COEFF must be in [0," << GG_N - 1 << "]\n";
        return 1;
    }

    const double sigma_base = GG_SIGMA_BASE;
    const double sigma_u = GG_SIGMA_U;

    PublicKey pk;
    SecretKey sk;
    keygen(pk, sk);

    const int k = (int)sk.s.size();

    // Collect empirical samples of y_i[coeff_index] for i=0..k-1.
    std::vector<std::vector<double>> samples((size_t)k, std::vector<double>((size_t)num_samples, 0.0));
    for (int t = 0; t < num_samples; ++t) {
        PolyIntVec y = sample_y_klein(sk, sigma_base, sigma_u);
        for (int i = 0; i < k; ++i) {
            samples[(size_t)i][(size_t)t] = (double)y[(size_t)i].a[(size_t)coeff_index];
        }
    }

    // Compute empirical mean and covariance among i=0..k-1 for the fixed coefficient index.
    std::vector<double> mean((size_t)k, 0.0);
    for (int i = 0; i < k; ++i) {
        double s = 0.0;
        for (int t = 0; t < num_samples; ++t) s += samples[(size_t)i][(size_t)t];
        mean[(size_t)i] = s / (double)num_samples;
    }

    std::vector<std::vector<double>> cov_emp((size_t)k, std::vector<double>((size_t)k, 0.0));
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            double acc = 0.0;
            for (int t = 0; t < num_samples; ++t) {
                const double di = samples[(size_t)i][(size_t)t] - mean[(size_t)i];
                const double dj = samples[(size_t)j][(size_t)t] - mean[(size_t)j];
                acc += di * dj;
            }
            cov_emp[(size_t)i][(size_t)j] = (num_samples > 1) ? (acc / (double)(num_samples - 1)) : 0.0;
        }
    }

    // Compute expected covariance using the same Σ construction used in GG_klein.cc:
    //   Σ = σ_base^2 * I - σ_u^2 * (S S^T)
    // where S is built from negacyclic circulant blocks from secret coefficients.
    std::vector<std::vector<double>> SSt((size_t)k * (size_t)GG_N, std::vector<double>((size_t)k * (size_t)GG_N, 0.0));
    for (int br = 0; br < k; ++br) {
        auto circ_i = circ_poly_expected(sk.s[(size_t)br]);
        for (int bc = 0; bc < k; ++bc) {
            auto circ_j = circ_poly_expected(sk.s[(size_t)bc]);
            for (int i = 0; i < GG_N; ++i) {
                for (int j = 0; j < GG_N; ++j) {
                    double sum = 0.0;
                    for (int t = 0; t < GG_N; ++t) {
                        // Match GG_klein.cc:
                        //   circ_j_T[t][j] = circ_j[j][t]
                        // so sum += circ_i[i][t] * circ_j_T[t][j] = circ_i[i][t] * circ_j[j][t]
                        sum += circ_i[(size_t)i][(size_t)t] * circ_j[(size_t)j][(size_t)t];
                    }
                    SSt[(size_t)br * (size_t)GG_N + (size_t)i][(size_t)bc * (size_t)GG_N + (size_t)j] = sum;
                }
            }
        }
    }

    auto expected_cov_entry = [&](int i, int j) -> double {
        const int idx_i = i * GG_N + coeff_index;
        const int idx_j = j * GG_N + coeff_index;
        double val = 0.0;
        if (idx_i == idx_j) val = sigma_base * sigma_base;
        // subtract sigma_u^2 * SSt
        val -= sigma_u * sigma_u * SSt[(size_t)idx_i][(size_t)idx_j];
        return val;
    };

    // Print results.
    std::cout << "[GG_klein_check] num_samples=" << num_samples
              << ", coeff_index=" << coeff_index
              << ", k=" << k << "\n";
    std::cout << "mean(y_i[coeff]) (empirical):\n";
    for (int i = 0; i < k; ++i) {
        std::cout << "  i=" << i << " mean=" << mean[(size_t)i] << "\n";
    }

    std::cout << "\nCovariance among y_i[coeff] (empirical vs expected):\n";
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < k; ++j) {
            const double ce = cov_emp[(size_t)i][(size_t)j];
            const double cx = expected_cov_entry(i, j);
            std::cout << "  cov(" << i << "," << j << ") emp=" << ce << ", exp=" << cx << "\n";
        }
    }
    return 0;
}

