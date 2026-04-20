#include "GG_klein.h"
#include "GG_scheme.h"  // for SecretKey definition
#include "GG_sampler.h"
#ifdef q0
#undef q0
#endif
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>

namespace GG {


static bool should_dump_circ() {
    const char* v = std::getenv("GG_DUMP_CIRC");//s의 circ 행렬이 찍히는데 모드 연산은 X
    return v != nullptr && v[0] != '\0' && v[0] != '0';
}

static bool should_dump_circ_full() {
    const char* v = std::getenv("GG_DUMP_CIRC_FULL");
    return v != nullptr && v[0] != '\0' && v[0] != '0';
}

static const char* dump_circ_full_path() {
    const char* v = std::getenv("GG_DUMP_CIRC_FULL_PATH");
    if (v != nullptr && v[0] != '\0') return v;
    return "S_full.csv";
}

static bool should_dump_chol() {
    const char* v = std::getenv("GG_DUMP_CHOL");
    return v != nullptr && v[0] != '\0' && v[0] != '0';
}

static bool env_flag_enabled(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

static void dump_matrix_preview(
    const std::vector<std::vector<double>>& M,
    const char* name,
    int max_rows = 16,
    int max_cols = 16
) {
    const int rows = (int)M.size();
    const int cols = rows ? (int)M[0].size() : 0;
    std::cerr << "[GG_DUMP_CIRC] " << name << " (" << rows << "x" << cols << ")\n";
    const int rlim = std::min(rows, max_rows);
    const int clim = std::min(cols, max_cols);
    std::cerr.setf(std::ios::fixed);
    std::cerr << std::setprecision(2);
    for (int i = 0; i < rlim; ++i) {
        for (int j = 0; j < clim; ++j) {
            std::cerr << std::setw(7) << M[i][j] << (j + 1 == clim ? "" : " ");
        }
        if (clim < cols) std::cerr << " ...";
        std::cerr << "\n";
    }
    if (rlim < rows) std::cerr << "...\n";
}

static void dump_matrix_csv(
    const std::vector<std::vector<double>>& M,
    const char* path
) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "[GG_DUMP_CIRC] failed to open CSV: " << path << "\n";
        return;
    }
    const int rows = (int)M.size();
    const int cols = rows ? (int)M[0].size() : 0;
    out.setf(std::ios::fixed);
    out << std::setprecision(10);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            out << M[i][j];
            if (j + 1 != cols) out << ',';
        }
        out << "\n";
    }
}

// 1. circ_poly: N x N skew-circulant (negacyclic) matrix from Poly
static std::vector<std::vector<double>> circ_poly(const Poly2Q& s) {
    std::vector<std::vector<double>> mat(
        GG_N, 
        std::vector<double>(GG_N, 0.0)
    );
    for (int row = 0; row < GG_N; ++row) {
        for (int col = 0; col < GG_N; ++col) {
            int idx = (row - col + GG_N) % GG_N;
            //const int32_t c = centered_coeff(s.a[idx]);
            const int32_t c = s.a[idx];
            double val = (double)c;
            if (row < col) val = -val;
            mat[row][col] = val;
        }
    }
    return mat;
}

// 2. circ_vec: kN x kN block skew-circulant matrix from PolyVec
static std::vector<std::vector<double>> circ_vec(const Poly2QVec& s_vec) {
    int k = s_vec.size();
    std::vector<std::vector<double>> mat(
        k * GG_N, 
        std::vector<double>(k * GG_N, 0.0)
    );
    for (int block_row = 0; block_row < k; ++block_row) {
        for (int block_col = 0; block_col < k; ++block_col) {
            int idx = (block_row - block_col + k) % k;
            auto circ = circ_poly(s_vec[idx]);
            for (int i = 0; i < GG_N; ++i) {
                for (int j = 0; j < GG_N; ++j) {
                    mat[block_row * GG_N + i][block_col * GG_N + j] = circ[i][j];
                }
            }
        }
    }

    if (should_dump_circ()) {
        // dump each unique circ_poly(s_vec[t]) once
        for (int t = 0; t < k; ++t) {
            auto Ct = circ_poly(s_vec[t]);
            std::string nm = std::string("circ_poly(s[") + std::to_string(t) + "])";
            dump_matrix_preview(Ct, nm.c_str());
        }
        dump_matrix_preview(mat, "circ_vec(s) (full S)", 80, 80);

        if (should_dump_circ_full()) {
            const char* path = dump_circ_full_path();
            dump_matrix_csv(mat, path);
            std::cerr << "[GG_DUMP_CIRC] wrote full S CSV: " << path << "\n";
        }
    }
    return mat;
}

// 3. compute_covariance_kn: kN x kN covariance matrix
// Σ = σ² I_{kN} - σ_u² S S^T, S = circ_vec(s)
static std::vector<std::vector<double>> compute_covariance_kn(
    const Poly2QVec& s_vec, double sigma, double sigma_u
) {
    int k = s_vec.size();//k일거임
    int dim = k * GG_N;
    /*
    auto S = circ_vec(s_vec); // PolyVec -> MatrixVec(kN x kN)
    // S 행렬 전체를 파일로 출력
    {
        int length = (int)S.size();
        std::ofstream s_log("S_full.log");
        s_log << "[DEBUG] S (circ_vec) full matrix (" << length << "x" << length << "):\n";
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                s_log << S[i][j];
                if (j + 1 != length) s_log << ' ';
            }
            s_log << "\n";
        }
        s_log.close();
    }*/
    // S S^T 계산 (블록 구조: 각 블록은 circ(s_i) * circ(s_j)^T)
    std::vector<std::vector<double>> SSt(dim, std::vector<double>(dim, 0.0));
    for (int block_row = 0; block_row < k; ++block_row) {
        for (int block_col = 0; block_col < k; ++block_col) {
            // circ(s_i), circ(s_j)
            auto circ_i = circ_poly(s_vec[block_row]);
            auto circ_j = circ_poly(s_vec[block_col]);
            // circ_j^T 계산
            std::vector<std::vector<double>> circ_j_T(GG_N, std::vector<double>(GG_N, 0.0));
            for (int r = 0; r < GG_N; ++r)
                for (int c = 0; c < GG_N; ++c)
                    circ_j_T[r][c] = circ_j[c][r];
            // 블록 곱셈 및 SSt에 할당
            for (int i = 0; i < GG_N; ++i) {
                for (int j = 0; j < GG_N; ++j) {
                    double sum = 0.0;
                    for (int t = 0; t < GG_N; ++t)
                        sum += circ_i[i][t] * circ_j_T[t][j];
                    SSt[block_row * GG_N + i][block_col * GG_N + j] = sum;
                }
            }
        }
    }
    // SSt 행렬 전체를 파일로 출력
    {
        int length = (int)SSt.size();
        std::ofstream sst_log("SSt_full.log");
        sst_log << "[DEBUG] SSt (S S^T) full matrix (" << length << "x" << length << "):\n";
        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j) {
                sst_log << SSt[i][j];
                if (j + 1 != length) sst_log << ' ';
            }
            sst_log << "\n";
        }
        sst_log.close();
    }
    // 디버깅용: S, S^T, SSt 일부 출력
    std::cerr << "[DEBUG] S S^T (SSt) preview (first 8x8):\n";
    for (int i = 0; i < std::min(8, dim); ++i) {
        for (int j = 0; j < std::min(8, dim); ++j) {
            std::cerr << std::setw(12) << SSt[i][j];
        }
        std::cerr << "\n";
    }


    // Σ = σ² I_{kN} - σ_u² S S^T
    std::vector<std::vector<double>> Sigma(dim, std::vector<double>(dim, 0.0));//kN x kN 행렬 Sigma 선언
    for (int i = 0; i < dim; ++i) {
        Sigma[i][i] = sigma * sigma;
        for (int j = 0; j < dim; ++j) {
            Sigma[i][j] -= sigma_u * sigma_u * SSt[i][j];
        }
    }
        // 디버깅용: Sigma 행렬 전체를 파일로 출력
    int length = (int)Sigma.size();
    std::ofstream sigma_log("Sigma_full.log");
    sigma_log << "[DEBUG] Sigma (covariance) full matrix (" << length << "x" << length << "):\n";
    for (int i = 0; i < length; ++i) {
        for (int j = 0; j < length; ++j) {
            sigma_log << Sigma[i][j];
            if (j + 1 != length) sigma_log << ' ';
        }
        sigma_log << "\n";
    }
    sigma_log.close();
    return Sigma;
}

//============================================================
// Linear algebra helpers (small-dim SPD)
//============================================================

// A = L L^T (L lower-triangular). Returns false if not (numerically) SPD.
static bool cholesky_lower(
    const std::vector<std::vector<double>>& A,
    std::vector<std::vector<double>>& L
) {
    int n = (int)A.size();
    Eigen::MatrixXd mat(n, n); //Eigen 라이브러리의 Eigen::MatrixXd 타입으로 타입변환
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            mat(i, j) = A[i][j];

    Eigen::LLT<Eigen::MatrixXd> llt(mat); // Cholesky 분해 시도
    if (llt.info() != Eigen::Success) { //분해 성공 여부를 확인하여
        //디버깅용이고,,
        std::cerr << "[GG_CHOL] Eigen Cholesky failed (matrix not SPD)\n";
        // Print A
        std::cerr << "[GG_CHOL] Matrix A (input) preview (first 8x8):\n";
        for (int i = 0; i < std::min(8, n); ++i) {
            for (int j = 0; j < std::min(8, n); ++j) {
                std::cerr << std::setw(12) << mat(i, j);
            }
            std::cerr << "\n";
        }
        // 성공하지 않았다면 false 반환하여 호출자에게 알림 + 밑에 코드 실행 안함
        return false;
    }
    //성공했다면 해당 코드를 통해 Lmat라는 행렬에 llt.matrixL()를 실행하여 하삼각 행렬을 반환
    Eigen::MatrixXd Lmat = llt.matrixL();
    L.assign(n, std::vector<double>(n, 0.0)); //하삼각행렬 Lmat를 
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j)
            L[i][j] = Lmat(i, j); // 2차원 벡터 L에 복사(Ls에 해당)

    // 이하 디버깅용
    // Always print A, L, L^T (first 8x8)
    std::cerr << "[GG_CHOL] Matrix A (input) preview (first 8x8):\n";
    for (int i = 0; i < std::min(8, n); ++i) {
        for (int j = 0; j < std::min(8, n); ++j) {
            std::cerr << std::setw(12) << mat(i, j);
        }
        std::cerr << "\n";
    }
    std::cerr << "[GG_CHOL] Cholesky L (first 8x8):\n";
    for (int i = 0; i < std::min(8, n); ++i) {
        for (int j = 0; j < std::min(8, n); ++j) {
            std::cerr << std::setw(12) << ((j <= i) ? L[i][j] : 0.0);
        }
        std::cerr << "\n";
    }
    std::cerr << "[GG_CHOL] Cholesky L^T (first 8x8):\n";
    for (int i = 0; i < std::min(8, n); ++i) {
        for (int j = 0; j < std::min(8, n); ++j) {
            std::cerr << std::setw(12) << ((i <= j) ? L[j][i] : 0.0);
        }
        std::cerr << "\n";
    }
    return true;
}

// Given L from Cholesky (A = L L^T), compute A^{-1}. <- 주석처리
static std::vector<std::vector<double>> spd_inverse_from_cholesky(
    const std::vector<std::vector<double>>& L
) {
    const int n = (int)L.size();

    // Solve L * Y = I
    std::vector<std::vector<double>> Y(n, std::vector<double>(n, 0.0));
    for (int col = 0; col < n; ++col) {
        for (int i = 0; i < n; ++i) {
            double rhs = (i == col) ? 1.0 : 0.0;
            for (int k = 0; k < i; ++k) rhs -= L[i][k] * Y[k][col];
            Y[i][col] = rhs / L[i][i];
        }
    }

    // Solve L^T * X = Y
    std::vector<std::vector<double>> X(n, std::vector<double>(n, 0.0));
    for (int col = 0; col < n; ++col) {
        for (int i = n - 1; i >= 0; --i) {
            double rhs = Y[i][col];
            for (int k = i + 1; k < n; ++k) rhs -= L[k][i] * X[k][col];
            X[i][col] = rhs / L[i][i];
        }
    }
    return X;
}

// Sample x ~ D_{Z^n, Sigma, 0} using conditionals from Q = Sigma^{-1} = R^T R.
static bool sample_discrete_gaussian_precision(
    std::vector<int32_t>& x,//메모리 주소를 받아서 여기에 샘플링된 값을 넣음
    const std::vector<std::vector<double>>& Sigma
) {
    static bool printed_mode = false;
    if (!printed_mode) {
        std::cerr << "[GG_KLEIN] sigma = sqrt(sigma2) (fixed)\n";
        printed_mode = true;
    }
    const int n = (int)Sigma.size();//n=kN
    if (n == 0) return false;
    // For each i, sample x[i] | x[i+1..n-1] using conditional mean/variance
    // Sigma is KN x KN covariance matrix
    x.assign(n, 0);
    std::vector<double> x_real(n, 0.0); // for conditional mean computation
    for (int i = n - 1; i >= 0; --i) {
        // Partition Sigma as:
        // [ a  b^T ]
        // [ b  C  ]
        // where a = Sigma[i][i], b = Sigma[i][i+1..n-1], C = Sigma[i+1..n-1][i+1..n-1]
        double a = Sigma[i][i];
        int m = n - i - 1;
        std::vector<double> b(m, 0.0);
        std::vector<std::vector<double>> C(m, std::vector<double>(m, 0.0));
        for (int j = 0; j < m; ++j) {
            b[j] = Sigma[i][i + 1 + j];
            for (int k = 0; k < m; ++k) {
                C[j][k] = Sigma[i + 1 + j][i + 1 + k];
            }
        }
        // Compute C^{-1} * (x_{i+1..n-1})
        std::vector<double> x_sub(m, 0.0);
        for (int j = 0; j < m; ++j) x_sub[j] = x_real[i + 1 + j];
        // Invert C (if m > 0)
        std::vector<std::vector<double>> Cinv(m, std::vector<double>(m, 0.0));
        if (m > 0) {
            // Use Eigen for inversion
            Eigen::MatrixXd Ceig(m, m);
            for (int j = 0; j < m; ++j)
                for (int k = 0; k < m; ++k)
                    Ceig(j, k) = C[j][k];
            Eigen::MatrixXd Cinv_eig = Ceig.inverse();
            for (int j = 0; j < m; ++j)
                for (int k = 0; k < m; ++k)
                    Cinv[j][k] = Cinv_eig(j, k);

            // Debug: Check if Ceig * Cinv_eig is close to identity
            Eigen::MatrixXd prod = Ceig * Cinv_eig;
            double max_abs_err = 0.0;
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < m; ++k) {
                    double expected = (j == k) ? 1.0 : 0.0; //대각행렬은 1, 나머지는 0이여야함.
                    double err = std::abs(prod(j, k) - expected);
                    if (err > max_abs_err) max_abs_err = err;
                }
            }
            //오차를 1X10^-8로 잡았는데, 이보다 크면 디버깅 메시지 출력
            if (max_abs_err > 1e-8) { 
                std::cerr << "[DEBUG] C * Cinv max abs error: " << max_abs_err << " (m=" << m << ")\n";
            }
        }
        // Compute conditional mean: mu = -b^T C^{-1} x_sub
        double mu = 0.0;
        if (m > 0) {
            for (int j = 0; j < m; ++j) {
                double tmp = 0.0;
                for (int k = 0; k < m; ++k) tmp += Cinv[j][k] * x_sub[k];
                mu += b[j] * tmp;//-=에서 +=로 수정(0407_1628)
            }
        }
        // Compute conditional variance: sigma2 = a - b^T C^{-1} b
        double sigma2 = a;
        if (m > 0) {
            for (int j = 0; j < m; ++j) {
                double tmp = 0.0;
                for (int k = 0; k < m; ++k) tmp += Cinv[j][k] * b[k];
                sigma2 -= b[j] * tmp;
            }
        }
        if (!(sigma2 > 0.0)) {
            std::cerr << "[COND_SAMPLER] Non-positive conditional variance at i=" << i << ": " << sigma2 << "\n";
            return false;
        }
        const double sigma = std::sqrt(sigma2);
        x[i] = sample_gaussian_centered(mu, sigma);
        x_real[i] = (double)x[i];
    }
    return true;
}





//============================================================
// Sample Full Masking Vector y with Klein Sampling
//============================================================

PolyIntVec sample_y_klein(const SecretKey& sk, double sigma_base, double sigma_u) {
    // 1) Build covariance Sigma (kN x kN)
    const auto Sigma = compute_covariance_kn(sk.s, sigma_base, sigma_u);
    const int k = (int)sk.s.size();
    const int N = GG_N;
    const int nk = k * N;

    // 2) Sample from D_{Z^{nk}, Sigma, 0} using precision factor conditionals
    std::vector<int32_t> y_vec; // nk-dimensional vector임
    if (!sample_discrete_gaussian_precision(y_vec, Sigma)) {
        const bool strict_klein = env_flag_enabled("GG_STRICT_KLEIN");
        if (strict_klein) {
            throw std::runtime_error("[GG_STRICT_KLEIN] discrete sampler failed; fallback disabled");
        }
        std::cerr << "Warning: discrete sampler failed; using spherical Gaussian instead\n";
        std::ofstream tty("/dev/tty");
        if (tty) tty << "Warning: discrete sampler failed; using spherical Gaussian instead\n";
        y_vec.assign(nk, 0);
        for (int i = 0; i < nk; ++i) y_vec[i] = sample_gaussian(sigma_base);
    }

    // 3) Convert to PolyVec (mod q)
    PolyIntVec y(k); // 다항식 k개로 구성된 vector
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j) {
            y[i].a[j] = y_vec[i * N + j];
        }
    }
    return y;
}

} // namespace GG
