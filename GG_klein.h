#ifndef GG_KLEIN_H
#define GG_KLEIN_H

#include "GG_params.h"
#include "GG_poly.h"
#include <vector>

namespace GG {

// Forward declaration
struct SecretKey;

//============================================================
// Klein-style Gaussian Sampling for G+G
// Samples from D_{Z^m, Σ(S)} where Σ(S) = σ² I_m - s² SS^T
//============================================================

// Covariance matrix structure for Klein sampling
struct CovarianceMatrix {
    double sigma_sq;           // σ²
    double sigma_u_sq;         // σ_u² (new scaling factor for S_mat S_mat^T)
    std::vector<double> S_mat; // bold S matrix (m × k), flattened row-major
    int m_dim;
    int k_dim;
    
    // Cholesky decomposition of Σ(S) (if needed)
    std::vector<double> chol_L; // Lower triangular, flattened
    bool chol_computed;
    
    CovarianceMatrix(int m, int k) 
            : sigma_sq(0), sigma_u_sq(0), S_mat(m*k, 0.0), 
                m_dim(m), k_dim(k), chol_L(m*m, 0.0), chol_computed(false) {}
};

// Compute Σ(S) = σ² I_m - σ_u² SS^T from secret key
CovarianceMatrix compute_covariance_from_secret(
    const SecretKey& sk, 
    double sigma_base,
    double sigma_u
);

// Klein-style sampling: sample y ~ D_{Z^m, Σ(S)}
// For each coefficient position, sample m-dimensional vector
std::vector<int32_t> sample_klein_coeff(
    const CovarianceMatrix& cov,
    int coeff_index  // which coefficient (0 to N-1)
);

// Sample full masking vector y with proper covariance
// Returns m polynomials, each of degree N
PolyIntVec sample_y_klein(const SecretKey& sk, double sigma_base, double sigma_u);

/*
// Helper: Compute ζs (component-wise multiplication by ζ)
PolyVec zeta_mul(const PolyVec& s_vec);
*/

// Helper: Compute Cholesky decomposition L where Σ = LL^T
bool compute_cholesky(CovarianceMatrix& cov);

// Helper: Sample from multivariate Gaussian using Cholesky
void sample_multivariate_gaussian(
    std::vector<int32_t>& output,
    const std::vector<double>& chol_L,
    int dim
);

} // namespace GG

#endif // GG_KLEIN_H
