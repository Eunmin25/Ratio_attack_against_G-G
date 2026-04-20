#ifndef GG_PARAMS_H
#define GG_PARAMS_H

#include <cstdint>
#include <cmath>

//============================================================
// G+G Signature Scheme Parameters (Toy Version)
//============================================================

// Ring dimension (toy; real would be 512, 1024, ...)
static const int GG_N = 16;

// Module dimensions
static const int GG_k = 7;  // dimension k
static const int GG_m = 3;  // dimension m (must satisfy k > m+1)

// Modulus
static const int32_t GG_q = 12289;  // prime modulus

// Secret key bound: χ_η = U({y ∈ R | ||y||_∞ ≤ η})
static const int32_t GG_ETA = 1;  // bound for uniform secret key distribution

// Gaussian widths for covariance matrix and u
static const double GG_SIGMA_BASE = 14;  // base sigma for Σ(S)
//static const double GG_SIGMA_FACTOR = 2.5;  // multiplier for Σ(s) width <- generic G+G에서는 X
static const double GG_SIGMA_U = 1;      // sigma_u for u ~ D_{R, sigma_u^2 I_n}

// G+G_Paper Sec. 4.2: γ = 1.01·√(n k)·σ + √(n m)·(1 + α/4)  (n=GG_N, k=GG_k, m=GG_m, σ=GG_SIGMA_BASE)
static const double GG_BZ_ALPHA = 256.0;
static const double GG_BZ_BOUND =
    1.01 * std::sqrt((double)GG_N * (double)GG_k) * GG_SIGMA_BASE
    + std::sqrt((double)GG_N * (double)GG_m) * (1.0 + GG_BZ_ALPHA / 4.0);

// Keygen: max attempts when rejection sampling is on (env GG_KEYGEN_MAX_RETRIES overrides)
static const int GG_KEYGEN_MAX_RETRIES = 10000;

// Target squared norm for secret key (when using Gaussian - legacy)
//static const double GG_TARGET_SQ_NORM_S = 20.0;

#endif // GG_PARAMS_H
