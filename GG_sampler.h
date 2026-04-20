#ifndef GG_SAMPLER_H
#define GG_SAMPLER_H

#include "GG_params.h"
#include "GG_poly.h"
#include "Sampling.h"  // IBE의 Sample0~Sample4 사용

namespace GG {

//============================================================
// Discrete Gaussian Sampler Wrapper
// Uses IBE project's Sampling.h functions
//============================================================

// State for efficient repeated sampling with fixed sigma
struct DGState {
    double sigma;
    unsigned long k;
    double coeff;
};

// Initialize sampler state for given sigma
DGState init_sampler(double sigma);

// Sample from discrete Gaussian D_{sigma} (center = 0)
int32_t sample_gaussian(double sigma);

// Sample from discrete Gaussian D_{c,sigma} (center = c)
int32_t sample_gaussian_centered(double center, double sigma);

// Sample from discrete Gaussian using pre-initialized state
int32_t sample_gaussian_stateful(const DGState& state);

// Sample a Gaussian polynomial (all coefficients independent)
Poly gaussian_poly(double sigma);

// Sample a Gaussian polynomial with pre-initialized state
Poly gaussian_poly_stateful(const DGState& state);

// Sample a Gaussian vector
PolyVec gaussian_vector(int dim, double sigma);

// Sample a polynomial with uniform coefficients in [-eta, eta]
Poly uniform_poly_bounded(int32_t eta);

// Sample a vector of polynomials with uniform bounded coefficients
PolyVec uniform_vector_bounded(int dim, int32_t eta);

// Sample an *integer* polynomial with uniform coefficients in [-eta, eta]
PolyInt uniform_polyint_bounded(int32_t eta);

// Sample a vector of integer polynomials with uniform bounded coefficients
PolyIntVec uniform_vectorint_bounded(int dim, int32_t eta);

} // namespace GG

#endif // GG_SAMPLER_H
