#include "GG_sampler.h"
#include <cmath>

namespace GG {

// IBE의 상수들 (Sampling.cc에서 사용)
extern "C" {
    // Sampling.cc에서 정의된 함수들 사용
    // Sample3, Sample4는 Sampling.h에 선언되어 있음
}

// sigma_1은 Sampling.cc에서 사용하는 상수
static const double SIGMA_1 = 1.4142135623730951;  // sqrt(2)

DGState init_sampler(double sigma) {
    DGState state;
    state.sigma = sigma;
    state.k = (unsigned long)std::ceil(sigma / SIGMA_1);
    state.coeff = 1.0 / (2 * sigma * sigma) 
                - 1.0 / (2 * state.k * state.k * SIGMA_1 * SIGMA_1);
    return state;
}

int32_t sample_gaussian(double sigma) {
    // IBE의 Sample3 사용
    return Sample3(sigma);
}

int32_t sample_gaussian_centered(double center, double sigma) {
    // IBE의 Sample4 사용
    return Sample4(center, sigma);
}

int32_t sample_gaussian_stateful(const DGState& state) {
    // Sample2는 Sampling.cc에 구현되어 있음
    signed int x;
    while(1) {
        x = Sample2(state.k);
        double alea  = ((double)rand()) / RAND_MAX;
        double borne = std::exp(-x * x * state.coeff);
        if(alea < borne) return x;
    }
}

Poly gaussian_poly(double sigma) {
    Poly p;
    for (int i = 0; i < GG_N; ++i) {
        int x = sample_gaussian(sigma);
        p.a[i] = mod_q(x);
    }
    return p;
}

Poly gaussian_poly_stateful(const DGState& state) {
    Poly p;
    for (int i = 0; i < GG_N; ++i) {
        int x = sample_gaussian_stateful(state);
        p.a[i] = mod_q(x);
    }
    return p;
}

PolyVec gaussian_vector(int dim, double sigma) {
    PolyVec v(dim);
    for (int i = 0; i < dim; ++i) {
        v[i] = gaussian_poly(sigma);
    }
    return v;
}

/*// Sample polynomial with uniform coefficients in [-eta, eta]
Poly uniform_poly_bounded(int32_t eta) {
    Poly p;
    std::uniform_int_distribution<int32_t> dist(-eta, eta);
    
    for (int i = 0; i < GG_N; ++i) {
        int32_t coeff = dist(rng);
        p.a[i] = mod_q(coeff);
    }
    return p;
}

// Sample vector of polynomials with uniform bounded coefficients
PolyVec uniform_vector_bounded(int dim, int32_t eta) {
    PolyVec v(dim);
    for (int i = 0; i < dim; ++i) {
        v[i] = uniform_poly_bounded(eta);
    }
    return v;
}*/

// Sample integer polynomial with uniform coefficients in [-eta, eta]
PolyInt uniform_polyint_bounded(int32_t eta) {
    PolyInt p = PolyInt::zero();
    std::uniform_int_distribution<int32_t> dist(-eta, eta);
    for (int i = 0; i < GG_N; ++i) {
        p.a[i] = (int64_t)dist(rng);
    }
    return p;
}

// Sample vector of integer polynomials with uniform bounded coefficients
PolyIntVec uniform_vectorint_bounded(int dim, int32_t eta) {
    PolyIntVec v(dim);
    for (int i = 0; i < dim; ++i) {
        v[i] = uniform_polyint_bounded(eta);
    }
    return v;
}

} // namespace GG
