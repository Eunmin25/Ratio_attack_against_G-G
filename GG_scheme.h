#include <string>
#ifndef GG_SCHEME_H
#define GG_SCHEME_H

#include "GG_params.h"
#include "GG_poly.h"
#include <string>

namespace GG {
// PolyVec 디버깅 출력 함수 선언
void print_polyvec(const PolyVec& v, const std::string& name);
void print_polyvec_full(const PolyVec& v, const std::string& name);
void print_polyintvec(const PolyIntVec& v, const std::string& name);
void print_polyintvec_full(const PolyIntVec& v, const std::string& name);

//============================================================
// G+G Signature Scheme Structures
//============================================================

struct SecretKey {
    Poly2QVec s;  // length k
};

struct PublicKey {
    Matrix2Q A;  // m x k
    //Matrix2Q A2q; // m x k (paper: lives mod 2q)
};

struct Signature {
    PolyIntVec z;  // length k (paper: z \in R^k)
    PolyInt c;     // challenge polynomial
};

// Paper-style signature where z is kept as integers (no modular reduction).
struct SignatureInt {
    std::vector<PolyInt> z; // length k
    Poly c;
};

//============================================================
// G+G Signature Scheme Functions
//============================================================
// Note on Gaussian Sampling:
// - In full G+G paper, y is sampled from D_{R^m, Σ(s)}
//   where Σ(s) is derived from the spectral norm σ₁(rot(ζs))
// - This toy implementation uses a simplified approximation:
//   y ~ D_{σ_y} with σ_y proportional to σ_s
//============================================================

// Key generation
void keygen(PublicKey& pk, SecretKey& sk);

// Sign a message
Signature sign(const PublicKey& pk, const SecretKey& sk, const std::string& message);

// Paper-style sign: v := A*y (mod 2q), c := H(v,m), z := y + (2u+c)s over integers.
SignatureInt sign_paper(const PublicKey& pk, const SecretKey& sk, const std::string& message);

// Verify a signature (toy version, may not be fully implemented)
bool verify(const PublicKey& pk, const std::string& message, const Signature& sig);

// Hash function: H(w, mu) -> challenge polynomial
PolyInt hash_to_challenge(const Poly2QVec& w, const std::string& message);

//============================================================
// Utility Functions for Printing
//============================================================

void print_poly(const Poly& p, const std::string& name = "Poly");
void print_polyint(const PolyInt& p, const std::string& name = "PolyInt");
void print_polyvec(const PolyVec& v, const std::string& name = "PolyVec");
void print_polyvec_full(const PolyVec& v, const std::string& name = "PolyVec(full)");
void print_polyintvec(const PolyIntVec& v, const std::string& name = "PolyIntVec");
void print_polyintvec_full(const PolyIntVec& v, const std::string& name = "PolyIntVec(full)");
void print_matrix(const MatrixR& M, const std::string& name = "Matrix");
void print_matrix(const Matrix2Q& M, const std::string& name = "Matrix2Q");
void print_public_key(const PublicKey& pk);
void print_secret_key(const SecretKey& sk);
void print_signature(const Signature& sig);
bool verify(const PublicKey& pk, const std::string& message, const Signature& sig);
} // namespace GG

#endif // GG_SCHEME_H
