#ifndef GG_POLY_H
#define GG_POLY_H

#include "GG_params.h"
#include <array>
#include <vector>
#include <cstdint>
#include <random>

//============================================================
// Polynomial arithmetic over Z_q[x]/(x^N + 1)
//============================================================

namespace GG {



// Modular reduction helper, q = always odd, prime
inline int32_t mod_q(int64_t x) {
    int32_t r = (int32_t)(x % GG_q);
    if (r < 0) r += GG_q;                  // r in [0, q)
    int32_t half = (int32_t)((GG_q - 1) / 2);
    if (r > half) r -= (int32_t)GG_q;      // r in [-(q-1)/2, (q-1)/2]
    return r;
}

// Modular reduction helper for Z_{2q}
inline int32_t mod_2q(int64_t x) {
    const int64_t Q2 = (int64_t)GG_q * 2;
    int32_t r = (int32_t)(x % Q2);
    if (r < 0) r += (int32_t)Q2;
    const int32_t half = (int32_t)GG_q;   // q
    if (r > half) r -= (int32_t)Q2;       // r in (-q, q]
    return r;
}


// Integer polynomial (no modular reduction).
struct PolyInt {
    std::array<int64_t, GG_N> a{};

    static PolyInt zero();
    PolyInt& operator+=(const PolyInt& other);
    PolyInt& mul_scalar(int64_t c);

    friend PolyInt operator+(PolyInt lhs, const PolyInt& rhs);
    friend PolyInt operator*(const PolyInt& f, const PolyInt& g);
};

using PolyIntVec = std::vector<PolyInt>;

struct MatrixInt {
    int rows, cols;
    std::vector<PolyInt> data;

    MatrixInt(int r=0, int c=0);
    PolyInt& at(int r, int c);
    const PolyInt& at(int r, int c) const;
};

PolyIntVec mat_vec_mul(const MatrixInt& A, const PolyIntVec& x);


// ---------------- mod 2q polynomial arithmetic ----------------
// Separate type from Poly to avoid mixing mod-q and mod-2q invariants.
struct Poly2Q {
    std::array<int32_t, GG_N> a{}; // coefficients reduced by mod_2q

    Poly2Q(bool zero = true);
    static Poly2Q zero();

    Poly2Q& operator+=(const Poly2Q& other);
    Poly2Q& operator-=(const Poly2Q& other);
    Poly2Q& mul_scalar(int32_t c);

    friend Poly2Q operator+(Poly2Q lhs, const Poly2Q& rhs);
    friend Poly2Q operator-(Poly2Q lhs, const Poly2Q& rhs);
    friend Poly2Q operator*(const Poly2Q& f, const Poly2Q& g);
};

using Poly2QVec = std::vector<Poly2Q>;

// Multiply a Poly2Q by an integer polynomial and reduce mod 2q.
Poly2Q poly2q_mul_mod2q(const Poly2Q& f, const PolyInt& g);

// Lift a mod-2q polynomial to an integer polynomial using centered representatives in (-q, q].
PolyInt polyint_from_poly2q_centered(const Poly2Q& p);

// Integer (non-modular) negacyclic convolution between PolyInt and Poly2Q.
// (Poly2Q coefficients are used as their stored representatives.)
PolyInt operator*(const PolyInt& f, const Poly2Q& g);
PolyInt operator*(const Poly2Q& f, const PolyInt& g);

struct Matrix2Q {
    int rows, cols;
    std::vector<Poly2Q> data;

    Matrix2Q(int r=0, int c=0);
    Poly2Q& at(int r, int c);
    const Poly2Q& at(int r, int c) const;
};


// Matrix-vector multiplication mod 2q (two variants used in the codebase).
Poly2QVec mat_vec_mul_mod2q(const Matrix2Q& A, const Poly2QVec& x);
Poly2QVec mat_vec_mul_mod2q(const Matrix2Q& A, const  PolyIntVec& x);


//다항식 한개
struct Poly {
    std::array<int32_t, GG_N> a{}; // coefficients mod q
    //기본값은 모든 계수를 0으로 초기화
    Poly(bool zero = true);
    //계수를 균등분포에서 샘플링한 다항식 생성
    static Poly random_uniform();
    //모든 계수가 0인 다항식 생성
    static Poly zero();
    
    //다항식 덧셈/뺄셈/스칼라 c 곱셈
    Poly& operator+=(const Poly& other); //a=a+b
    Poly& operator-=(const Poly& other);
    Poly& mul_scalar(int32_t c);
    
    //두 피연산자 모두 자유롭게 받을 수 있음
    friend Poly operator+(Poly lhs, const Poly& rhs);//c=a+b
    friend Poly operator-(Poly lhs, const Poly& rhs);
    friend Poly operator*(const Poly& f, const Poly& g);
};
//여러개의 다항식으로 이루어진 벡터, std::vector는 동적 배열 컨테이너
using PolyVec = std::vector<Poly>;


struct MatrixR {
    int rows, cols;
    std::vector<Poly> data;
    //rxc 크기의 2차원 행렬로 각 원소가 다항삭
    MatrixR(int r=0, int c=0);
    //접근할 때 A.at(i,j)로 poly에 접근
    Poly& at(int r, int c);
    const Poly& at(int r, int c) const;
    //균등분포에에서 
    static MatrixR random_uniform(int r, int c);
};

// Matrix-vector multiplication 함수 선언으로 반환 타입이 다항식 벡터(PolyVec)임
PolyVec mat_vec_mul(const MatrixR& A, const PolyVec& x);

// ---------- Conversions between domains ----------
// Convert integer polynomials to mod-q polynomials (coefficient-wise mod_q).
Poly poly_from_int_modq(const PolyInt& p);
PolyVec polyvec_from_intvec_modq(const PolyIntVec& v);

// Convert integer polynomials to mod-2q polynomials (coefficient-wise mod_2q).
Poly2Q poly2q_from_int_mod2q(const PolyInt& p);
Poly2QVec poly2qvec_from_intvec_mod2q(const PolyIntVec& v);

// Convert mod-2q polynomials to mod-q polynomials (coefficient-wise mod_q).
Poly poly_from_poly2q_modq(const Poly2Q& p);
PolyVec polyvec_from_poly2qvec_modq(const Poly2QVec& v);

// Convert mod-q polynomials to mod-2q polynomials (coefficient-wise mod_2q).
Poly2Q poly2q_from_poly_mod2q(const Poly& p);
Poly2QVec poly2qvec_from_polyvec_mod2q(const PolyVec& v);

// Random number generator
extern std::mt19937_64 rng;

bool polyint_equal(const PolyInt& a, const PolyInt& b);

} // namespace GG

#endif // GG_POLY_H
