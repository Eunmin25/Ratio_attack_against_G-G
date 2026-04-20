#include "GG_poly.h"
#include "GG_scheme.h"
#include <random>

namespace GG {

static std::random_device rd;
std::mt19937_64 rng(rd());  // Non-static for external use

// ---------- PolyInt implementation ----------

PolyInt PolyInt::zero() {
    PolyInt p;
    for (int i = 0; i < GG_N; ++i) p.a[i] = 0;
    return p;
}

PolyInt& PolyInt::operator+=(const PolyInt& other) {
    for (int i = 0; i < GG_N; ++i) a[i] += other.a[i];
    return *this;
}

PolyInt& PolyInt::mul_scalar(int64_t c) {
    for (int i = 0; i < GG_N; ++i) a[i] *= c;
    return *this;
}

PolyInt operator+(PolyInt lhs, const PolyInt& rhs) {
    lhs += rhs;
    return lhs;
}

PolyInt operator*(const PolyInt& f, const PolyInt& g) {
    PolyInt res = PolyInt::zero();
    for (int i = 0; i < GG_N; ++i) {
        for (int j = 0; j < GG_N; ++j) {
            int idx = i + j;
            int64_t sign = 1;
            if (idx >= GG_N) {
                idx -= GG_N;
                sign = -1;
            }
            __int128 prod = (__int128)f.a[i] * (__int128)g.a[j] * (__int128)sign;
            res.a[idx] += (int64_t)prod;
        }
    }
    return res;
}

PolyInt polyint_from_poly2q_centered(const Poly2Q& p) {
    PolyInt out = PolyInt::zero();
    for (int i = 0; i < GG_N; ++i) {
        // Normalize to centered representative even if invariants drift.
        out.a[i] = (int64_t)mod_2q((int64_t)p.a[i]);
    }
    return out;
}

PolyInt operator*(const PolyInt& f, const Poly2Q& g) {
    return f * polyint_from_poly2q_centered(g);
}

PolyInt operator*(const Poly2Q& f, const PolyInt& g) {
    return g * f;
}

// ---------- MatrixInt implementation ----------

MatrixInt::MatrixInt(int r, int c) : rows(r), cols(c), data((size_t)r * (size_t)c) {}

PolyInt& MatrixInt::at(int r, int c) {
    return data[(size_t)r * (size_t)cols + (size_t)c];
}

const PolyInt& MatrixInt::at(int r, int c) const {
    return data[(size_t)r * (size_t)cols + (size_t)c];
}

PolyIntVec mat_vec_mul(const MatrixInt& A, const PolyIntVec& x) {
    PolyIntVec res((size_t)A.rows, PolyInt::zero());
    for (int r = 0; r < A.rows; ++r) {
        PolyInt acc = PolyInt::zero();
        for (int c = 0; c < A.cols; ++c) {
            acc += A.at(r, c) * x[(size_t)c];
        }
        res[(size_t)r] = acc;
    }
    return res;
}

Poly poly_from_int_modq(const PolyInt& p) {
    Poly out = Poly::zero();
    for (int i = 0; i < GG_N; ++i) out.a[i] = mod_q(p.a[i]);
    return out;
}

PolyVec polyvec_from_intvec_modq(const PolyIntVec& v) {
    PolyVec out;
    out.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = poly_from_int_modq(v[i]);
    return out;
}

Poly2Q poly2q_from_int_mod2q(const PolyInt& p) {
    Poly2Q out = Poly2Q::zero();
    for (int i = 0; i < GG_N; ++i) out.a[i] = mod_2q(p.a[i]);
    return out;
}

Poly2QVec poly2qvec_from_intvec_mod2q(const PolyIntVec& v) {
    Poly2QVec out;
    out.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = poly2q_from_int_mod2q(v[i]);
    return out;
}

Poly poly_from_poly2q_modq(const Poly2Q& p) {
    Poly out = Poly::zero();
    for (int i = 0; i < GG_N; ++i) out.a[i] = mod_q((int64_t)p.a[i]);
    return out;
}

PolyVec polyvec_from_poly2qvec_modq(const Poly2QVec& v) {
    PolyVec out;
    out.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = poly_from_poly2q_modq(v[i]);
    return out;
}

Poly2Q poly2q_from_poly_mod2q(const Poly& p) {
    Poly2Q out = Poly2Q::zero();
    for (int i = 0; i < GG_N; ++i) out.a[i] = mod_2q((int64_t)p.a[i]);
    return out;
}

Poly2QVec poly2qvec_from_polyvec_mod2q(const PolyVec& v) {
    Poly2QVec out;
    out.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i) out[i] = poly2q_from_poly_mod2q(v[i]);
    return out;
}


// ---------- Poly implementation ----------

Poly::Poly(bool zero) {
    if (zero) {
        for (int i = 0; i < GG_N; ++i) a[i] = 0;
    }
}
//uniform_poly_bounded or uniform_vector_bounded와 mod q 연산과 연관 , A0 행렬 생성에 사용
Poly Poly::random_uniform() {
    Poly p;
    std::uniform_int_distribution<int32_t> dist(-((GG_q-1)/2), (GG_q-1)/2);
    for (int i = 0; i < GG_N; ++i) p.a[i] = mod_q(dist(rng));
    return p;
}

Poly Poly::zero() { 
    return Poly(true); 
}
//자기 자신에 other를 더하거나 빼거나 
Poly& Poly::operator+=(const Poly& other) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_q((int64_t)a[i] + other.a[i]);
    }
    return *this;
}

Poly& Poly::operator-=(const Poly& other) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_q((int64_t)a[i] - other.a[i]);
    }
    return *this;
}
// 자가 자신에서 모든 계수에 c를 곱함.
Poly& Poly::mul_scalar(int32_t c) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_q((int64_t)a[i] * c);
    }
    return *this;
}

//lhs에 rhs를 더하거나 빼거나 곱한 결과를 새 Poly로 반환
Poly operator+(Poly lhs, const Poly& rhs) {
    //내부적으로는 += 연산자 사용
    lhs += rhs; 
    return lhs;
}

Poly operator-(Poly lhs, const Poly& rhs) {
    lhs -= rhs; 
    return lhs;
}

Poly operator*(const Poly& f, const Poly& g) {
    Poly res;
    for (int i = 0; i < GG_N; ++i) res.a[i] = 0;
    for (int i = 0; i < GG_N; ++i) {
        for (int j = 0; j < GG_N; ++j) {
            int idx = i + j;
            //부호 선언 및 초기화
            int32_t sign = 1;
            // x^GG_N = -1 적용
            if (idx >= GG_N) { 
                idx -= GG_N; 
                sign = -1; 
            }
            int64_t val = (int64_t)f.a[i] * g.a[j] * sign;
            res.a[idx] = mod_q((int64_t)res.a[idx] + val);
        }
    }
    return res;
}
// ---------- MatrixR implementation ----------

MatrixR::MatrixR(int r, int c) : rows(r), cols(c), data(r*c) {}

// (r,c) 위치의 Poly 원소에 접근
Poly& MatrixR::at(int r, int c) { 
    return data[r*cols + c]; 
}

const Poly& MatrixR::at(int r, int c) const { 
    return data[r*cols + c]; 
}

MatrixR MatrixR::random_uniform(int r, int c) {
    MatrixR M(r, c);
    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            M.at(i,j) = Poly::random_uniform();
    return M;
}


// ---------- Poly2Q implementation ----------

Poly2Q::Poly2Q(bool zero) {
    if (zero) {
        for (int i = 0; i < GG_N; ++i) a[i] = 0;
    }
}

Poly2Q Poly2Q::zero() {
    return Poly2Q(true);
}

Poly2Q& Poly2Q::operator+=(const Poly2Q& other) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_2q((int64_t)a[i] + other.a[i]);
    }
    return *this;
}

Poly2Q& Poly2Q::operator-=(const Poly2Q& other) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_2q((int64_t)a[i] - other.a[i]);
    }
    return *this;
}

Poly2Q& Poly2Q::mul_scalar(int32_t c) {
    for (int i = 0; i < GG_N; ++i) {
        a[i] = mod_2q((int64_t)a[i] * (int64_t)c);
    }
    return *this;
}

Poly2Q operator+(Poly2Q lhs, const Poly2Q& rhs) {
    lhs += rhs;
    return lhs;
}

Poly2Q operator-(Poly2Q lhs, const Poly2Q& rhs) {
    lhs -= rhs;
    return lhs;
}

Poly2Q operator*(const Poly2Q& f, const Poly2Q& g) {
    Poly2Q res = Poly2Q::zero();
    for (int i = 0; i < GG_N; ++i) {
        for (int j = 0; j < GG_N; ++j) {
            int idx = i + j;
            int64_t sign = 1;
            if (idx >= GG_N) {
                idx -= GG_N;
                sign = -1;
            }
            __int128 prod = (__int128)f.a[i] * (__int128)g.a[j] * (__int128)sign;
            int64_t acc = (int64_t)res.a[idx] + (int64_t)prod;
            res.a[idx] = mod_2q(acc);
        }
    }
    return res;
}



// ---------- Matrix2Q implementation ----------

Matrix2Q::Matrix2Q(int r, int c) : rows(r), cols(c), data((size_t)r * (size_t)c) {}

Poly2Q& Matrix2Q::at(int r, int c) {
    return data[(size_t)r * (size_t)cols + (size_t)c];
}

const Poly2Q& Matrix2Q::at(int r, int c) const {
    return data[(size_t)r * (size_t)cols + (size_t)c];
}

Poly2Q poly2q_mul_mod2q(const Poly2Q& f, const PolyInt& g) {
    Poly2Q res = Poly2Q::zero();
    for (int i = 0; i < GG_N; ++i) {
        for (int j = 0; j < GG_N; ++j) {
            int idx = i + j;
            int64_t sign = 1;
            if (idx >= GG_N) {
                idx -= GG_N;
                sign = -1;
            }
            __int128 prod = (__int128)f.a[i] * (__int128)g.a[j] * (__int128)sign;
            int64_t acc = (int64_t)res.a[idx] + (int64_t)prod;
            res.a[idx] = mod_2q(acc);
        }
    }
    return res;
}

Poly2QVec mat_vec_mul_mod2q(const Matrix2Q& A, const Poly2QVec& x) {
    Poly2QVec res((size_t)A.rows, Poly2Q::zero());
    for (int r = 0; r < A.rows; ++r) {
        Poly2Q acc = Poly2Q::zero();
        for (int c = 0; c < A.cols; ++c) {
            acc += A.at(r, c) * x[(size_t)c];
        }
        res[(size_t)r] = acc;
    }
    return res;
}

Poly2QVec mat_vec_mul_mod2q(const Matrix2Q& A, const PolyIntVec& x) {
    Poly2QVec res((size_t)A.rows, Poly2Q::zero());
    for (int r = 0; r < A.rows; ++r) {
        Poly2Q acc = Poly2Q::zero();
        for (int c = 0; c < A.cols; ++c) {
            Poly2Q term = poly2q_mul_mod2q(A.at(r, c), x[(size_t)c]);
            acc += term;
        }
        res[(size_t)r] = acc;
    }
    return res;
}

bool polyint_equal(const PolyInt& a, const PolyInt& b) {
    for (int i = 0; i < GG_N; ++i) {
        if (a.a[i] != b.a[i]) return false;
    }
    return true;
}

// ---------- Matrix-vector multiplication ----------

PolyVec mat_vec_mul(const MatrixR& A, const PolyVec& x) {
    PolyVec res(A.rows, Poly::zero());
    for (int i = 0; i < A.rows; ++i) {
        Poly acc = Poly::zero();
        for (int j = 0; j < A.cols; ++j) {
            acc += A.at(i,j) * x[j];
        }
        res[i] = acc;
    }
    return res;
}

} // namespace GG
