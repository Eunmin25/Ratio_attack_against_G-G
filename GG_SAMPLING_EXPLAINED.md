sk.s[0] = polynomial "1"  // [1, 0, 0, ..., 0]
s[0] = [1, 0, 0, 0, ...]           // polynomial "1"
s[1] = [1, 0, 0, 0, ...]           // small Gaussian coeffs
s[3] = [0, -1, 0, -1, ...]
s[4] = [-1, -1, 1, -1, ...]
double sigma_base = GG_SIGMA_S * GG_SIGMA_FACTOR;  // 1.2 × 2.5 = 3.0
y[0] = [3, -2, 1, 4, ...]    // 각 계수가 Σ(S)의 
y[1] = [-1, 3, 2, -2, ...]   // 공분산을 가진
y[2] = [2, 1, -3, 1, ...]    // 다변량 가우시안
w = A · y (mod q)
w[i] = Σ_{j=0}^{k-1} A[i][j] * y[j % m]  (mod q, in R_q)
c = [697, 12183, 10465, 4108, 5670, 9486, 9807, 1184, ...]
z[0] = y[0] + c · s[0]
z[1] = y[1] + c · s[1]
z[2] = y[2] + c · s[2]
z[0] = [700, 12183, 10465, 4108, ...]
z[1] = [856, 10597, 4469, 6820, ...]
z[2] = [4420, 5313, 3238, 3515, ...]

# GG Signature Scheme: Sampling & Parameter Flow (최신 코드 기준)

## 1. 파라미터 요약

| 기호 | 설명 | Toy Example | 논문(실제) |
|------|------|-------------|------------|
| $n$ | 다항식 차수 | 16 | 256 |
| $q$ | 모듈러스 | 12289 | 64513 등 |
| $m, k$ | 행렬 차원 | 3, 5 | (3,4), (4,5) 등 |
| $\eta$ | 비밀키 계수 범위 | 2 | 1 |
| $\sigma$ | 가우시안 표준편차 | 3.0 | 640~728 |
| $S$ | 비밀키 노름 제한 | 8 | 79~90 |
| $\gamma$ | 서명 노름 제한 | 20 | 31972~38437 |

*Toy 예제는 구조·실험용, 실제 보안은 논문 파라미터 사용*

---

## 2. 주요 샘플링 방식

### 비밀키 $s$
- **분포:** $[-\eta, \eta]$에서 균등 (uniform)
- **코드:** `uniform_poly_bounded` (GG_sampler)
- **체크:** $\ell_2$(rot($\zeta s$)) $< S$ (리젝션)

### 마스킹 $y$
- **분포:** $D_{\mathcal{R}_k, \Sigma(s)}$ (Klein 샘플링)
- **코드:** `klein_sample` (GG_klein)
- **공분산:** $\Sigma(s) = \sigma^2 I_m - s^2 SS^T$

### 챌린지 $c$
- **분포:** $c = H(w, \mu)$ (SHA-256)

### 응답 $z$
- **계산:** $z = y + (\zeta u + c)s$
- **체크:** $\|z\| \leq \gamma$

---

## 3. 전체 흐름 요약

**KeyGen**
1. $A, A_0 \leftarrow$ uniform mod $q$
2. $s_1, s_2 \leftarrow$ uniform in $[-\eta, \eta]$
3. $b \leftarrow a + A_0 s_1 + s_2 \pmod{q}$
4. $s \leftarrow (1|s_1^T|s_2^T - b_0^T)$
5. $\sigma_1$(rot($\zeta s$)) $\geq S$면 리젝션
6. $(A, s)$ 반환

**Sign**
1. $y \leftarrow D_{\mathcal{R}_k, \Sigma(s)}$
2. $w \leftarrow Ay \bmod 2q$
3. $c \leftarrow H(w, \mu)$
4. $u \leftarrow D_{\mathcal{R}, s, -\zeta^* c/2}$
5. $z \leftarrow y + (\zeta u + c)s$
6. $(z, c)$ 반환

**Verify**
1. $w \leftarrow Az - qcj \bmod 2q$
2. $c = H(w, \mu)$, $\|z\| \leq \gamma$면 accept

---

## 4. 구현 팁

- 모든 샘플링 함수는 GG_sampler, GG_klein에 모듈화됨
- 비밀키는 항상 균등분포 (가우시안 아님)
- 파라미터는 GG_params.h에서 쉽게 변경
- toy/실제 파라미터 교체만으로 실험/보안 모두 가능

---

## 5. 실제 파라미터 적용법

1. GG_params.h에서 $n$, $q$, $m$, $k$, $\eta$, $\sigma$, $S$, $\gamma$를 논문 값으로 변경
2. 빌드 & 실행
