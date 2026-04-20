    
#include "GG_scheme.h"
#include "GG_params.h"
#include <vector>
#include <iostream>
#include <numeric>
#include <cmath>
#include <fstream> 
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <utility>
using namespace GG;


static bool file_is_empty_or_missing(const char* path) {
    std::ifstream in(path, std::ios::in | std::ios::binary);
    if (!in.good()) return true;
    in.seekg(0, std::ios::end);
    return in.tellg() == 0;
}

static bool file_header_is(const char* path, const std::string& expected_header) {
    std::ifstream in(path);
    if (!in.good()) return false;
    std::string line;
    if (!std::getline(in, line)) return false;
    return line == expected_header;
}

static bool rename_file(const char* from, const char* to) {
    return std::rename(from, to) == 0;
}

static bool env_flag_enabled(const char* name) {
    const char* value = std::getenv(name);
    return value != nullptr && value[0] != '\0' && value[0] != '0';
}

static int env_int_or_default(const char* name, int default_value) {
    const char* value = std::getenv(name);
    if (value == nullptr || value[0] == '\0') return default_value;
    return std::atoi(value);
}

int main() {
    // One keypair, then for every (i,j) != (0,0) compute mean of Z_{i,j}/Z_{0,0} and
    // check whether the mean is inside the paper interval.
    // 자동 N_SAMPLES 계산: 논문식 N_{i,j,L}를 사용

    // 파라미터 정의
    const double sigma = GG_SIGMA_BASE;
    const double sigma_u = GG_SIGMA_U;
    const double omega= 2.326;

    // alpha_hat을 논문 수식대로 정의: (sigma_W^2 - sigma_u^2) / sigma_Z00^2
    const double sigma_W_sq = 4.0 * sigma_u * sigma_u;
    const double sigma_Z00= std::sqrt(sigma * sigma + 3.0 * sigma_u * sigma_u);
    const double alpha_hat_theory = (sigma_W_sq - sigma_u * sigma_u) / (sigma_Z00 * sigma_Z00);
    // For data-driven calibration (1번): we'll override `alpha_hat` after sampling.
    double alpha_hat = alpha_hat_theory;
    
    // 논문 조건: N_{i,j,L} * exp(-l^2/2) < 1 인 l을 자동으로 선택
    // sigma_z는 Z_{i,j}의 표준편차로 논문 수식대로 정의
    const double sigma_Zij = std::sqrt(sigma * sigma + sigma_u * sigma_u * GG_ETA * (GG_ETA + 1) * GG_N);
    double l = 2.1;
    double L = 0;
    bool found_l = false;
    const bool ratio_stdout_only = env_flag_enabled("GG_RATIO_STDOUT_ONLY");
    const int target_poly = env_int_or_default("GG_TARGET_POLY", -1);
    const int target_coeff = env_int_or_default("GG_TARGET_COEFF", -1);
    const bool has_ratio_target = target_poly >= 0 && target_coeff >= 0;
    const int max_samples_override = env_int_or_default("GG_MAX_SAMPLES", -1);//GG_MAX_SAMPLES=5000 

    std::ofstream raw_stdout;
    std::ofstream null_stdout;
    std::streambuf* saved_cout_buf = nullptr;
    if (ratio_stdout_only) {
        if (!has_ratio_target) {
            std::cerr << "[GG_RATIO_STDOUT_ONLY] set GG_TARGET_POLY and GG_TARGET_COEFF\n";
            return 1;
        }
        raw_stdout.open("/dev/stdout");
        null_stdout.open("/dev/null");
        if (null_stdout.is_open()) {
            saved_cout_buf = std::cout.rdbuf(null_stdout.rdbuf());
        }
    }

    // Generate a single keypair, and use the actual secret coefficient s_{i,j} (no forcing).
    PublicKey pk;
    SecretKey sk;
    keygen(pk, sk);

    std::cout << "[config] GG_Z_RATIO_COUNTERMEASURE="
              << (env_flag_enabled("GG_Z_RATIO_COUNTERMEASURE") ? "1" : "0")
              << " (1: sign adds (2u'+c')s' + (2u''+c'')bar{s}' to z; verify may fail)\n";

    // 개별 서명값 로깅용 CSV 파일
    const char* sig_log_path = "signature_values.csv";
    std::ofstream sig_log(sig_log_path);
    sig_log << "sample,poly,coeff,s_ij,z00,zij,ratio" << "\n";

    // Precompute secret coefficients
    std::vector<std::vector<int32_t>> s_val(GG_k, std::vector<int32_t>(GG_N, 0));
    for (int poly = 0; poly < GG_k; ++poly) {
        for (int coeff = 0; coeff < GG_N; ++coeff) {
            s_val[poly][coeff] = sk.s[poly].a[coeff];
        }
    }
    /*
    // Pick a few (poly, coeff) pairs to validate the paper's truncation event:
    //   good <=> |U - s_{i,j} * alpha_hat| <= L
    // where U is approximated by the ratio zij/z00 in this experiment.
    std::vector<std::pair<int,int>> trunc_stats_pairs;
    if (has_ratio_target && target_poly < GG_k && target_coeff < GG_N && !(target_poly == 0 && target_coeff == 0)) {
        trunc_stats_pairs.emplace_back(target_poly, target_coeff);
    } else {
        int found_pos = 0, found_neg = 0, found_zero = 0;
        for (int poly = 0; poly < GG_k; ++poly) {
            for (int coeff = 0; coeff < GG_N; ++coeff) {
                if (poly == 0 && coeff == 0) continue;
                const int32_t s_ij = s_val[poly][coeff];
                if (s_ij == 1 && !found_pos) {
                    trunc_stats_pairs.emplace_back(poly, coeff);
                    found_pos = 1;
                } else if (s_ij == -1 && !found_neg) {
                    trunc_stats_pairs.emplace_back(poly, coeff);
                    found_neg = 1;
                } else if (s_ij == 0 && !found_zero) {
                    trunc_stats_pairs.emplace_back(poly, coeff);
                    found_zero = 1;
                }
                if (found_pos + found_neg + found_zero >= 2) break;
            }
            if (found_pos + found_neg + found_zero >= 2) break;
        }
    }
    std::vector<std::vector<int>> trunc_stats_id(GG_k, std::vector<int>(GG_N, -1));
    for (int t = 0; t < (int)trunc_stats_pairs.size(); ++t) {
        const int p = trunc_stats_pairs[t].first;
        const int c = trunc_stats_pairs[t].second;
        trunc_stats_id[p][c] = t;
    }
    std::vector<long long> trunc_good_counts(trunc_stats_pairs.size(), 0);
    std::vector<long long> trunc_bad_counts(trunc_stats_pairs.size(), 0);
    std::vector<long long> trunc_checked_counts(trunc_stats_pairs.size(), 0);
    std::vector<long long> trunc_narrow_good_counts(trunc_stats_pairs.size(), 0); // |U - s*alpha| <= alpha/2*/
   
    // N_{i,j,L} 계산 
    std::vector<std::vector<int>> N_ij_L_arr(GG_k, std::vector<int>(GG_N, 0));
    while (!found_l && l < 10.0) {
        L = std::ceil(l * sigma_Zij);
        bool all_ok = true;
        for (int poly = 0; poly < GG_k; ++poly) {
            for (int coeff = 0; coeff < GG_N; ++coeff) {
                
                int32_t s_ij = s_val[poly][coeff];
                                
                // 논문 수식대로 beta_{i,j} 계산
                double num = s_ij * (sigma_W_sq - sigma_u * sigma_u);//여기서의 s_ij때문에 N_ijL을 i,j에 대해서 계산
                double denom = sigma_Z00 * sigma_Zij;
                double frac = (denom != 0.0) ? (num / denom) : 0.0;
                double beta_ij = (sigma_Zij / sigma_Z00) * std::sqrt(1.0 - frac * frac);//여기서 완성
                
                // sigma_{i,j,L} 계산
                double t = L / beta_ij;
                double atan_t = std::atan(t);
                double inner = (beta_ij * L) / atan_t - (beta_ij * beta_ij);
                double sigma_ij_L = (inner > 0.0) ? std::sqrt(inner) : 0.0;

                //N_{i,j,L} 계산
                double target = std::abs(alpha_hat) / 2.0;
                int N_ij_L = (target > 0.0) ? (int)std::ceil((omega * sigma_ij_L / target) * (omega * sigma_ij_L / target)) : 0;
                N_ij_L_arr[poly][coeff] = N_ij_L;
                
                double cond = N_ij_L * std::exp(-l * l / 2.0);
                if (cond >= 1.0) all_ok = false;
            }
        }
        if (all_ok) {
            found_l = true;
        } else {
            l += 0.1;
        }
    }

    // 최종 L, l, N_ij_L 출력
    L = std::ceil(l * sigma_Zij);
    std::cout << "[auto l] l = " << l << ", L = " << L << std::endl;
    int max_N_samples = 0;
    for (int poly = 0; poly < GG_k; ++poly) {
        for (int coeff = 0; coeff < GG_N; ++coeff) {
            if (N_ij_L_arr[poly][coeff] > max_N_samples) max_N_samples = N_ij_L_arr[poly][coeff];//최댓값을 찾아서 그만큼 서명을 생성
            std::cout << "N_ij_L[" << poly << "][" << coeff << "] = " << N_ij_L_arr[poly][coeff] << std::endl;
        }
    }
    //최댓값을 5% 증가시켜서 서명을 생성(너무 큰 서명은 reject시키기 때문에 이론적으로 필요한 서명개수보다 적게 서명을 생성하는 일을 방지)
    const int N_SAMPLES = static_cast<int>(std::ceil(max_N_samples * 1.05)); 
    //const int N_SAMPLES = 5592716;
    std::cout << "[alpha_hat] alpha_* = " << alpha_hat
              << " (sigma=" << sigma << ", sigma_u=" << sigma_u << ")\n";
    std::cout << "[N_SAMPLES] 자동 계산된 N_SAMPLES = " << N_SAMPLES << "\n";
    const int N_SAMPLES_EFF = (max_samples_override > 0) ? std::min(N_SAMPLES, max_samples_override) : N_SAMPLES;
    if (N_SAMPLES_EFF != N_SAMPLES) {
        std::cout << "[GG_MAX_SAMPLES] overriding: N_SAMPLES_EFF = " << N_SAMPLES_EFF << "\n";
    }
    const std::string corr_header = "i,j,s_ij,mean_ratio,alpha_hat,lower,upper,in_range,used,target_used";
    const char* corr_path = "correlation_data.csv";
    const char* corr_backup_path = "correlation_data_old.csv";

    // If an existing file has a different schema (e.g., in_range_rate), move it aside.
    if (!file_is_empty_or_missing(corr_path) && !file_header_is(corr_path, corr_header)) {
        // Best-effort rename; if it fails we will still append to the existing file.
        (void)rename_file(corr_path, corr_backup_path);
    }

    const bool write_corr_header = file_is_empty_or_missing(corr_path);
    std::ofstream fout(corr_path, std::ios::out | std::ios::app);
    if (write_corr_header) {
        fout << corr_header << "\n";
    }

    // Precompute interval bounds
    std::vector<std::vector<double>> lower(GG_k, std::vector<double>(GG_N, 0.0));
    std::vector<std::vector<double>> upper(GG_k, std::vector<double>(GG_N, 0.0));
    for (int poly = 0; poly < GG_k; ++poly) {
        for (int coeff = 0; coeff < GG_N; ++coeff) {
            const int32_t sv = s_val[poly][coeff];
            double lo = (double)sv * alpha_hat - alpha_hat / 2.0;
            double hi = (double)sv * alpha_hat + alpha_hat / 2.0;
            if (lo > hi) std::swap(lo, hi);
            lower[poly][coeff] = lo;
            upper[poly][coeff] = hi;
        }
    }

    std::vector<std::vector<double>> sum_ratio(GG_k, std::vector<double>(GG_N, 0.0));
    std::vector<std::vector<int>> used_arr(GG_k, std::vector<int>(GG_N, 0));

    for (int n = 0; n < N_SAMPLES_EFF; ++n) {
        Signature sig = sign(pk, sk, "test message");
        const double z00 = (double)sig.z[0].a[0];
        if (std::abs(z00) < 1e-9) continue;//바로 다음 샘플로 넘어감.
        for (int poly = 0; poly < GG_k; ++poly) {
            for (int coeff = 0; coeff < GG_N; ++coeff) {
                if (poly == 0 && coeff == 0) continue;//z00은 앞에서 계산했음
                const double zij = (double)sig.z[poly].a[coeff];//[0][1]부터 시작
                const double ratio = zij / z00;
                //디버깅용 출력
                if (ratio_stdout_only && poly == target_poly && coeff == target_coeff) {
                    raw_stdout << ratio << "\n";
                    raw_stdout.flush();
                }
                // 디버깅용; 개별 서명값 로깅(signature_values.csv 파일에 저장)
                sig_log << n << "," << poly << "," << coeff << "," << s_val[poly][coeff] << "," << z00 << "," << zij << "," << ratio << "\n";
                
                //mean_ratio 계산
                sum_ratio[poly][coeff] += ratio;
                //실제로 몇번 샘플링했는지 카운트
                used_arr[poly][coeff]++;

                /*//실제로 잘린 꼬리쪽에서 튄 샘플이 있는지 확인.
                // L-based truncation event check for selected pairs only.
                const int tid = trunc_stats_id[poly][coeff];
                if (tid >= 0) {
                    trunc_checked_counts[tid]++;
                    const double center = (double)s_val[poly][coeff] * alpha_hat_theory; // alpha_{i,j} (theory)
                    if (std::abs(ratio - center) <= L) {
                        trunc_good_counts[tid]++;
                    } else {
                        trunc_bad_counts[tid]++;
                    }
                    if (std::abs(ratio - center) <= (alpha_hat_theory / 2.0)) {
                        trunc_narrow_good_counts[tid]++;
                    }
                }*/
            }
        }
    }
    // 샘플링 루프 끝난 후 파일 닫기
    sig_log.close();

    // ============================================================
    // (1) Data-driven alpha calibration
    // Estimate alpha from empirical means: alpha_emp ≈ mean_ratio / s_ij
    // Then rebuild [lower, upper] using the empirical alpha.
    // ============================================================
    const int alpha_calibration_mode = env_int_or_default("GG_ALPHA_CALIBRATION_MODE", 1); // 1=global median, 0=disable
    if (alpha_calibration_mode != 0) {
        std::vector<double> alpha_emp;
        alpha_emp.reserve((size_t)GG_k * (size_t)GG_N);
        for (int poly = 0; poly < GG_k; ++poly) {
            for (int coeff = 0; coeff < GG_N; ++coeff) {
                if (poly == 0 && coeff == 0) continue;
                const int32_t s_ij = s_val[poly][coeff];
                const int used = used_arr[poly][coeff];
                if (s_ij == 0 || used <= 0) continue;
                const double mean = sum_ratio[poly][coeff] / (double)used;
                const double a = mean / (double)s_ij; // should be positive
                if (a > 0.0) alpha_emp.push_back(a);
            }
        }

        if (!alpha_emp.empty()) {
            std::sort(alpha_emp.begin(), alpha_emp.end());
            const double alpha_median = alpha_emp[alpha_emp.size() / 2];
            alpha_hat = alpha_median;

            // Rebuild interval bounds with calibrated alpha.
            for (int poly = 0; poly < GG_k; ++poly) {
                for (int coeff = 0; coeff < GG_N; ++coeff) {
                    const int32_t sv = s_val[poly][coeff];
                    double lo = (double)sv * alpha_hat - alpha_hat / 2.0;
                    double hi = (double)sv * alpha_hat + alpha_hat / 2.0;
                    if (lo > hi) std::swap(lo, hi);
                    lower[poly][coeff] = lo;
                    upper[poly][coeff] = hi;
                }
            }

            std::cout << "\n[Alpha calibration] alpha_hat_theory=" << alpha_hat_theory
                      << ", alpha_hat_empirical(median)=" << alpha_hat
                      << ", count=" << alpha_emp.size() << "\n";
        } else {
            std::cout << "\n[Alpha calibration] skipped (no valid alpha_emp candidates). Using alpha_hat_theory=" << alpha_hat_theory << "\n";
            alpha_hat = alpha_hat_theory;
        }
    }
    

    std::cout << "\n[N_SAMPLES vs 실제 used 샘플수]" << std::endl;
    for (int poly = 0; poly < GG_k; ++poly) {
        for (int coeff = 0; coeff < GG_N; ++coeff) {
            if (poly == 0 && coeff == 0) continue;

            int used = used_arr[poly][coeff];//분모에 해당
            const int N_target = N_SAMPLES_EFF;
            if (used > N_target) used = N_target;// 중복 안전장치.
            double mean = 0.0;
            if (used > 0) {
                mean = sum_ratio[poly][coeff] / (double)used;
            }
            // User-defined: in_range = 0 if mean_ratio is inside [lower, upper], else 1.
            const int in_range = (mean >= lower[poly][coeff] && mean <= upper[poly][coeff]) ? 0 : 1;

            fout << poly << "," << coeff << "," << s_val[poly][coeff]
                 << "," << mean
                 << "," << alpha_hat
                 << "," << lower[poly][coeff]
                 << "," << upper[poly][coeff]
                 << "," << in_range
                 << "," << used
                 << "," << N_ij_L_arr[poly][coeff]
                 << "\n";
            fout.flush();
        }
    }

    /*// Print truncation-event empirical rates for the chosen pairs.
    if (!trunc_stats_pairs.empty()) {
        std::cout << "\n[Truncation event check] good <=> |zij/z00 - s_ij*alpha_hat| <= L\n";
        std::cout << "alpha_hat_theory=" << alpha_hat_theory << ", L=" << L << "\n";
        for (int t = 0; t < (int)trunc_stats_pairs.size(); ++t) {
            const int poly = trunc_stats_pairs[t].first;
            const int coeff = trunc_stats_pairs[t].second;
            const int32_t s_ij = s_val[poly][coeff];
            const long long good = trunc_good_counts[t];
            const long long bad = trunc_bad_counts[t];
            const long long checked = trunc_checked_counts[t];
            const int used = used_arr[poly][coeff];
            const double good_rate = (checked > 0) ? (double)good / (double)checked : 0.0;
            const long long narrow_good = trunc_narrow_good_counts[t];
            const double narrow_good_rate = (checked > 0) ? (double)narrow_good / (double)checked : 0.0;
            double mean = 0.0;
            if (used > 0) mean = sum_ratio[poly][coeff] / (double)used;
            const int in_range = (mean >= lower[poly][coeff] && mean <= upper[poly][coeff]) ? 0 : 1;
            double alpha_emp = 0.0;
            bool has_alpha_emp = (s_ij != 0);
            if (has_alpha_emp) alpha_emp = mean / (double)s_ij;
            std::cout << "pair(poly=" << poly << ", coeff=" << coeff << "), s_ij=" << s_ij
                      << ", used=" << used
                      << ", checked=" << checked
                      << ", good=" << good << ", bad=" << bad
                      << ", good_rate=" << good_rate
                      << ", narrow_good=" << narrow_good << ", narrow_good_rate=" << narrow_good_rate
                      << ", mean_ratio=" << mean
                      << ", alpha_empirical=" << (has_alpha_emp ? alpha_emp : 0.0)
                      << ", in_range=" << in_range
                      << "\n";
        }
    }*/

    //std::cout << "Done. used=" << used << " (skipped " << (N_SAMPLES - used) << " where z00==0)." << std::endl;
    fout.close();
    std::cout << "\n결과가 correlation_data.csv 파일에 저장되었습니다." << std::endl;
    if (saved_cout_buf != nullptr) {
        std::cout.rdbuf(saved_cout_buf);
    }
    return 0;
}