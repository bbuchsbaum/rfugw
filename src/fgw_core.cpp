// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <limits>
#include <numeric>
#include <queue>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

constexpr double kTiny = 1e-300;
constexpr double kExpClip = 700.0;
constexpr float kTinyF = 1e-20f;
constexpr float kExpClipF = 80.0f;
constexpr arma::uword kMatvecBlockRows = 128;
constexpr arma::uword kMatvecBlockCols = 64;
constexpr arma::uword kMatvecBlockedMinWorkD = 25000;   // ~160x160
constexpr arma::uword kMatvecBlockedMinWorkF = 160000;  // ~400x400
constexpr arma::uword kGemvMinWorkD = 120000;           // ~346x346
constexpr arma::uword kGemvMinWorkF = 25000;            // ~160x160
constexpr arma::uword kKernelOmpMinWork = 200000;
#ifdef _OPENMP
constexpr arma::uword kMatvecOmpMinWork = 200000;
#endif

inline bool can_use_inner_omp(arma::uword work, arma::uword min_work) {
#ifdef _OPENMP
  return work >= min_work && omp_get_max_threads() > 1 && omp_in_parallel() == 0;
#else
  (void)work;
  (void)min_work;
  return false;
#endif
}

inline arma::uword parse_uword_env(const char* name, arma::uword fallback) {
  const char* raw = std::getenv(name);
  if (raw == nullptr || raw[0] == '\0') {
    return fallback;
  }
  char* end = nullptr;
  const unsigned long long v = std::strtoull(raw, &end, 10);
  if (end == raw || (end != nullptr && *end != '\0')) {
    return fallback;
  }
  return static_cast<arma::uword>(v);
}

inline double parse_double_env(const char* name, double fallback) {
  const char* raw = std::getenv(name);
  if (raw == nullptr || raw[0] == '\0') {
    return fallback;
  }
  char* end = nullptr;
  const double v = std::strtod(raw, &end);
  if (end == raw || (end != nullptr && *end != '\0') || !std::isfinite(v) || v <= 0.0) {
    return fallback;
  }
  return v;
}

inline bool fugw_enable_warm_start() {
  static const bool enabled = (parse_uword_env("RFUGW_ENABLE_WARM_START", 0u) != 0u);
  return enabled;
}

inline bool fugw_enable_adaptive_inner_tol() {
  static const bool enabled = (parse_uword_env("RFUGW_ADAPTIVE_INNER_TOL", 0u) != 0u);
  return enabled;
}

inline double fugw_adaptive_inner_tol_stage1() {
  static const double v = parse_double_env("RFUGW_ADAPTIVE_INNER_TOL_STAGE1", 1e-6);
  return v;
}

inline double fugw_adaptive_inner_tol_stage2() {
  static const double v = parse_double_env("RFUGW_ADAPTIVE_INNER_TOL_STAGE2", 5e-7);
  return v;
}

inline double fugw_effective_inner_tol(double base_tol_ot, double outer_err, double outer_tol) {
  if (!fugw_enable_adaptive_inner_tol()) {
    return base_tol_ot;
  }
  if (!(base_tol_ot > 0.0) || !std::isfinite(base_tol_ot)) {
    return base_tol_ot;
  }

  const double stage1 = std::max(base_tol_ot, fugw_adaptive_inner_tol_stage1());
  const double stage2 = std::max(base_tol_ot, fugw_adaptive_inner_tol_stage2());
  if (!std::isfinite(outer_err)) {
    return stage1;
  }

  const double gate1 = std::max(100.0 * outer_tol, 1e-6);
  const double gate2 = std::max(10.0 * outer_tol, 1e-7);
  if (outer_err > gate1) {
    return stage1;
  }
  if (outer_err > gate2) {
    return stage2;
  }
  return base_tol_ot;
}

inline double fugw_sinkhorn_target_tol_abs_mult() {
  static const double v = parse_double_env("RFUGW_SINKHORN_TARGET_TOL_ABS_MULT", 100.0);
  return v;
}

inline double fugw_sinkhorn_target_tol_mult() {
  static const double v = parse_double_env("RFUGW_SINKHORN_TARGET_TOL_MULT", 52.0);
  return v;
}

inline double fugw_sinkhorn_target_tol_cap_d() {
  static const double v = parse_double_env("RFUGW_SINKHORN_TARGET_TOL_CAP_D", 1.5e-4);
  return v;
}

inline double fugw_sinkhorn_target_tol_cap_f() {
  static const double v = parse_double_env("RFUGW_SINKHORN_TARGET_TOL_CAP_F", 2e-4);
  return v;
}

inline double fugw_sinkhorn_rel_tol_mult() {
  static const double v = parse_double_env("RFUGW_SINKHORN_REL_TOL_MULT", 200.0);
  return v;
}

inline double fugw_sinkhorn_rel_tol_floor_d() {
  static const double v = parse_double_env("RFUGW_SINKHORN_REL_TOL_FLOOR_D", 4e-5);
  return v;
}

inline double fugw_sinkhorn_rel_tol_floor_f() {
  static const double v = parse_double_env("RFUGW_SINKHORN_REL_TOL_FLOOR_F", 2e-4);
  return v;
}

inline double fugw_sinkhorn_col_rel_tol_mult() {
  static const double v = parse_double_env("RFUGW_SINKHORN_COL_REL_TOL_MULT", 1000.0);
  return v;
}

inline double fugw_sinkhorn_col_rel_tol_floor_d() {
  static const double v = parse_double_env("RFUGW_SINKHORN_COL_REL_TOL_FLOOR_D", 1.5e-4);
  return v;
}

inline double fugw_sinkhorn_col_rel_tol_floor_f() {
  static const double v = parse_double_env("RFUGW_SINKHORN_COL_REL_TOL_FLOOR_F", 2e-4);
  return v;
}

inline void dgemm_nn(
    const arma::mat& A,
    const arma::mat& B,
    arma::mat& C,
    double alpha = 1.0,
    double beta = 0.0) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_cols);
  const char transN = 'N';
  C.set_size(A.n_rows, B.n_cols);
  arma::blas::gemm<double>(
    &transN, &transN,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &k,
    &beta,
    C.memptr(), &n
  );
}

inline void dgemm_nt(
    const arma::mat& A,
    const arma::mat& B,
    arma::mat& C,
    double alpha = 1.0,
    double beta = 0.0) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_rows);
  const char transN = 'N';
  const char transT = 'T';
  C.set_size(A.n_rows, B.n_rows);
  arma::blas::gemm<double>(
    &transN, &transT,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &m,
    &beta,
    C.memptr(), &n
  );
}

inline void dgemm_nn_accum(
    const arma::mat& A,
    const arma::mat& B,
    arma::mat& C,
    double alpha,
    double beta) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_cols);
  const char transN = 'N';
  if (C.n_rows != A.n_rows || C.n_cols != B.n_cols) {
    C.zeros(A.n_rows, B.n_cols);
    beta = 0.0;
  }
  arma::blas::gemm<double>(
    &transN, &transN,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &k,
    &beta,
    C.memptr(), &n
  );
}

inline void dgemm_nt_accum(
    const arma::mat& A,
    const arma::mat& B,
    arma::mat& C,
    double alpha,
    double beta) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_rows);
  const char transN = 'N';
  const char transT = 'T';
  if (C.n_rows != A.n_rows || C.n_cols != B.n_rows) {
    C.zeros(A.n_rows, B.n_rows);
    beta = 0.0;
  }
  arma::blas::gemm<double>(
    &transN, &transT,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &m,
    &beta,
    C.memptr(), &n
  );
}

inline void sgemm_nn(
    const arma::fmat& A,
    const arma::fmat& B,
    arma::fmat& C,
    float alpha = 1.0f,
    float beta = 0.0f) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_cols);
  const char transN = 'N';
  C.set_size(A.n_rows, B.n_cols);
  arma::blas::gemm<float>(
    &transN, &transN,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &k,
    &beta,
    C.memptr(), &n
  );
}

inline void sgemm_nn_accum(
    const arma::fmat& A,
    const arma::fmat& B,
    arma::fmat& C,
    float alpha,
    float beta) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_cols);
  const char transN = 'N';
  if (C.n_rows != A.n_rows || C.n_cols != B.n_cols) {
    C.zeros(A.n_rows, B.n_cols);
    beta = 0.0f;
  }
  arma::blas::gemm<float>(
    &transN, &transN,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &k,
    &beta,
    C.memptr(), &n
  );
}

inline void sgemm_nt_accum(
    const arma::fmat& A,
    const arma::fmat& B,
    arma::fmat& C,
    float alpha,
    float beta) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_rows);
  const char transN = 'N';
  const char transT = 'T';
  if (C.n_rows != A.n_rows || C.n_cols != B.n_rows) {
    C.zeros(A.n_rows, B.n_rows);
    beta = 0.0f;
  }
  arma::blas::gemm<float>(
    &transN, &transT,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &m,
    &beta,
    C.memptr(), &n
  );
}

inline void sgemm_nt(
    const arma::fmat& A,
    const arma::fmat& B,
    arma::fmat& C,
    float alpha = 1.0f,
    float beta = 0.0f) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int k = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int m = static_cast<arma::blas_int>(B.n_rows);
  const char transN = 'N';
  const char transT = 'T';
  C.set_size(A.n_rows, B.n_rows);
  arma::blas::gemm<float>(
    &transN, &transT,
    &n, &m, &k,
    &alpha,
    A.memptr(), &n,
    B.memptr(), &m,
    &beta,
    C.memptr(), &n
  );
}

inline void dgemv_n(
    const arma::mat& A,
    const arma::vec& x,
    arma::vec& y,
    double alpha = 1.0,
    double beta = 0.0) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int m = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int inc = 1;
  const char transN = 'N';
  y.set_size(A.n_rows);
  arma::blas::gemv<double>(
    &transN,
    &n, &m,
    &alpha,
    A.memptr(), &n,
    x.memptr(), &inc,
    &beta,
    y.memptr(), &inc
  );
}

inline void dgemv_t(
    const arma::mat& A,
    const arma::vec& x,
    arma::vec& y,
    double alpha = 1.0,
    double beta = 0.0) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int m = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int inc = 1;
  const char transT = 'T';
  y.set_size(A.n_cols);
  arma::blas::gemv<double>(
    &transT,
    &n, &m,
    &alpha,
    A.memptr(), &n,
    x.memptr(), &inc,
    &beta,
    y.memptr(), &inc
  );
}

inline void sgemv_n(
    const arma::fmat& A,
    const arma::fvec& x,
    arma::fvec& y,
    float alpha = 1.0f,
    float beta = 0.0f) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int m = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int inc = 1;
  const char transN = 'N';
  y.set_size(A.n_rows);
  arma::blas::gemv<float>(
    &transN,
    &n, &m,
    &alpha,
    A.memptr(), &n,
    x.memptr(), &inc,
    &beta,
    y.memptr(), &inc
  );
}

inline void sgemv_t(
    const arma::fmat& A,
    const arma::fvec& x,
    arma::fvec& y,
    float alpha = 1.0f,
    float beta = 0.0f) {
  const arma::blas_int n = static_cast<arma::blas_int>(A.n_rows);
  const arma::blas_int m = static_cast<arma::blas_int>(A.n_cols);
  const arma::blas_int inc = 1;
  const char transT = 'T';
  y.set_size(A.n_cols);
  arma::blas::gemv<float>(
    &transT,
    &n, &m,
    &alpha,
    A.memptr(), &n,
    x.memptr(), &inc,
    &beta,
    y.memptr(), &inc
  );
}

template <typename T>
inline void matvec_colmajor(const arma::Mat<T>& A, const arma::Col<T>& x, arma::Col<T>& y) {
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  y.zeros(n);
  const T* ap = A.memptr();
  const T* xp = x.memptr();
  T* yp = y.memptr();
  for (arma::uword j = 0; j < m; ++j) {
    const T xj = xp[j];
    const T* col = ap + (j * n);
    for (arma::uword i = 0; i < n; ++i) {
      yp[i] += col[i] * xj;
    }
  }
}

template <typename T>
inline void tmatvec_colmajor(const arma::Mat<T>& A, const arma::Col<T>& x, arma::Col<T>& y) {
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  y.set_size(m);
  const T* ap = A.memptr();
  const T* xp = x.memptr();
  T* yp = y.memptr();
  for (arma::uword j = 0; j < m; ++j) {
    const T* col = ap + (j * n);
    T acc = 0;
    for (arma::uword i = 0; i < n; ++i) {
      acc += col[i] * xp[i];
    }
    yp[j] = acc;
  }
}

template <typename T>
inline void matvec_colmajor_blocked(
    const arma::Mat<T>& A,
    const arma::Col<T>& x,
    arma::Col<T>& y,
    arma::uword block_rows = kMatvecBlockRows,
    arma::uword block_cols = kMatvecBlockCols) {
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  y.zeros(n);
  const T* ap = A.memptr();
  const T* xp = x.memptr();
  T* yp = y.memptr();

#ifdef _OPENMP
  const bool use_omp = can_use_inner_omp(n * m, kMatvecOmpMinWork);
  if (use_omp) {
#pragma omp parallel for schedule(static)
    for (arma::sword ibs = 0; ibs < static_cast<arma::sword>(n); ibs += static_cast<arma::sword>(block_rows)) {
      const arma::uword ib = static_cast<arma::uword>(ibs);
      const arma::uword iend = std::min(n, ib + block_rows);
      const arma::uword len = iend - ib;
      T* yblk = yp + ib;
      for (arma::uword jb = 0; jb < m; jb += block_cols) {
        const arma::uword jend = std::min(m, jb + block_cols);
        for (arma::uword j = jb; j < jend; ++j) {
          const T xj = xp[j];
          const T* col = ap + (j * n) + ib;
#ifdef _OPENMP
#pragma omp simd
#endif
          for (arma::uword i = 0; i < len; ++i) {
            yblk[i] += col[i] * xj;
          }
        }
      }
    }
    return;
  }
#endif

  for (arma::uword jb = 0; jb < m; jb += block_cols) {
    const arma::uword jend = std::min(m, jb + block_cols);
    for (arma::uword ib = 0; ib < n; ib += block_rows) {
      const arma::uword iend = std::min(n, ib + block_rows);
      const arma::uword len = iend - ib;
      T* yblk = yp + ib;
      for (arma::uword j = jb; j < jend; ++j) {
        const T xj = xp[j];
        const T* col = ap + (j * n) + ib;
#ifdef _OPENMP
#pragma omp simd
#endif
        for (arma::uword i = 0; i < len; ++i) {
          yblk[i] += col[i] * xj;
        }
      }
    }
  }
}

template <typename T>
inline void tmatvec_colmajor_blocked(
    const arma::Mat<T>& A,
    const arma::Col<T>& x,
    arma::Col<T>& y,
    arma::uword block_rows = kMatvecBlockRows) {
  const arma::uword n = A.n_rows;
  const arma::uword m = A.n_cols;
  y.set_size(m);
  const T* ap = A.memptr();
  const T* xp = x.memptr();
  T* yp = y.memptr();

#ifdef _OPENMP
  const bool use_omp = can_use_inner_omp(n * m, kMatvecOmpMinWork);
  if (use_omp) {
#pragma omp parallel for schedule(static)
    for (arma::sword js = 0; js < static_cast<arma::sword>(m); ++js) {
      const arma::uword j = static_cast<arma::uword>(js);
      const T* col = ap + (j * n);
      T acc = 0;
      for (arma::uword ib = 0; ib < n; ib += block_rows) {
        const arma::uword iend = std::min(n, ib + block_rows);
#ifdef _OPENMP
#pragma omp simd reduction(+:acc)
#endif
        for (arma::uword i = ib; i < iend; ++i) {
          acc += col[i] * xp[i];
        }
      }
      yp[j] = acc;
    }
    return;
  }
#endif

  for (arma::uword j = 0; j < m; ++j) {
    const T* col = ap + (j * n);
    T acc = 0;
    for (arma::uword ib = 0; ib < n; ib += block_rows) {
      const arma::uword iend = std::min(n, ib + block_rows);
#ifdef _OPENMP
#pragma omp simd reduction(+:acc)
#endif
      for (arma::uword i = ib; i < iend; ++i) {
        acc += col[i] * xp[i];
      }
    }
    yp[j] = acc;
  }
}

struct MatvecDispatchThresholds {
  arma::uword blocked_min_d;
  arma::uword gemv_min_d;
  arma::uword blocked_min_f;
  arma::uword gemv_min_f;
};

template <typename Fn>
inline double bench_min_ms(Fn&& fn, int reps = 4) {
  fn();  // warm-up
  double best = std::numeric_limits<double>::infinity();
  for (int i = 0; i < reps; ++i) {
    const auto t0 = std::chrono::steady_clock::now();
    fn();
    const auto t1 = std::chrono::steady_clock::now();
    const std::chrono::duration<double, std::milli> dt = t1 - t0;
    best = std::min(best, dt.count());
  }
  return best;
}

inline void fill_tune_matrix(arma::mat& A, arma::vec& x) {
  const arma::uword n = A.n_rows;
  for (arma::uword j = 0; j < n; ++j) {
    for (arma::uword i = 0; i < n; ++i) {
      const double v = static_cast<double>(((i + 3u) * (j + 5u)) % 97u);
      A(i, j) = (v + 1.0) / 98.0;
    }
  }
  for (arma::uword i = 0; i < n; ++i) {
    const double v = static_cast<double>((i + 11u) % 31u);
    x(i) = (v + 1.0) / 32.0;
  }
}

inline void fill_tune_matrix_f(arma::fmat& A, arma::fvec& x) {
  const arma::uword n = A.n_rows;
  for (arma::uword j = 0; j < n; ++j) {
    for (arma::uword i = 0; i < n; ++i) {
      const float v = static_cast<float>(((i + 7u) * (j + 9u)) % 89u);
      A(i, j) = (v + 1.0f) / 90.0f;
    }
  }
  for (arma::uword i = 0; i < n; ++i) {
    const float v = static_cast<float>((i + 13u) % 29u);
    x(i) = (v + 1.0f) / 30.0f;
  }
}

inline void autotune_matvec_thresholds_double(arma::uword& blocked_min, arma::uword& gemv_min) {
  const std::vector<arma::uword> sizes = {48u, 64u, 80u, 96u, 128u, 160u, 192u, 256u};
  bool blocked_set = false;
  bool gemv_set = false;
  volatile double sink = 0.0;

  for (arma::uword n : sizes) {
    arma::mat A(n, n);
    arma::vec x(n);
    arma::vec y(n);
    arma::vec z(n);
    fill_tune_matrix(A, x);

    const double t_naive = bench_min_ms([&]() {
      matvec_colmajor(A, x, y);
      tmatvec_colmajor(A, y, z);
      sink += z(0);
    });
    const double t_blocked = bench_min_ms([&]() {
      matvec_colmajor_blocked(A, x, y);
      tmatvec_colmajor_blocked(A, y, z);
      sink += z(0);
    });
    const double t_blas = bench_min_ms([&]() {
      dgemv_n(A, x, y);
      dgemv_t(A, y, z);
      sink += z(0);
    });

    const arma::uword work = n * n;
    if (!blocked_set && t_blocked <= (0.98 * t_naive)) {
      blocked_min = work;
      blocked_set = true;
    }
    if (!gemv_set && t_blas <= (0.98 * std::min(t_naive, t_blocked))) {
      gemv_min = work;
      gemv_set = true;
    }
  }
  (void)sink;
}

inline void autotune_matvec_thresholds_float(arma::uword& blocked_min, arma::uword& gemv_min) {
  const std::vector<arma::uword> sizes = {48u, 64u, 80u, 96u, 128u, 160u, 192u, 256u};
  bool blocked_set = false;
  bool gemv_set = false;
  volatile float sink = 0.0f;

  for (arma::uword n : sizes) {
    arma::fmat A(n, n);
    arma::fvec x(n);
    arma::fvec y(n);
    arma::fvec z(n);
    fill_tune_matrix_f(A, x);

    const double t_naive = bench_min_ms([&]() {
      matvec_colmajor(A, x, y);
      tmatvec_colmajor(A, y, z);
      sink += z(0);
    });
    const double t_blocked = bench_min_ms([&]() {
      matvec_colmajor_blocked(A, x, y);
      tmatvec_colmajor_blocked(A, y, z);
      sink += z(0);
    });
    const double t_blas = bench_min_ms([&]() {
      sgemv_n(A, x, y);
      sgemv_t(A, y, z);
      sink += z(0);
    });

    const arma::uword work = n * n;
    if (!blocked_set && t_blocked <= (0.98 * t_naive)) {
      blocked_min = work;
      blocked_set = true;
    }
    if (!gemv_set && t_blas <= (0.98 * std::min(t_naive, t_blocked))) {
      gemv_min = work;
      gemv_set = true;
    }
  }
  (void)sink;
}

inline MatvecDispatchThresholds init_matvec_dispatch_thresholds() {
  MatvecDispatchThresholds t{
    kMatvecBlockedMinWorkD,
    kGemvMinWorkD,
    kMatvecBlockedMinWorkF,
    kGemvMinWorkF
  };

  t.blocked_min_d = parse_uword_env("RFUGW_MATVEC_BLOCKED_MIN_WORK_D", t.blocked_min_d);
  t.gemv_min_d = parse_uword_env("RFUGW_MATVEC_GEMV_MIN_WORK_D", t.gemv_min_d);
  t.blocked_min_f = parse_uword_env("RFUGW_MATVEC_BLOCKED_MIN_WORK_F", t.blocked_min_f);
  t.gemv_min_f = parse_uword_env("RFUGW_MATVEC_GEMV_MIN_WORK_F", t.gemv_min_f);

  if (parse_uword_env("RFUGW_AUTOTUNE_MATVEC", 1u) != 0u) {
    autotune_matvec_thresholds_double(t.blocked_min_d, t.gemv_min_d);
    autotune_matvec_thresholds_float(t.blocked_min_f, t.gemv_min_f);
  }

  if (t.gemv_min_d < t.blocked_min_d) {
    t.gemv_min_d = t.blocked_min_d;
  }
  if (t.gemv_min_f < t.blocked_min_f) {
    t.gemv_min_f = t.blocked_min_f;
  }

  return t;
}

inline const MatvecDispatchThresholds& matvec_dispatch_thresholds() {
  static const MatvecDispatchThresholds t = init_matvec_dispatch_thresholds();
  return t;
}

inline arma::uword gemv_min_work_d() {
  return matvec_dispatch_thresholds().gemv_min_d;
}

inline arma::uword matvec_blocked_min_work_d() {
  return matvec_dispatch_thresholds().blocked_min_d;
}

inline arma::uword gemv_min_work_f() {
  return matvec_dispatch_thresholds().gemv_min_f;
}

inline arma::uword matvec_blocked_min_work_f() {
  return matvec_dispatch_thresholds().blocked_min_f;
}

inline void build_kernel_from_cost(
    const arma::mat& cost,
    double inv_epsilon,
    arma::mat& K) {
  K.set_size(cost.n_rows, cost.n_cols);
  const arma::uword n = cost.n_elem;
  const double* cp = cost.memptr();
  double* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      double z = -inv_epsilon * cp[static_cast<arma::uword>(i)];
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[static_cast<arma::uword>(i)] = std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      double z = -inv_epsilon * cp[i];
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[i] = std::exp(z);
    }
  }
}

inline void build_kernel_from_cost_f(
    const arma::fmat& cost,
    float inv_epsilon,
    arma::fmat& K) {
  K.set_size(cost.n_rows, cost.n_cols);
  const arma::uword n = cost.n_elem;
  const float* cp = cost.memptr();
  float* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      float z = -inv_epsilon * cp[static_cast<arma::uword>(i)];
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[static_cast<arma::uword>(i)] = std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      float z = -inv_epsilon * cp[i];
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[i] = std::exp(z);
    }
  }
}

inline void build_kernel_from_cost_mul(
    const arma::mat& cost,
    const arma::mat& c,
    double inv_epsilon,
    arma::mat& K) {
  K.set_size(cost.n_rows, cost.n_cols);
  const arma::uword n = cost.n_elem;
  const double* cp = cost.memptr();
  const double* wp = c.memptr();
  double* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      double z = -inv_epsilon * cp[iu];
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[iu] = wp[iu] * std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      double z = -inv_epsilon * cp[i];
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[i] = wp[i] * std::exp(z);
    }
  }
}

inline void build_kernel_from_cost_mul_f(
    const arma::fmat& cost,
    const arma::fmat& c,
    float inv_epsilon,
    arma::fmat& K) {
  K.set_size(cost.n_rows, cost.n_cols);
  const arma::uword n = cost.n_elem;
  const float* cp = cost.memptr();
  const float* wp = c.memptr();
  float* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      float z = -inv_epsilon * cp[iu];
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[iu] = wp[iu] * std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      float z = -inv_epsilon * cp[i];
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[i] = wp[i] * std::exp(z);
    }
  }
}

inline void build_kernel_from_fused_cost_mul(
    const arma::mat& quad_cost,
    const arma::mat& lin_cost,
    const arma::mat& c,
    double inv_epsilon,
    double quad_scale,
    double lin_scale,
    arma::mat& K) {
  K.set_size(quad_cost.n_rows, quad_cost.n_cols);
  const arma::uword n = quad_cost.n_elem;
  const double* qp = quad_cost.memptr();
  const double* lp = lin_cost.memptr();
  const double* wp = c.memptr();
  double* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      double z = -inv_epsilon * (quad_scale * qp[iu] + lin_scale * lp[iu]);
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[iu] = wp[iu] * std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      double z = -inv_epsilon * (quad_scale * qp[i] + lin_scale * lp[i]);
      z = std::max(-kExpClip, std::min(kExpClip, z));
      kp[i] = wp[i] * std::exp(z);
    }
  }
}

inline void build_kernel_from_fused_cost_mul_f(
    const arma::fmat& quad_cost,
    const arma::fmat& lin_cost,
    const arma::fmat& c,
    float inv_epsilon,
    float quad_scale,
    float lin_scale,
    arma::fmat& K) {
  K.set_size(quad_cost.n_rows, quad_cost.n_cols);
  const arma::uword n = quad_cost.n_elem;
  const float* qp = quad_cost.memptr();
  const float* lp = lin_cost.memptr();
  const float* wp = c.memptr();
  float* kp = K.memptr();
  if (can_use_inner_omp(n, kKernelOmpMinWork)) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
    for (arma::sword i = 0; i < static_cast<arma::sword>(n); ++i) {
      const arma::uword iu = static_cast<arma::uword>(i);
      float z = -inv_epsilon * (quad_scale * qp[iu] + lin_scale * lp[iu]);
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[iu] = wp[iu] * std::exp(z);
    }
#endif
  } else {
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      float z = -inv_epsilon * (quad_scale * qp[i] + lin_scale * lp[i]);
      z = std::max(-kExpClipF, std::min(kExpClipF, z));
      kp[i] = wp[i] * std::exp(z);
    }
  }
}

inline void init_matrices_square(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    arma::mat& constC,
    arma::mat& hC1,
    arma::mat& hC2) {
  const arma::mat fC1 = C1 % C1;
  const arma::mat fC2 = C2 % C2;
  hC1 = C1;
  hC2 = 2.0 * C2;

  const arma::vec left = fC1 * p;
  const arma::rowvec right = q.t() * fC2.t();
  constC.set_size(C1.n_rows, C2.n_rows);
  for (arma::uword j = 0; j < C2.n_rows; ++j) {
    constC.col(j) = left;
    constC.col(j) += right[j];
  }
}

inline void init_matrices_square_f(
    const arma::fmat& C1,
    const arma::fmat& C2,
    const arma::fvec& p,
    const arma::fvec& q,
    arma::fmat& constC,
    arma::fmat& hC1,
    arma::fmat& hC2) {
  const arma::fmat fC1 = C1 % C1;
  const arma::fmat fC2 = C2 % C2;
  hC1 = C1;
  hC2 = 2.0f * C2;

  const arma::fvec left = fC1 * p;
  const arma::frowvec right = q.t() * fC2.t();
  constC.set_size(C1.n_rows, C2.n_rows);
  for (arma::uword j = 0; j < C2.n_rows; ++j) {
    constC.col(j) = left;
    constC.col(j) += right[j];
  }
}

inline arma::mat tensor_product(
    const arma::mat& constC,
    const arma::mat& hC1,
    const arma::mat& hC2,
    const arma::mat& T) {
  arma::mat out, scratch;
  dgemm_nn(hC1, T, scratch);
  out = constC;
  dgemm_nt_accum(scratch, hC2, out, -1.0, 1.0);
  return out;
}

inline void tensor_product_blas_scaled(
    const arma::mat& base,
    const arma::mat& hC1,
    const arma::mat& hC2,
    const arma::mat& T,
    double cross_scale,
    arma::mat& out,
    arma::mat& scratch) {
  dgemm_nn(hC1, T, scratch);
  out = base;
  dgemm_nt_accum(scratch, hC2, out, cross_scale, 1.0);
}

inline void tensor_product_asym_blas_scaled(
    const arma::mat& base,
    const arma::mat& hC1,
    const arma::mat& hC2,
    const arma::mat& hC1t,
    const arma::mat& hC2_asym,
    const arma::mat& T,
    double cross_scale,
    arma::mat& out,
    arma::mat& scratch,
    arma::mat& scratch_asym) {
  dgemm_nn(hC1, T, scratch);
  dgemm_nn(hC1t, T, scratch_asym);
  out = base;
  dgemm_nt_accum(scratch, hC2, out, cross_scale, 1.0);
  dgemm_nt_accum(scratch_asym, hC2_asym, out, cross_scale, 1.0);
}

inline void tensor_product_blas_scaled_f(
    const arma::fmat& base,
    const arma::fmat& hC1,
    const arma::fmat& hC2,
    const arma::fmat& T,
    float cross_scale,
    arma::fmat& out,
    arma::fmat& scratch) {
  sgemm_nn(hC1, T, scratch);
  out = base;
  sgemm_nt_accum(scratch, hC2, out, cross_scale, 1.0f);
}

inline void tensor_product_asym_blas_scaled_f(
    const arma::fmat& base,
    const arma::fmat& hC1,
    const arma::fmat& hC2,
    const arma::fmat& hC1t,
    const arma::fmat& hC2_asym,
    const arma::fmat& T,
    float cross_scale,
    arma::fmat& out,
    arma::fmat& scratch,
    arma::fmat& scratch_asym) {
  sgemm_nn(hC1, T, scratch);
  sgemm_nn(hC1t, T, scratch_asym);
  out = base;
  sgemm_nt_accum(scratch, hC2, out, cross_scale, 1.0f);
  sgemm_nt_accum(scratch_asym, hC2_asym, out, cross_scale, 1.0f);
}

inline arma::mat cross_feature_cost_matrix(
    const arma::mat& X,
    const arma::mat& Y,
    const arma::vec& y2,
    bool use_euclidean,
    bool rescale01) {
  arma::mat gram;
  dgemm_nt(X, Y, gram);
  const arma::vec x2 = arma::sum(arma::square(X), 1);
  arma::mat out(X.n_rows, Y.n_rows);
  for (arma::uword j = 0; j < Y.n_rows; ++j) {
    const double y2j = y2[j];
    const double* gcol = gram.colptr(j);
    double* ocol = out.colptr(j);
    for (arma::uword i = 0; i < X.n_rows; ++i) {
      double v = x2[i] + y2j - 2.0 * gcol[i];
      if (v < 0.0) {
        v = 0.0;
      }
      ocol[i] = use_euclidean ? std::sqrt(v) : v;
    }
  }
  if (rescale01) {
    const double mx = out.max();
    if (std::isfinite(mx) && mx > 0.0) {
      out /= mx;
    }
  }
  return out;
}

inline bool svd_lowrank_factors(
    const arma::mat& C,
    int rank,
    arma::mat& A_scaled,
    arma::mat& Bt) {
  if (rank <= 0) {
    return false;
  }
  arma::mat U, V;
  arma::vec s;
  if (!arma::svd_econ(U, s, V, C)) {
    return false;
  }
  if (s.n_elem == 0) {
    return false;
  }
  const arma::uword r = std::min<arma::uword>(static_cast<arma::uword>(rank), s.n_elem);
  if (r == 0) {
    return false;
  }
  const arma::mat Ur = U.cols(0, r - 1);
  const arma::mat Vr = V.cols(0, r - 1);
  const arma::rowvec sr = s.subvec(0, r - 1).t();
  A_scaled = Ur.each_row() % sr;
  Bt = Vr.t();
  return true;
}

inline bool svd_lowrank_factors_f(
    const arma::fmat& C,
    int rank,
    arma::fmat& A_scaled,
    arma::fmat& Bt) {
  if (rank <= 0) {
    return false;
  }
  arma::fmat U, V;
  arma::fvec s;
  if (!arma::svd_econ(U, s, V, C)) {
    return false;
  }
  if (s.n_elem == 0) {
    return false;
  }
  const arma::uword r = std::min<arma::uword>(static_cast<arma::uword>(rank), s.n_elem);
  if (r == 0) {
    return false;
  }
  const arma::fmat Ur = U.cols(0, r - 1);
  const arma::fmat Vr = V.cols(0, r - 1);
  const arma::frowvec sr = s.subvec(0, r - 1).t();
  A_scaled = Ur.each_row() % sr;
  Bt = Vr.t();
  return true;
}

inline void tensor_product_lowrank_scaled(
    const arma::mat& base,
    const arma::mat& A1_scaled,
    const arma::mat& B1t,
    const arma::mat& A2t_scaled,
    const arma::mat& B2t,
    const arma::mat& T,
    double cross_scale,
    arma::mat& out,
  arma::mat& tmp_r1_nt,
  arma::mat& tmp_r1_r2,
  arma::mat& tmp_ns_r2) {
  dgemm_nn(B1t, T, tmp_r1_nt);
  dgemm_nt(tmp_r1_nt, B2t, tmp_r1_r2);
  dgemm_nn(A1_scaled, tmp_r1_r2, tmp_ns_r2);
  out = base;
  dgemm_nn_accum(tmp_ns_r2, A2t_scaled, out, cross_scale, 1.0);
}

inline void tensor_product_lowrank_scaled_f(
    const arma::fmat& base,
    const arma::fmat& A1_scaled,
    const arma::fmat& B1t,
    const arma::fmat& A2t_scaled,
    const arma::fmat& B2t,
    const arma::fmat& T,
    float cross_scale,
    arma::fmat& out,
  arma::fmat& tmp_r1_nt,
  arma::fmat& tmp_r1_r2,
  arma::fmat& tmp_ns_r2) {
  sgemm_nn(B1t, T, tmp_r1_nt);
  sgemm_nt(tmp_r1_nt, B2t, tmp_r1_r2);
  sgemm_nn(A1_scaled, tmp_r1_r2, tmp_ns_r2);
  out = base;
  sgemm_nn_accum(tmp_ns_r2, A2t_scaled, out, cross_scale, 1.0f);
}

inline double gw_loss(
    const arma::mat& constC,
    const arma::mat& hC1,
    const arma::mat& hC2,
    const arma::mat& T) {
  const arma::mat tens = tensor_product(constC, hC1, hC2, T);
  return arma::accu(tens % T);
}

struct SinkhornBalancedResult {
  arma::mat plan;
  arma::vec u;
  arma::vec v;
  arma::vec f;
  arma::vec g;
  int iters;
  double err;
};

struct SinkhornBalancedResultF {
  arma::fmat plan;
  arma::fvec u;
  arma::fvec v;
  int iters;
  float err;
};

struct LowRankC2CacheD {
  bool valid = false;
  arma::mat A2t_scaled;
  arma::mat B2t;
};

struct LowRankC2CacheF {
  bool valid = false;
  arma::fmat A2t_scaled;
  arma::fmat B2t;
};

struct SquareC2CacheD {
  bool valid = false;
  arma::mat hC2;
  arma::rowvec right_term;
};

struct SquareC2CacheF {
  bool valid = false;
  arma::fmat C2;
  arma::fvec q;
  arma::fmat hC2;
  arma::frowvec right_term;
};

inline LowRankC2CacheD build_lowrank_c2_cache_d(
    const arma::mat& C2,
    int approx_rank) {
  LowRankC2CacheD cache;
  if (approx_rank <= 0) {
    return cache;
  }
  arma::mat A2_scaled;
  if (!svd_lowrank_factors(C2, approx_rank, A2_scaled, cache.B2t)) {
    return cache;
  }
  cache.A2t_scaled = A2_scaled.t();
  cache.valid = true;
  return cache;
}

inline LowRankC2CacheF build_lowrank_c2_cache_f(
    const arma::fmat& C2,
    int approx_rank) {
  LowRankC2CacheF cache;
  if (approx_rank <= 0) {
    return cache;
  }
  arma::fmat A2_scaled;
  if (!svd_lowrank_factors_f(C2, approx_rank, A2_scaled, cache.B2t)) {
    return cache;
  }
  cache.A2t_scaled = A2_scaled.t();
  cache.valid = true;
  return cache;
}

inline bool lowrank_cache_compatible(
    const LowRankC2CacheD& cache,
    const arma::mat& C2) {
  return cache.valid &&
    cache.A2t_scaled.n_cols == C2.n_rows &&
    cache.B2t.n_cols == C2.n_cols &&
    cache.A2t_scaled.n_rows == cache.B2t.n_rows;
}

inline bool lowrank_cache_compatible_f(
    const LowRankC2CacheF& cache,
    const arma::fmat& C2) {
  return cache.valid &&
    cache.A2t_scaled.n_cols == C2.n_rows &&
    cache.B2t.n_cols == C2.n_cols &&
    cache.A2t_scaled.n_rows == cache.B2t.n_rows;
}

inline bool lowrank_c1_cache_compatible(
    const arma::mat& A1_scaled,
    const arma::mat& B1t,
    const arma::mat& C1) {
  return A1_scaled.n_rows == C1.n_rows &&
    B1t.n_cols == C1.n_cols &&
    A1_scaled.n_cols == B1t.n_rows &&
    A1_scaled.n_cols > 0;
}

inline bool lowrank_c1_cache_compatible_f(
    const arma::fmat& A1_scaled,
    const arma::fmat& B1t,
    const arma::fmat& C1) {
  return A1_scaled.n_rows == C1.n_rows &&
    B1t.n_cols == C1.n_cols &&
    A1_scaled.n_cols == B1t.n_rows &&
    A1_scaled.n_cols > 0;
}

inline SquareC2CacheD build_square_c2_cache_d(
    const arma::mat& C2,
    const arma::vec& q) {
  SquareC2CacheD cache;
  if (C2.n_rows == 0 || C2.n_rows != C2.n_cols || q.n_elem != C2.n_rows) {
    return cache;
  }
  const arma::mat fC2 = C2 % C2;
  cache.hC2 = 2.0 * C2;
  cache.right_term = q.t() * fC2.t();
  cache.valid = true;
  return cache;
}

inline SquareC2CacheF build_square_c2_cache_f(
    const arma::fmat& C2,
    const arma::fvec& q) {
  SquareC2CacheF cache;
  if (C2.n_rows == 0 || C2.n_rows != C2.n_cols || q.n_elem != C2.n_rows) {
    return cache;
  }
  const arma::fmat fC2 = C2 % C2;
  cache.C2 = C2;
  cache.q = q;
  cache.hC2 = 2.0f * C2;
  cache.right_term = q.t() * fC2.t();
  cache.valid = true;
  return cache;
}

inline bool square_cache_compatible(
    const SquareC2CacheD& cache,
    const arma::mat& C2,
    const arma::vec& q) {
  return cache.valid &&
    C2.n_rows == C2.n_cols &&
    q.n_elem == C2.n_rows &&
    cache.hC2.n_rows == C2.n_rows &&
    cache.hC2.n_cols == C2.n_cols &&
    cache.right_term.n_elem == C2.n_rows;
}

inline bool square_cache_compatible_f(
    const SquareC2CacheF& cache,
    const arma::mat& C2,
    const arma::vec& q) {
  return cache.valid &&
    C2.n_rows == C2.n_cols &&
    q.n_elem == C2.n_rows &&
    cache.C2.n_rows == C2.n_rows &&
    cache.C2.n_cols == C2.n_cols &&
    cache.q.n_elem == C2.n_rows &&
    cache.hC2.n_rows == C2.n_rows &&
    cache.hC2.n_cols == C2.n_cols &&
    cache.right_term.n_elem == C2.n_rows;
}

inline void init_matrices_square_from_cache_d(
    const arma::mat& C1,
    const arma::vec& p,
    const SquareC2CacheD& cache,
    arma::mat& constC,
    arma::mat& hC1) {
  const arma::mat fC1 = C1 % C1;
  hC1 = C1;
  const arma::vec left = fC1 * p;
  constC.set_size(C1.n_rows, cache.hC2.n_rows);
  for (arma::uword j = 0; j < cache.hC2.n_rows; ++j) {
    constC.col(j) = left;
    constC.col(j) += cache.right_term[j];
  }
}

inline void init_matrices_square_from_cache_f(
    const arma::fmat& C1,
    const arma::fvec& p,
    const SquareC2CacheF& cache,
    arma::fmat& constC,
    arma::fmat& hC1) {
  const arma::fmat fC1 = C1 % C1;
  hC1 = C1;
  const arma::fvec left = fC1 * p;
  constC.set_size(C1.n_rows, cache.hC2.n_rows);
  for (arma::uword j = 0; j < cache.hC2.n_rows; ++j) {
    constC.col(j) = left;
    constC.col(j) += cache.right_term[j];
  }
}

inline double logsumexp_vec(const arma::vec& x) {
  const double m = x.max();
  if (!std::isfinite(m)) {
    return m;
  }
  return m + std::log(arma::accu(arma::exp(x - m)));
}

inline SinkhornBalancedResult sinkhorn_balanced(
    const arma::vec& p,
    const arma::vec& q,
    const arma::mat& cost,
    double epsilon,
    int max_iter,
    double tol,
    arma::vec u,
    arma::vec v) {
  arma::mat K;
  build_kernel_from_cost(cost, 1.0 / epsilon, K);
  const bool use_blas = (K.n_elem >= gemv_min_work_d());
  const bool use_blocked = (!use_blas && K.n_elem >= matvec_blocked_min_work_d());
  const int check_interval = (K.n_elem <= 160000) ? 5 : 10;
  if (u.n_elem != p.n_elem) {
    u = arma::ones<arma::vec>(p.n_elem);
  }
  if (v.n_elem != q.n_elem) {
    v = arma::ones<arma::vec>(q.n_elem);
  }

  arma::vec Kv(p.n_elem);
  arma::vec Ktu(q.n_elem);
  if (use_blas) {
    dgemv_t(K, u, Ktu);
  } else if (use_blocked) {
    tmatvec_colmajor_blocked(K, u, Ktu);
  } else {
    tmatvec_colmajor(K, u, Ktu);
  }
  Ktu += kTiny;

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if (use_blas) {
      dgemv_n(K, v, Kv);
    } else if (use_blocked) {
      matvec_colmajor_blocked(K, v, Kv);
    } else {
      matvec_colmajor(K, v, Kv);
    }
    Kv += kTiny;

    if ((it + 1) % check_interval == 0) {
      const arma::vec row_marg = u % Kv;
      const arma::vec col_marg = v % Ktu;
      const double er = arma::max(arma::abs(row_marg - p));
      const double ec = arma::max(arma::abs(col_marg - q));
      err = std::max(er, ec);
      if (err < tol) {
        ++it;
        break;
      }
    }

    u = p / Kv;
    if (use_blas) {
      dgemv_t(K, u, Ktu);
    } else if (use_blocked) {
      tmatvec_colmajor_blocked(K, u, Ktu);
    } else {
      tmatvec_colmajor(K, u, Ktu);
    }
    Ktu += kTiny;
    v = q / Ktu;
  }

  SinkhornBalancedResult out;
  out.plan = (u * v.t()) % K;
  out.u = std::move(u);
  out.v = std::move(v);
  out.f = arma::vec();
  out.g = arma::vec();
  out.iters = it;
  out.err = err;
  return out;
}

inline SinkhornBalancedResultF sinkhorn_balanced_f(
    const arma::fvec& p,
    const arma::fvec& q,
    const arma::fmat& cost,
    float epsilon,
    int max_iter,
    float tol,
    arma::fvec u,
    arma::fvec v) {
  arma::fmat K;
  build_kernel_from_cost_f(cost, 1.0f / epsilon, K);
  const bool use_blas = (K.n_elem >= gemv_min_work_f());
  const bool use_blocked = (!use_blas && K.n_elem >= matvec_blocked_min_work_f());
  const int check_interval = (K.n_elem <= 160000) ? 5 : 10;
  if (u.n_elem != p.n_elem) {
    u = arma::ones<arma::fvec>(p.n_elem);
  }
  if (v.n_elem != q.n_elem) {
    v = arma::ones<arma::fvec>(q.n_elem);
  }

  arma::fvec Kv(p.n_elem);
  arma::fvec Ktu(q.n_elem);
  if (use_blas) {
    sgemv_t(K, u, Ktu);
  } else if (use_blocked) {
    tmatvec_colmajor_blocked(K, u, Ktu);
  } else {
    tmatvec_colmajor(K, u, Ktu);
  }
  Ktu += kTinyF;

  float err = std::numeric_limits<float>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if (use_blas) {
      sgemv_n(K, v, Kv);
    } else if (use_blocked) {
      matvec_colmajor_blocked(K, v, Kv);
    } else {
      matvec_colmajor(K, v, Kv);
    }
    Kv += kTinyF;

    if ((it + 1) % check_interval == 0) {
      const arma::fvec row_marg = u % Kv;
      const arma::fvec col_marg = v % Ktu;
      const float er = arma::max(arma::abs(row_marg - p));
      const float ec = arma::max(arma::abs(col_marg - q));
      err = std::max(er, ec);
      if (err < tol) {
        ++it;
        break;
      }
    }

    u = p / Kv;
    if (use_blas) {
      sgemv_t(K, u, Ktu);
    } else if (use_blocked) {
      tmatvec_colmajor_blocked(K, u, Ktu);
    } else {
      tmatvec_colmajor(K, u, Ktu);
    }
    Ktu += kTinyF;
    v = q / Ktu;
  }

  SinkhornBalancedResultF out;
  out.plan = (u * v.t()) % K;
  out.u = std::move(u);
  out.v = std::move(v);
  out.iters = it;
  out.err = err;
  return out;
}

inline SinkhornBalancedResult sinkhorn_balanced_log(
    const arma::vec& p,
    const arma::vec& q,
    const arma::mat& cost,
    double epsilon,
    int max_iter,
    double tol,
    arma::vec f,
    arma::vec g) {
  const arma::uword ns = p.n_elem;
  const arma::uword nt = q.n_elem;
  if (f.n_elem != ns) {
    f = arma::zeros<arma::vec>(ns);
  }
  if (g.n_elem != nt) {
    g = arma::zeros<arma::vec>(nt);
  }

  const arma::vec logp = arma::log(p + kTiny);
  const arma::vec logq = arma::log(q + kTiny);

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    for (arma::uword i = 0; i < ns; ++i) {
      const arma::vec z = (g - cost.row(i).t()) / epsilon;
      f(i) = epsilon * (logp(i) - logsumexp_vec(z));
    }
    for (arma::uword j = 0; j < nt; ++j) {
      const arma::vec z = (f - cost.col(j)) / epsilon;
      g(j) = epsilon * (logq(j) - logsumexp_vec(z));
    }

    if ((it + 1) % 10 == 0) {
      arma::vec row_marg(ns);
      arma::vec col_marg(nt);
      for (arma::uword i = 0; i < ns; ++i) {
        const arma::vec z = (g - cost.row(i).t()) / epsilon;
        row_marg(i) = std::exp((f(i) / epsilon) + logsumexp_vec(z));
      }
      for (arma::uword j = 0; j < nt; ++j) {
        const arma::vec z = (f - cost.col(j)) / epsilon;
        col_marg(j) = std::exp((g(j) / epsilon) + logsumexp_vec(z));
      }
      const double er = arma::max(arma::abs(row_marg - p));
      const double ec = arma::max(arma::abs(col_marg - q));
      err = std::max(er, ec);
      if (err < tol) {
        ++it;
        break;
      }
    }
  }

  arma::mat plan(ns, nt);
  for (arma::uword i = 0; i < ns; ++i) {
    for (arma::uword j = 0; j < nt; ++j) {
      const double z = (f(i) + g(j) - cost(i, j)) / epsilon;
      const double zc = std::max(-kExpClip, std::min(kExpClip, z));
      plan(i, j) = std::exp(zc);
    }
  }

  SinkhornBalancedResult out;
  out.plan = std::move(plan);
  out.u = arma::vec();
  out.v = arma::vec();
  out.f = std::move(f);
  out.g = std::move(g);
  out.iters = it;
  out.err = err;
  return out;
}

inline double kl_div(const arma::mat& p, const arma::mat& q, bool mass) {
  const arma::mat ratio = (p + kTiny) / (q + kTiny);
  double value = arma::accu(p % arma::log(ratio));
  if (mass) {
    value += arma::accu(q - p);
  }
  return value;
}

inline double kl_div(const arma::vec& p, const arma::vec& q, bool mass) {
  const arma::vec ratio = (p + kTiny) / (q + kTiny);
  double value = arma::accu(p % arma::log(ratio));
  if (mass) {
    value += arma::accu(q - p);
  }
  return value;
}

inline double div_between_product_kl(
    const arma::vec& mu,
    const arma::vec& nu,
    const arma::vec& alpha,
    const arma::vec& beta) {
  const double m_mu = arma::accu(mu);
  const double m_nu = arma::accu(nu);
  const double m_alpha = arma::accu(alpha);
  const double m_beta = arma::accu(beta);
  const double cst = (m_mu - m_alpha) * (m_nu - m_beta);
  return m_nu * kl_div(mu, alpha, true) + m_mu * kl_div(nu, beta, true) + cst;
}

inline double div_between_product_kl(
    const arma::mat& mu,
    const arma::mat& nu,
    const arma::mat& alpha,
    const arma::mat& beta) {
  const double m_mu = arma::accu(mu);
  const double m_nu = arma::accu(nu);
  const double m_alpha = arma::accu(alpha);
  const double m_beta = arma::accu(beta);
  const double cst = (m_mu - m_alpha) * (m_nu - m_beta);
  return m_nu * kl_div(mu, alpha, true) + m_mu * kl_div(nu, beta, true) + cst;
}

inline double div_to_product_kl(
    const arma::mat& pi,
    const arma::vec& a,
    const arma::vec& b,
    const arma::vec& pi1,
    const arma::vec& pi2,
    bool mass) {
  double res = arma::accu(pi % arma::log(pi + kTiny)) -
    arma::accu(pi1 % arma::log(a + kTiny)) -
    arma::accu(pi2 % arma::log(b + kTiny));
  if (mass) {
    res += -arma::accu(pi1) + arma::accu(a) * arma::accu(b);
  }
  return res;
}

inline double kl_scalar_shift_from_marginals_and_plan(
    const arma::mat& pi,
    const arma::vec& pi1,
    const arma::vec& pi2,
    const arma::vec& log_a,
    const arma::vec& log_b,
    double rho_x,
    double rho_y,
    double eps) {
  double pi_log_pi = 0.0;
  const double* pip = pi.memptr();
  const arma::uword n_elem = pi.n_elem;
  for (arma::uword idx = 0; idx < n_elem; ++idx) {
    const double v = pip[idx];
    pi_log_pi += v * std::log(v + kTiny);
  }

  double kl_x = 0.0;
  double kl_y = 0.0;
  double pi1_log_a = 0.0;
  double pi2_log_b = 0.0;
  const double* p1p = pi1.memptr();
  const double* p2p = pi2.memptr();
  const double* lap = log_a.memptr();
  const double* lbp = log_b.memptr();

  for (arma::uword i = 0; i < pi1.n_elem; ++i) {
    const double v = p1p[i];
    const double lv = std::log(v + kTiny);
    const double la = lap[i];
    kl_x += v * (lv - la);
    pi1_log_a += v * la;
  }

  for (arma::uword j = 0; j < pi2.n_elem; ++j) {
    const double v = p2p[j];
    const double lv = std::log(v + kTiny);
    const double lb = lbp[j];
    kl_y += v * (lv - lb);
    pi2_log_b += v * lb;
  }

  double scalar_shift = 0.0;
  if (std::isfinite(rho_x) && rho_x != 0.0) {
    scalar_shift += rho_x * kl_x;
  }
  if (std::isfinite(rho_y) && rho_y != 0.0) {
    scalar_shift += rho_y * kl_y;
  }
  if (eps > 0.0) {
    scalar_shift += eps * (pi_log_pi - pi1_log_a - pi2_log_b);
  }
  return scalar_shift;
}

inline void finalize_uot_cost_affine_inplace(
    arma::mat& uot_cost,
    const arma::mat& M,
    const arma::vec& A,
    const arma::vec& B,
    double scalar_shift) {
  const arma::uword n = uot_cost.n_rows;
  const arma::uword m = uot_cost.n_cols;
  const double* mp = M.memptr();
  const double* ap = A.memptr();
  const double* bp = B.memptr();
  double* up = uot_cost.memptr();

  for (arma::uword j = 0; j < m; ++j) {
    const double bj = bp[j] + scalar_shift;
    const arma::uword col_offset = j * n;
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      const arma::uword idx = col_offset + i;
      up[idx] = -2.0 * up[idx] + mp[idx] + ap[i] + bj;
    }
  }
}

inline void uot_cost_matrix_kl_joint_inplace(
    const arma::mat& X_sqr,
    const arma::mat& Y_sqr,
    const arma::mat& X,
    const arma::mat& Yt,
    const arma::mat& M,
    const arma::mat& pi,
    const arma::vec& a,
    const arma::vec& b,
    const arma::vec& log_a,
    const arma::vec& log_b,
    double rho_x,
    double rho_y,
    double eps,
    arma::mat& uot_cost,
    arma::mat& scratch,
    arma::vec& pi1,
    arma::vec& pi2,
    arma::vec& A,
    arma::vec& B) {
  pi1 = arma::sum(pi, 1);
  pi2 = arma::sum(pi, 0).t();
  A = X_sqr * pi1;
  B = Y_sqr * pi2;
  dgemm_nn(X, pi, scratch);
  dgemm_nn(scratch, Yt, uot_cost);

  (void)a;
  (void)b;
  const double scalar_shift = kl_scalar_shift_from_marginals_and_plan(
    pi, pi1, pi2, log_a, log_b, rho_x, rho_y, eps
  );

  finalize_uot_cost_affine_inplace(uot_cost, M, A, B, scalar_shift);
}

inline float kl_scalar_shift_from_marginals_and_plan_f(
    const arma::fmat& pi,
    const arma::fvec& pi1,
    const arma::fvec& pi2,
    const arma::fvec& log_a,
    const arma::fvec& log_b,
    float rho_x,
    float rho_y,
    float eps) {
  float pi_log_pi = 0.0f;
  const float* pip = pi.memptr();
  const arma::uword n_elem = pi.n_elem;
  for (arma::uword idx = 0; idx < n_elem; ++idx) {
    const float v = pip[idx];
    pi_log_pi += v * std::log(v + kTinyF);
  }

  float kl_x = 0.0f;
  float kl_y = 0.0f;
  float pi1_log_a = 0.0f;
  float pi2_log_b = 0.0f;
  const float* p1p = pi1.memptr();
  const float* p2p = pi2.memptr();
  const float* lap = log_a.memptr();
  const float* lbp = log_b.memptr();

  for (arma::uword i = 0; i < pi1.n_elem; ++i) {
    const float v = p1p[i];
    const float lv = std::log(v + kTinyF);
    const float la = lap[i];
    kl_x += v * (lv - la);
    pi1_log_a += v * la;
  }

  for (arma::uword j = 0; j < pi2.n_elem; ++j) {
    const float v = p2p[j];
    const float lv = std::log(v + kTinyF);
    const float lb = lbp[j];
    kl_y += v * (lv - lb);
    pi2_log_b += v * lb;
  }

  float scalar_shift = 0.0f;
  if (std::isfinite(rho_x) && rho_x != 0.0f) {
    scalar_shift += rho_x * kl_x;
  }
  if (std::isfinite(rho_y) && rho_y != 0.0f) {
    scalar_shift += rho_y * kl_y;
  }
  if (eps > 0.0f) {
    scalar_shift += eps * (pi_log_pi - pi1_log_a - pi2_log_b);
  }
  return scalar_shift;
}

inline void finalize_uot_cost_affine_inplace_f(
    arma::fmat& uot_cost,
    const arma::fmat& M,
    const arma::fvec& A,
    const arma::fvec& B,
    float scalar_shift) {
  const arma::uword n = uot_cost.n_rows;
  const arma::uword m = uot_cost.n_cols;
  const float* mp = M.memptr();
  const float* ap = A.memptr();
  const float* bp = B.memptr();
  float* up = uot_cost.memptr();

  for (arma::uword j = 0; j < m; ++j) {
    const float bj = bp[j] + scalar_shift;
    const arma::uword col_offset = j * n;
#ifdef _OPENMP
#pragma omp simd
#endif
    for (arma::uword i = 0; i < n; ++i) {
      const arma::uword idx = col_offset + i;
      up[idx] = -2.0f * up[idx] + mp[idx] + ap[i] + bj;
    }
  }
}

inline void uot_cost_matrix_kl_joint_inplace_f(
    const arma::fmat& X_sqr,
    const arma::fmat& Y_sqr,
    const arma::fmat& X,
    const arma::fmat& Yt,
    const arma::fmat& M,
    const arma::fmat& pi,
    const arma::fvec& a,
    const arma::fvec& b,
    const arma::fvec& log_a,
    const arma::fvec& log_b,
    float rho_x,
    float rho_y,
    float eps,
    arma::fmat& uot_cost,
    arma::fmat& scratch,
    arma::fvec& pi1,
    arma::fvec& pi2,
    arma::fvec& A,
    arma::fvec& B) {
  pi1 = arma::sum(pi, 1);
  pi2 = arma::sum(pi, 0).t();
  A = X_sqr * pi1;
  B = Y_sqr * pi2;
  sgemm_nn(X, pi, scratch);
  sgemm_nn(scratch, Yt, uot_cost);

  (void)a;
  (void)b;
  const float scalar_shift = kl_scalar_shift_from_marginals_and_plan_f(
    pi, pi1, pi2, log_a, log_b, rho_x, rho_y, eps
  );

  finalize_uot_cost_affine_inplace_f(uot_cost, M, A, B, scalar_shift);
}

struct SinkhornUnbalancedResult {
  arma::mat plan;
  int iters;
  double err;
  bool warm_started;
  bool warm_fallback;
};

struct SinkhornUnbalancedResultF {
  arma::fmat plan;
  int iters;
  float err;
  bool warm_started;
  bool warm_fallback;
};

struct SinkhornUnbalancedWorkspace {
  arma::mat K;
  arma::vec u;
  arma::vec v;
  arma::vec Kv;
  arma::vec Ktu;
  arma::vec K1;
  arma::vec KT1;
  arma::vec ones_a;
  arma::vec ones_b;
  bool has_scaling;
  double last_err;
  int last_iters;
  int warm_rejects;

  SinkhornUnbalancedWorkspace()
      : has_scaling(false),
        last_err(std::numeric_limits<double>::infinity()),
        last_iters(0),
        warm_rejects(0) {}
};

struct SinkhornUnbalancedWorkspaceF {
  arma::fmat K;
  arma::fvec u;
  arma::fvec v;
  arma::fvec Kv;
  arma::fvec Ktu;
  arma::fvec K1;
  arma::fvec KT1;
  arma::fvec ones_a;
  arma::fvec ones_b;
  bool has_scaling;
  float last_err;
  int last_iters;
  int warm_rejects;

  SinkhornUnbalancedWorkspaceF()
      : has_scaling(false),
        last_err(std::numeric_limits<float>::infinity()),
        last_iters(0),
        warm_rejects(0) {}
};

inline void ensure_sinkhorn_workspace(
    const arma::mat& cost,
    const arma::vec& a,
    const arma::vec& b,
    SinkhornUnbalancedWorkspace& ws) {
  const bool same_dims =
    ws.K.n_rows == cost.n_rows && ws.K.n_cols == cost.n_cols &&
    ws.u.n_elem == a.n_elem && ws.v.n_elem == b.n_elem;

  if (!same_dims) {
    ws.has_scaling = false;
    ws.last_err = std::numeric_limits<double>::infinity();
    ws.last_iters = 0;
    ws.warm_rejects = 0;
    ws.K.set_size(cost.n_rows, cost.n_cols);
    ws.u.ones(a.n_elem);
    ws.v.ones(b.n_elem);
    ws.Kv.set_size(a.n_elem);
    ws.Ktu.set_size(b.n_elem);
    ws.K1.set_size(a.n_elem);
    ws.KT1.set_size(b.n_elem);
    ws.ones_a.ones(a.n_elem);
    ws.ones_b.ones(b.n_elem);
  } else {
    if (ws.ones_a.n_elem != a.n_elem) {
      ws.ones_a.ones(a.n_elem);
    }
    if (ws.ones_b.n_elem != b.n_elem) {
      ws.ones_b.ones(b.n_elem);
    }
  }
}

inline void ensure_sinkhorn_workspace_f(
    const arma::fmat& cost,
    const arma::fvec& a,
    const arma::fvec& b,
    SinkhornUnbalancedWorkspaceF& ws) {
  const bool same_dims =
    ws.K.n_rows == cost.n_rows && ws.K.n_cols == cost.n_cols &&
    ws.u.n_elem == a.n_elem && ws.v.n_elem == b.n_elem;

  if (!same_dims) {
    ws.has_scaling = false;
    ws.last_err = std::numeric_limits<float>::infinity();
    ws.last_iters = 0;
    ws.warm_rejects = 0;
    ws.K.set_size(cost.n_rows, cost.n_cols);
    ws.u.ones(a.n_elem);
    ws.v.ones(b.n_elem);
    ws.Kv.set_size(a.n_elem);
    ws.Ktu.set_size(b.n_elem);
    ws.K1.set_size(a.n_elem);
    ws.KT1.set_size(b.n_elem);
    ws.ones_a.ones(a.n_elem);
    ws.ones_b.ones(b.n_elem);
  } else {
    if (ws.ones_a.n_elem != a.n_elem) {
      ws.ones_a.ones(a.n_elem);
    }
    if (ws.ones_b.n_elem != b.n_elem) {
      ws.ones_b.ones(b.n_elem);
    }
  }
}

inline void update_unbalanced_scaling(
    const arma::vec& numer,
    const arma::vec& denom,
    double tau,
    bool do_check,
    arma::vec& out,
    double& max_diff) {
  const arma::uword n = out.n_elem;
  const double* np = numer.memptr();
  const double* dp = denom.memptr();
  double* op = out.memptr();
  const bool tau_one = std::abs(tau - 1.0) <= 1e-14;
  const bool tau_zero = std::abs(tau) <= 1e-14;
  max_diff = 0.0;

  if (!do_check) {
    if (tau_one) {
      out = numer / denom;
      return;
    }
    if (tau_zero) {
      out.ones();
      return;
    }
    out = arma::pow(numer / denom, tau);
    return;
  }

  for (arma::uword i = 0; i < n; ++i) {
    const double ratio = np[i] / dp[i];
    const double old_val = op[i];
    const double new_val = tau_one ? ratio : (tau_zero ? 1.0 : std::pow(ratio, tau));
    op[i] = new_val;
    if (do_check) {
      const double d = std::abs(new_val - old_val);
      if (d > max_diff) {
        max_diff = d;
      }
    }
  }
}

inline void update_unbalanced_scaling_f(
    const arma::fvec& numer,
    const arma::fvec& denom,
    float tau,
    bool do_check,
    arma::fvec& out,
    float& max_diff) {
  const arma::uword n = out.n_elem;
  const float* np = numer.memptr();
  const float* dp = denom.memptr();
  float* op = out.memptr();
  const bool tau_one = std::abs(tau - 1.0f) <= 1e-6f;
  const bool tau_zero = std::abs(tau) <= 1e-6f;
  max_diff = 0.0f;

  if (!do_check) {
    if (tau_one) {
      out = numer / denom;
      return;
    }
    if (tau_zero) {
      out.ones();
      return;
    }
    out = arma::pow(numer / denom, tau);
    return;
  }

  for (arma::uword i = 0; i < n; ++i) {
    const float ratio = np[i] / dp[i];
    const float old_val = op[i];
    const float new_val = tau_one ? ratio : (tau_zero ? 1.0f : std::pow(ratio, tau));
    op[i] = new_val;
    if (do_check) {
      const float d = std::abs(new_val - old_val);
      if (d > max_diff) {
        max_diff = d;
      }
    }
  }
}

inline double scaling_update_residual(
    const arma::vec& numer,
    const arma::vec& denom,
    double tau,
    const arma::vec& current) {
  const arma::uword n = current.n_elem;
  const double* np = numer.memptr();
  const double* dp = denom.memptr();
  const double* cp = current.memptr();
  const bool tau_one = std::abs(tau - 1.0) <= 1e-14;
  const bool tau_zero = std::abs(tau) <= 1e-14;

  double max_diff = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double ratio = np[i] / dp[i];
    const double next = tau_one ? ratio : (tau_zero ? 1.0 : std::pow(ratio, tau));
    const double d = std::abs(next - cp[i]);
    if (d > max_diff) {
      max_diff = d;
    }
  }
  return max_diff;
}

inline float scaling_update_residual_f(
    const arma::fvec& numer,
    const arma::fvec& denom,
    float tau,
    const arma::fvec& current) {
  const arma::uword n = current.n_elem;
  const float* np = numer.memptr();
  const float* dp = denom.memptr();
  const float* cp = current.memptr();
  const bool tau_one = std::abs(tau - 1.0f) <= 1e-6f;
  const bool tau_zero = std::abs(tau) <= 1e-6f;

  float max_diff = 0.0f;
  for (arma::uword i = 0; i < n; ++i) {
    const float ratio = np[i] / dp[i];
    const float next = tau_one ? ratio : (tau_zero ? 1.0f : std::pow(ratio, tau));
    const float d = std::abs(next - cp[i]);
    if (d > max_diff) {
      max_diff = d;
    }
  }
  return max_diff;
}

inline double max_abs_vec(const arma::vec& x) {
  const double* xp = x.memptr();
  const arma::uword n = x.n_elem;
  double vmax = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double ax = std::abs(xp[i]);
    if (ax > vmax) {
      vmax = ax;
    }
  }
  return vmax;
}

inline float max_abs_vec_f(const arma::fvec& x) {
  const float* xp = x.memptr();
  const arma::uword n = x.n_elem;
  float vmax = 0.0f;
  for (arma::uword i = 0; i < n; ++i) {
    const float ax = std::abs(xp[i]);
    if (ax > vmax) {
      vmax = ax;
    }
  }
  return vmax;
}

inline double max_abs_diff_vec(const arma::vec& x, const arma::vec& y) {
  const double* xp = x.memptr();
  const double* yp = y.memptr();
  const arma::uword n = x.n_elem;
  double vmax = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double ad = std::abs(xp[i] - yp[i]);
    if (ad > vmax) {
      vmax = ad;
    }
  }
  return vmax;
}

inline float max_abs_diff_vec_f(const arma::fvec& x, const arma::fvec& y) {
  const float* xp = x.memptr();
  const float* yp = y.memptr();
  const arma::uword n = x.n_elem;
  float vmax = 0.0f;
  for (arma::uword i = 0; i < n; ++i) {
    const float ad = std::abs(xp[i] - yp[i]);
    if (ad > vmax) {
      vmax = ad;
    }
  }
  return vmax;
}

inline SinkhornUnbalancedResult sinkhorn_unbalanced_kl(
    const arma::mat& cost,
    const arma::vec& a,
    const arma::vec& b,
    const arma::mat& c,
    double rho1,
    double rho2,
    double eps,
    int max_iter,
    double tol,
    const arma::mat& init_plan,
    SinkhornUnbalancedWorkspace& ws,
    bool allow_warm_start = true) {
  ensure_sinkhorn_workspace(cost, a, b, ws);
  build_kernel_from_cost_mul(cost, c, 1.0 / eps, ws.K);
  const bool use_blas = (ws.K.n_elem >= gemv_min_work_d());
  const bool use_blocked = (!use_blas && ws.K.n_elem >= matvec_blocked_min_work_d());
  const double tau1 = std::isinf(rho1) ? 1.0 : (rho1 <= 0.0 ? 0.0 : rho1 / (rho1 + eps));
  const double tau2 = std::isinf(rho2) ? 1.0 : (rho2 <= 0.0 ? 0.0 : rho2 / (rho2 + eps));

  bool warm_started = allow_warm_start && ws.has_scaling;
  bool warm_fallback = false;
  bool use_warm = warm_started;

  if (use_warm) {
    if (use_blas) {
      dgemv_n(ws.K, ws.v, ws.Kv);
      dgemv_t(ws.K, ws.u, ws.Ktu);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.v, ws.Kv);
      tmatvec_colmajor_blocked(ws.K, ws.u, ws.Ktu);
    } else {
      matvec_colmajor(ws.K, ws.v, ws.Kv);
      tmatvec_colmajor(ws.K, ws.u, ws.Ktu);
    }
    ws.Kv += kTiny;
    ws.Ktu += kTiny;
    const double warm_ru = scaling_update_residual(a, ws.Kv, tau1, ws.u);
    const double warm_rv = scaling_update_residual(b, ws.Ktu, tau2, ws.v);
    const double warm_resid = std::max(warm_ru, warm_rv);
    const double baseline = std::isfinite(ws.last_err) ? std::max(ws.last_err, tol) : tol;
    const double guard = std::max(32.0 * tol, 8.0 * baseline);
    if (!std::isfinite(warm_resid) || warm_resid > guard) {
      use_warm = false;
      warm_fallback = true;
      ++ws.warm_rejects;
    }
  }

  if (!use_warm) {
    ws.u.ones();
    ws.v.ones();
  }

  if (!use_warm && init_plan.n_rows == cost.n_rows && init_plan.n_cols == cost.n_cols) {
    const arma::vec init_r = arma::sum(init_plan, 1) + kTiny;
    const arma::vec init_c = arma::sum(init_plan, 0).t() + kTiny;
    if (use_blas) {
      dgemv_n(ws.K, ws.ones_b, ws.K1);
      dgemv_t(ws.K, ws.ones_a, ws.KT1);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.ones_b, ws.K1);
      tmatvec_colmajor_blocked(ws.K, ws.ones_a, ws.KT1);
    } else {
      matvec_colmajor(ws.K, ws.ones_b, ws.K1);
      tmatvec_colmajor(ws.K, ws.ones_a, ws.KT1);
    }
    ws.u = arma::sqrt(init_r / (ws.K1 + kTiny));
    ws.v = arma::sqrt(init_c / (ws.KT1 + kTiny));
  }

  double err = std::numeric_limits<double>::infinity();
  arma::vec col_prev;
  bool have_col_prev = false;
  const double max_b = std::max(max_abs_vec(b), kTiny);
  const double target_abs_mult = fugw_sinkhorn_target_tol_abs_mult();
  const double target_mult = fugw_sinkhorn_target_tol_mult();
  const double target_cap = fugw_sinkhorn_target_tol_cap_d();
  const double rel_mult = fugw_sinkhorn_rel_tol_mult();
  const double rel_floor = fugw_sinkhorn_rel_tol_floor_d();
  const double col_rel_mult = fugw_sinkhorn_col_rel_tol_mult();
  const double col_rel_floor = fugw_sinkhorn_col_rel_tol_floor_d();
  const double target_tol = std::max(
    target_abs_mult * tol,
    std::min(target_cap, target_mult * std::sqrt(std::max(tol, 1e-16)) * max_b)
  );
  const double rel_tol = std::max(rel_mult * tol, rel_floor);
  const double col_rel_tol = std::max(col_rel_mult * tol, col_rel_floor);
  int it = 0;
  for (; it < max_iter; ++it) {
    const bool do_check = ((it + 1) % 5 == 0);

    if (use_blas) {
      dgemv_n(ws.K, ws.v, ws.Kv);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.v, ws.Kv);
    } else {
      matvec_colmajor(ws.K, ws.v, ws.Kv);
    }
    ws.Kv += kTiny;
    double eu = 0.0;
    update_unbalanced_scaling(a, ws.Kv, tau1, do_check, ws.u, eu);

    if (use_blas) {
      dgemv_t(ws.K, ws.u, ws.Ktu);
    } else if (use_blocked) {
      tmatvec_colmajor_blocked(ws.K, ws.u, ws.Ktu);
    } else {
      tmatvec_colmajor(ws.K, ws.u, ws.Ktu);
    }
    ws.Ktu += kTiny;
    double ev = 0.0;
    update_unbalanced_scaling(b, ws.Ktu, tau2, do_check, ws.v, ev);

    if (do_check) {
      const arma::vec col_cur = ws.v % ws.Ktu;
      const double target_abs = max_abs_diff_vec(col_cur, b);
      const double update_abs = std::max(eu, ev);
      const double max_uv = std::max(max_abs_vec(ws.u), max_abs_vec(ws.v));
      const double update_rel = update_abs / (max_uv + kTiny);
      if (!have_col_prev) {
        col_prev = col_cur;
        have_col_prev = true;
        err = std::max(eu, ev);
      } else {
        const double delta_abs = max_abs_diff_vec(col_cur, col_prev);
        err = delta_abs;
        col_prev = col_cur;
        const double col_scale = std::max(std::max(max_abs_vec(col_cur), max_b), kTiny);
        const double delta_rel = delta_abs / col_scale;
        const bool legacy_stop = (delta_abs < tol);
        const bool robust_stop = (target_abs < target_tol) &&
          (update_rel < rel_tol || delta_rel < col_rel_tol);
        if (legacy_stop || robust_stop) {
          ++it;
          break;
        }
      }
    }
  }

  ws.has_scaling = true;
  ws.last_err = err;
  ws.last_iters = it;
  SinkhornUnbalancedResult out;
  out.plan = (ws.u * ws.v.t()) % ws.K;
  out.iters = it;
  out.err = err;
  out.warm_started = use_warm;
  out.warm_fallback = warm_fallback;
  return out;
}

inline SinkhornUnbalancedResultF sinkhorn_unbalanced_kl_f(
    const arma::fmat& cost,
    const arma::fvec& a,
    const arma::fvec& b,
    const arma::fmat& c,
    float rho1,
    float rho2,
    float eps,
    int max_iter,
    float tol,
    const arma::fmat& init_plan,
    SinkhornUnbalancedWorkspaceF& ws,
    bool allow_warm_start = true) {
  ensure_sinkhorn_workspace_f(cost, a, b, ws);
  build_kernel_from_cost_mul_f(cost, c, 1.0f / eps, ws.K);
  const bool use_blas = (ws.K.n_elem >= gemv_min_work_f());
  const bool use_blocked = (!use_blas && ws.K.n_elem >= matvec_blocked_min_work_f());
  const float tau1 = std::isinf(rho1) ? 1.0f : (rho1 <= 0.0f ? 0.0f : rho1 / (rho1 + eps));
  const float tau2 = std::isinf(rho2) ? 1.0f : (rho2 <= 0.0f ? 0.0f : rho2 / (rho2 + eps));

  bool warm_started = allow_warm_start && ws.has_scaling;
  bool warm_fallback = false;
  bool use_warm = warm_started;

  if (use_warm) {
    if (use_blas) {
      sgemv_n(ws.K, ws.v, ws.Kv);
      sgemv_t(ws.K, ws.u, ws.Ktu);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.v, ws.Kv);
      tmatvec_colmajor_blocked(ws.K, ws.u, ws.Ktu);
    } else {
      matvec_colmajor(ws.K, ws.v, ws.Kv);
      tmatvec_colmajor(ws.K, ws.u, ws.Ktu);
    }
    ws.Kv += kTinyF;
    ws.Ktu += kTinyF;
    const float warm_ru = scaling_update_residual_f(a, ws.Kv, tau1, ws.u);
    const float warm_rv = scaling_update_residual_f(b, ws.Ktu, tau2, ws.v);
    const float warm_resid = std::max(warm_ru, warm_rv);
    const float baseline = std::isfinite(ws.last_err) ? std::max(ws.last_err, tol) : tol;
    const float guard = std::max(32.0f * tol, 8.0f * baseline);
    if (!std::isfinite(warm_resid) || warm_resid > guard) {
      use_warm = false;
      warm_fallback = true;
      ++ws.warm_rejects;
    }
  }

  if (!use_warm) {
    ws.u.ones();
    ws.v.ones();
  }

  if (!use_warm && init_plan.n_rows == cost.n_rows && init_plan.n_cols == cost.n_cols) {
    const arma::fvec init_r = arma::sum(init_plan, 1) + kTinyF;
    const arma::fvec init_c = arma::sum(init_plan, 0).t() + kTinyF;
    if (use_blas) {
      sgemv_n(ws.K, ws.ones_b, ws.K1);
      sgemv_t(ws.K, ws.ones_a, ws.KT1);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.ones_b, ws.K1);
      tmatvec_colmajor_blocked(ws.K, ws.ones_a, ws.KT1);
    } else {
      matvec_colmajor(ws.K, ws.ones_b, ws.K1);
      tmatvec_colmajor(ws.K, ws.ones_a, ws.KT1);
    }
    ws.u = arma::sqrt(init_r / (ws.K1 + kTinyF));
    ws.v = arma::sqrt(init_c / (ws.KT1 + kTinyF));
  }

  float err = std::numeric_limits<float>::infinity();
  arma::fvec col_prev;
  bool have_col_prev = false;
  const float max_b = std::max(max_abs_vec_f(b), kTinyF);
  const float target_abs_mult = static_cast<float>(fugw_sinkhorn_target_tol_abs_mult());
  const float target_mult = static_cast<float>(fugw_sinkhorn_target_tol_mult());
  const float target_cap = static_cast<float>(fugw_sinkhorn_target_tol_cap_f());
  const float rel_mult = static_cast<float>(fugw_sinkhorn_rel_tol_mult());
  const float rel_floor = static_cast<float>(fugw_sinkhorn_rel_tol_floor_f());
  const float col_rel_mult = static_cast<float>(fugw_sinkhorn_col_rel_tol_mult());
  const float col_rel_floor = static_cast<float>(fugw_sinkhorn_col_rel_tol_floor_f());
  const float target_tol = std::max(
    target_abs_mult * tol,
    std::min(target_cap, target_mult * std::sqrt(std::max(tol, 1e-12f)) * max_b)
  );
  const float rel_tol = std::max(rel_mult * tol, rel_floor);
  const float col_rel_tol = std::max(col_rel_mult * tol, col_rel_floor);
  int it = 0;
  for (; it < max_iter; ++it) {
    const bool do_check = ((it + 1) % 5 == 0);

    if (use_blas) {
      sgemv_n(ws.K, ws.v, ws.Kv);
    } else if (use_blocked) {
      matvec_colmajor_blocked(ws.K, ws.v, ws.Kv);
    } else {
      matvec_colmajor(ws.K, ws.v, ws.Kv);
    }
    ws.Kv += kTinyF;
    float eu = 0.0f;
    update_unbalanced_scaling_f(a, ws.Kv, tau1, do_check, ws.u, eu);

    if (use_blas) {
      sgemv_t(ws.K, ws.u, ws.Ktu);
    } else if (use_blocked) {
      tmatvec_colmajor_blocked(ws.K, ws.u, ws.Ktu);
    } else {
      tmatvec_colmajor(ws.K, ws.u, ws.Ktu);
    }
    ws.Ktu += kTinyF;
    float ev = 0.0f;
    update_unbalanced_scaling_f(b, ws.Ktu, tau2, do_check, ws.v, ev);

    if (do_check) {
      const arma::fvec col_cur = ws.v % ws.Ktu;
      const float target_abs = max_abs_diff_vec_f(col_cur, b);
      const float update_abs = std::max(eu, ev);
      const float max_uv = std::max(max_abs_vec_f(ws.u), max_abs_vec_f(ws.v));
      const float update_rel = update_abs / (max_uv + kTinyF);
      if (!have_col_prev) {
        col_prev = col_cur;
        have_col_prev = true;
        err = std::max(eu, ev);
      } else {
        const float delta_abs = max_abs_diff_vec_f(col_cur, col_prev);
        err = delta_abs;
        col_prev = col_cur;
        const float col_scale = std::max(std::max(max_abs_vec_f(col_cur), max_b), kTinyF);
        const float delta_rel = delta_abs / col_scale;
        const bool legacy_stop = (delta_abs < tol);
        const bool robust_stop = (target_abs < target_tol) &&
          (update_rel < rel_tol || delta_rel < col_rel_tol);
        if (legacy_stop || robust_stop) {
          ++it;
          break;
        }
      }
    }
  }

  ws.has_scaling = true;
  ws.last_err = err;
  ws.last_iters = it;
  SinkhornUnbalancedResultF out;
  out.plan = (ws.u * ws.v.t()) % ws.K;
  out.iters = it;
  out.err = err;
  out.warm_started = use_warm;
  out.warm_fallback = warm_fallback;
  return out;
}

inline std::pair<double, double> fused_unbalanced_cost_square_joint(
    const arma::mat& Cx,
    const arma::mat& Cy,
    const arma::mat& Cx_sqr,
    const arma::mat& Cy_sqr,
    const arma::mat& M_samp,
    const arma::mat& M_feat,
    const arma::vec& wx,
    const arma::vec& wy,
    const arma::mat& wxy,
    const arma::mat& pi_samp,
    const arma::mat& pi_feat,
    double rho_x,
    double rho_y,
    double eps) {
  const arma::vec pi1_samp = arma::sum(pi_samp, 1);
  const arma::vec pi2_samp = arma::sum(pi_samp, 0).t();
  const arma::vec pi1_feat = arma::sum(pi_feat, 1);
  const arma::vec pi2_feat = arma::sum(pi_feat, 0).t();

  const double A_sqr = arma::dot(Cx_sqr * pi1_feat, pi1_samp);
  const double B_sqr = arma::dot(Cy_sqr * pi2_feat, pi2_samp);
  const arma::mat AB = (Cx * pi_feat * Cy.t()) % pi_samp;
  const double linear_cost = A_sqr + B_sqr - 2.0 * arma::accu(AB);

  double ucoot_cost = linear_cost;
  ucoot_cost += arma::accu(pi_samp % M_samp);
  ucoot_cost += arma::accu(pi_feat % M_feat);

  if (std::isfinite(rho_x) && rho_x != 0.0) {
    ucoot_cost += rho_x * div_between_product_kl(pi1_samp, pi1_feat, wx, wx);
  }
  if (std::isfinite(rho_y) && rho_y != 0.0) {
    ucoot_cost += rho_y * div_between_product_kl(pi2_samp, pi2_feat, wy, wy);
  }
  if (eps != 0.0) {
    ucoot_cost += eps * div_between_product_kl(pi_samp, pi_feat, wxy, wxy);
  }

  return {linear_cost, ucoot_cost};
}

inline std::size_t basis_index(int i, int j, int m) {
  return static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + static_cast<std::size_t>(j);
}

inline void init_transport_northwest(
    const arma::vec& p,
    const arma::vec& q,
    arma::mat& x,
    std::vector<unsigned char>& basis,
    double eps) {
  const int n = static_cast<int>(p.n_elem);
  const int m = static_cast<int>(q.n_elem);
  x.zeros(n, m);
  basis.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(m), 0);

  arma::vec supply = p;
  arma::vec demand = q;
  int i = 0;
  int j = 0;
  int nbas = 0;

  while (i < n && j < m) {
    const double val = std::min(supply(i), demand(j));
    x(i, j) = val;
    const std::size_t ij = basis_index(i, j, m);
    if (!basis[ij]) {
      basis[ij] = 1;
      ++nbas;
    }
    supply(i) -= val;
    demand(j) -= val;

    const bool row_done = (supply(i) <= eps);
    const bool col_done = (demand(j) <= eps);
    if (row_done) {
      supply(i) = 0.0;
    }
    if (col_done) {
      demand(j) = 0.0;
    }

    if (row_done && col_done) {
      if (i + 1 < n && nbas < n + m - 1) {
        const std::size_t zi = basis_index(i + 1, j, m);
        if (!basis[zi]) {
          basis[zi] = 1;
          ++nbas;
        }
      } else if (j + 1 < m && nbas < n + m - 1) {
        const std::size_t zj = basis_index(i, j + 1, m);
        if (!basis[zj]) {
          basis[zj] = 1;
          ++nbas;
        }
      }
      ++i;
      ++j;
    } else if (row_done) {
      ++i;
    } else if (col_done) {
      ++j;
    } else {
      if (supply(i) < demand(j)) {
        ++i;
      } else {
        ++j;
      }
    }
  }
}

inline void compute_transport_potentials(
    const arma::mat& cost,
    const std::vector<unsigned char>& basis,
    int n,
    int m,
    std::vector<double>& u,
    std::vector<double>& v) {
  u.assign(n, std::numeric_limits<double>::quiet_NaN());
  v.assign(m, std::numeric_limits<double>::quiet_NaN());
  std::deque<int> dq;

  for (int seed = 0; seed < n; ++seed) {
    if (std::isfinite(u[seed])) {
      continue;
    }
    u[seed] = 0.0;
    dq.push_back(seed);

    while (!dq.empty()) {
      const int node = dq.front();
      dq.pop_front();
      if (node < n) {
        const int i = node;
        for (int j = 0; j < m; ++j) {
          if (!basis[basis_index(i, j, m)]) {
            continue;
          }
          if (!std::isfinite(v[j])) {
            v[j] = cost(i, j) - u[i];
            dq.push_back(n + j);
          }
        }
      } else {
        const int j = node - n;
        for (int i = 0; i < n; ++i) {
          if (!basis[basis_index(i, j, m)]) {
            continue;
          }
          if (!std::isfinite(u[i])) {
            u[i] = cost(i, j) - v[j];
            dq.push_back(i);
          }
        }
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(u[i])) {
      u[i] = 0.0;
    }
  }
  for (int j = 0; j < m; ++j) {
    if (!std::isfinite(v[j])) {
      v[j] = 0.0;
    }
  }
}

inline bool find_basis_path(
    const std::vector<unsigned char>& basis,
    int n,
    int m,
    int row_start,
    int col_target,
    std::vector<int>& parent) {
  const int target = n + col_target;
  parent.assign(n + m, -1);
  std::queue<int> q;
  parent[row_start] = row_start;
  q.push(row_start);

  while (!q.empty() && parent[target] == -1) {
    const int node = q.front();
    q.pop();
    if (node < n) {
      const int i = node;
      for (int j = 0; j < m; ++j) {
        if (!basis[basis_index(i, j, m)]) {
          continue;
        }
        const int nxt = n + j;
        if (parent[nxt] == -1) {
          parent[nxt] = node;
          q.push(nxt);
        }
      }
    } else {
      const int j = node - n;
      for (int i = 0; i < n; ++i) {
        if (!basis[basis_index(i, j, m)]) {
          continue;
        }
        const int nxt = i;
        if (parent[nxt] == -1) {
          parent[nxt] = node;
          q.push(nxt);
        }
      }
    }
  }

  return parent[target] != -1;
}

struct TransportSimplexResult {
  arma::mat plan;
  int iterations;
  bool converged;
};

inline double solve_1d_linesearch_quad(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) {
    return 0.0;
  }
  if (a > 0.0) {
    return -b / (2.0 * a);
  }
  if ((a + b) < 0.0) {
    return 1.0;
  }
  return 0.0;
}

inline bool is_uniform_prob(const arma::vec& w, double tol = 1e-12) {
  if (w.n_elem == 0) {
    return false;
  }
  const double target = 1.0 / static_cast<double>(w.n_elem);
  return (std::abs(arma::accu(w) - 1.0) <= tol) && (arma::max(arma::abs(w - target)) <= tol);
}

inline arma::uvec weighted_sample_replace(const arma::vec& prob, int k) {
  const arma::uword n = prob.n_elem;
  arma::uvec out(static_cast<arma::uword>(std::max(0, k)), arma::fill::zeros);
  if (n == 0 || k <= 0) {
    return out;
  }

  arma::vec cdf(n, arma::fill::zeros);
  double accum = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double w = std::max(0.0, prob[i]);
    accum += w;
    cdf[i] = accum;
  }
  if (!(accum > 0.0) || !std::isfinite(accum)) {
    for (int t = 0; t < k; ++t) {
      out[static_cast<arma::uword>(t)] = static_cast<arma::uword>(std::floor(R::runif(0.0, static_cast<double>(n))));
    }
    return out;
  }
  cdf /= accum;

  for (int t = 0; t < k; ++t) {
    const double u = std::min(1.0 - 1e-15, std::max(0.0, R::runif(0.0, 1.0)));
    arma::uword lo = 0;
    arma::uword hi = n - 1;
    while (lo < hi) {
      const arma::uword mid = lo + (hi - lo) / 2;
      if (cdf[mid] >= u) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    out[static_cast<arma::uword>(t)] = lo;
  }
  return out;
}

inline arma::uvec weighted_sample_no_replace(const arma::vec& prob, int k) {
  const arma::uword n = prob.n_elem;
  const int kk = std::max(0, std::min(k, static_cast<int>(n)));
  arma::uvec out(static_cast<arma::uword>(kk), arma::fill::zeros);
  if (n == 0 || kk == 0) {
    return out;
  }
  if (kk == static_cast<int>(n)) {
    for (arma::uword i = 0; i < n; ++i) {
      out[i] = i;
    }
    return out;
  }

  arma::vec w = prob;
  w.transform([](double v) { return (std::isfinite(v) && v > 0.0) ? v : 0.0; });
  double total = arma::accu(w);
  std::vector<unsigned char> picked(static_cast<std::size_t>(n), 0);
  int filled = 0;

  for (int d = 0; d < kk; ++d) {
    if (!(total > 0.0) || !std::isfinite(total)) {
      break;
    }
    const double u = R::runif(0.0, total);
    double csum = 0.0;
    arma::uword chosen = n - 1;
    for (arma::uword i = 0; i < n; ++i) {
      const double wi = w[i];
      if (wi <= 0.0) {
        continue;
      }
      csum += wi;
      if (csum >= u) {
        chosen = i;
        break;
      }
    }
    out[static_cast<arma::uword>(filled)] = chosen;
    ++filled;
    picked[static_cast<std::size_t>(chosen)] = 1;
    total -= w[chosen];
    w[chosen] = 0.0;
  }

  if (filled < kk) {
    for (arma::uword i = 0; i < n && filled < kk; ++i) {
      if (!picked[static_cast<std::size_t>(i)]) {
        out[static_cast<arma::uword>(filled)] = i;
        ++filled;
      }
    }
  }
  return out;
}

inline arma::uvec deterministic_topk_indices(const arma::vec& prob, int k) {
  const arma::uword n = prob.n_elem;
  const int kk = std::max(0, std::min(k, static_cast<int>(n)));
  arma::uvec out(static_cast<arma::uword>(kk), arma::fill::zeros);
  if (n == 0 || kk == 0) {
    return out;
  }
  std::vector<arma::uword> idx(n);
  std::iota(idx.begin(), idx.end(), static_cast<arma::uword>(0));
  std::stable_sort(
    idx.begin(),
    idx.end(),
    [&](arma::uword a, arma::uword b) {
      const double wa = (std::isfinite(prob[a]) && prob[a] > 0.0) ? prob[a] : 0.0;
      const double wb = (std::isfinite(prob[b]) && prob[b] > 0.0) ? prob[b] : 0.0;
      if (wa == wb) {
        return a < b;
      }
      return wa > wb;
    }
  );
  for (int i = 0; i < kk; ++i) {
    out[static_cast<arma::uword>(i)] = idx[static_cast<std::size_t>(i)];
  }
  return out;
}

inline arma::uvec hungarian_assignment(const arma::mat& cost) {
  const std::size_t n = cost.n_rows;
  arma::uvec assignment(n, arma::fill::zeros);
  std::vector<double> u(n + 1u, 0.0);
  std::vector<double> v(n + 1u, 0.0);
  std::vector<int> p(n + 1u, 0);
  std::vector<int> way(n + 1u, 0);

  for (std::size_t i = 1; i <= n; ++i) {
    p[0] = static_cast<int>(i);
    int j0 = 0;
    std::vector<double> minv(n + 1u, std::numeric_limits<double>::infinity());
    std::vector<unsigned char> used(n + 1u, 0);

    do {
      used[static_cast<std::size_t>(j0)] = 1;
      const int i0 = p[static_cast<std::size_t>(j0)];
      double delta = std::numeric_limits<double>::infinity();
      int j1 = 0;

      for (std::size_t j = 1; j <= n; ++j) {
        if (used[j]) {
          continue;
        }
        const double cur = cost(static_cast<arma::uword>(i0 - 1), static_cast<arma::uword>(j - 1)) -
          u[static_cast<std::size_t>(i0)] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = static_cast<int>(j);
        }
      }

      for (std::size_t j = 0; j <= n; ++j) {
        if (used[j]) {
          u[static_cast<std::size_t>(p[j])] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[static_cast<std::size_t>(j0)] != 0);

    do {
      const int j1 = way[static_cast<std::size_t>(j0)];
      p[static_cast<std::size_t>(j0)] = p[static_cast<std::size_t>(j1)];
      j0 = j1;
    } while (j0 != 0);
  }

  for (std::size_t j = 1; j <= n; ++j) {
    const int i = p[j];
    if (i > 0) {
      assignment(static_cast<arma::uword>(i - 1)) = static_cast<arma::uword>(j - 1);
    }
  }
  return assignment;
}

inline TransportSimplexResult transport_simplex_solve(
    const arma::mat& cost,
    const arma::vec& p,
    const arma::vec& q,
    int max_iter = 20000,
    double tol_opt = 1e-12,
    double eps = 1e-14) {
  const int n = static_cast<int>(p.n_elem);
  const int m = static_cast<int>(q.n_elem);

  if (n == m && cost.n_rows == cost.n_cols && is_uniform_prob(p, 1e-9) && is_uniform_prob(q, 1e-9)) {
    arma::mat plan(n, n, arma::fill::zeros);
    const arma::uvec assign = hungarian_assignment(cost);
    for (int i = 0; i < n; ++i) {
      plan(static_cast<arma::uword>(i), assign(static_cast<arma::uword>(i))) = p(static_cast<arma::uword>(i));
    }
    TransportSimplexResult out;
    out.plan = std::move(plan);
    out.iterations = 1;
    out.converged = true;
    return out;
  }

  arma::mat x(n, m, arma::fill::zeros);
  std::vector<unsigned char> basis;
  init_transport_northwest(p, q, x, basis, eps);

  std::vector<double> u, v;
  std::vector<int> parent;
  std::vector<int> path_nodes;
  std::vector<std::pair<int, int>> cycle;

  int it = 0;
  for (; it < max_iter; ++it) {
    compute_transport_potentials(cost, basis, n, m, u, v);

    int enter_i = -1;
    int enter_j = -1;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        if (basis[basis_index(i, j, m)]) {
          continue;
        }
        const double rc = cost(i, j) - u[i] - v[j];
        if (rc < -tol_opt) {
          enter_i = i;
          enter_j = j;
          break;
        }
      }
      if (enter_i >= 0) {
        break;
      }
    }

    if (enter_i < 0) {
      break;
    }

    if (!find_basis_path(basis, n, m, enter_i, enter_j, parent)) {
      break;
    }

    path_nodes.clear();
    int cur = n + enter_j;
    while (cur != parent[cur]) {
      path_nodes.push_back(cur);
      cur = parent[cur];
    }
    path_nodes.push_back(enter_i);
    std::reverse(path_nodes.begin(), path_nodes.end());

    cycle.clear();
    cycle.emplace_back(enter_i, enter_j);
    for (std::size_t k = 0; k + 1 < path_nodes.size(); ++k) {
      const int a = path_nodes[k];
      const int b = path_nodes[k + 1];
      if (a < n && b >= n) {
        cycle.emplace_back(a, b - n);
      } else if (a >= n && b < n) {
        cycle.emplace_back(b, a - n);
      }
    }

    double theta = std::numeric_limits<double>::infinity();
    for (std::size_t k = 1; k < cycle.size(); k += 2) {
      const auto ij = cycle[k];
      theta = std::min(theta, x(ij.first, ij.second));
    }
    if (!std::isfinite(theta)) {
      break;
    }

    int leave_i = -1;
    int leave_j = -1;
    std::size_t leave_idx = std::numeric_limits<std::size_t>::max();
    for (std::size_t k = 1; k < cycle.size(); k += 2) {
      const auto ij = cycle[k];
      const double val = x(ij.first, ij.second);
      if (val <= theta + eps) {
        const std::size_t idx = basis_index(ij.first, ij.second, m);
        if (idx < leave_idx) {
          leave_idx = idx;
          leave_i = ij.first;
          leave_j = ij.second;
        }
      }
    }
    if (leave_i < 0) {
      break;
    }

    for (std::size_t k = 0; k < cycle.size(); ++k) {
      const auto ij = cycle[k];
      if ((k % 2) == 0) {
        x(ij.first, ij.second) += theta;
      } else {
        x(ij.first, ij.second) -= theta;
        if (x(ij.first, ij.second) < eps) {
          x(ij.first, ij.second) = 0.0;
        }
      }
    }

    basis[basis_index(enter_i, enter_j, m)] = 1;
    basis[basis_index(leave_i, leave_j, m)] = 0;
  }

  TransportSimplexResult out;
  out.plan = std::move(x);
  out.iterations = it;
  out.converged = (it < max_iter);
  return out;
}

inline Rcpp::List fgw_entropic_square_mixed_impl(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    int check_every,
    int approx_rank,
    const arma::mat& init_plan) {
  const float alpha_f = static_cast<float>(alpha);
  const float epsilon_f = static_cast<float>(epsilon);
  const float sinkhorn_tol_f = static_cast<float>(std::max(sinkhorn_tol, 1e-6));
  const float tol_f = static_cast<float>(std::max(tol, 1e-6));

  arma::fmat Mf = arma::conv_to<arma::fmat>::from(M);
  arma::fmat C1f = arma::conv_to<arma::fmat>::from(C1);
  arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
  arma::fvec pf = arma::conv_to<arma::fvec>::from(p);
  arma::fvec qf = arma::conv_to<arma::fvec>::from(q);

  arma::fmat constC, hC1, hC2;
  init_matrices_square_f(C1f, C2f, pf, qf, constC, hC1, hC2);
  arma::fmat constCt, hC1t, hC2_asym;
  if (!symmetric) {
    init_matrices_square_f(C1f.t(), C2f.t(), pf, qf, constCt, hC1t, hC2_asym);
  }
  const arma::fmat M_lin = (1.0f - alpha_f) * Mf;
  bool use_lowrank = false;
  arma::fmat lr_A1, lr_B1t, lr_A2, lr_A2t, lr_B2t;
  if (symmetric && approx_rank > 0) {
    use_lowrank = svd_lowrank_factors_f(C1f, approx_rank, lr_A1, lr_B1t) &&
      svd_lowrank_factors_f(C2f, approx_rank, lr_A2, lr_B2t);
    if (use_lowrank) {
      lr_A2t = lr_A2.t();
    }
  }
  const float symmetric_cross_scale = use_lowrank ? (-4.0f * alpha_f) : (-2.0f * alpha_f);
  const float asym_cross_scale = -alpha_f;
  arma::fmat base_sym;
  arma::fmat base_asym;
  if (symmetric) {
    base_sym = (2.0f * alpha_f) * constC;
    base_sym += M_lin;
  } else {
    base_asym = alpha_f * (constC + constCt);
    base_asym += M_lin;
  }

  arma::fmat T;
  if (init_plan.n_rows == p.n_elem && init_plan.n_cols == q.n_elem) {
    T = arma::conv_to<arma::fmat>::from(init_plan);
    T.for_each([](float& v) {
      if (!std::isfinite(v) || v < 0.0f) {
        v = 0.0f;
      }
    });
    const float m = arma::accu(T);
    if (m > 0.0f) {
      T /= m;
    } else {
      T = pf * qf.t();
    }
  } else {
    T = pf * qf.t();
  }
  arma::fmat T_prev_check = T;
  arma::fvec u_ws = arma::ones<arma::fvec>(pf.n_elem);
  arma::fvec v_ws = arma::ones<arma::fvec>(qf.n_elem);
  Rcpp::NumericVector err_trace;
  arma::fmat tens(C1f.n_rows, C2f.n_rows);
  arma::fmat scratch(C1f.n_rows, C2f.n_rows);
  arma::fmat scratch_asym;
  if (!symmetric) {
    scratch_asym.set_size(C1f.n_rows, C2f.n_rows);
  }
  arma::fmat lr_tmp_r1_nt;
  arma::fmat lr_tmp_r1_r2;
  arma::fmat lr_tmp_ns_r2;
  if (use_lowrank) {
    lr_tmp_r1_nt.set_size(lr_B1t.n_rows, C2f.n_rows);
    lr_tmp_r1_r2.set_size(lr_B1t.n_rows, lr_B2t.n_rows);
    lr_tmp_ns_r2.set_size(C1f.n_rows, lr_B2t.n_rows);
  }

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if ((it + 1) % check_every == 0) {
      T_prev_check = T;
    }

    if (symmetric) {
      if (use_lowrank) {
        tensor_product_lowrank_scaled_f(
          base_sym, lr_A1, lr_B1t, lr_A2t, lr_B2t, T, symmetric_cross_scale,
          tens, lr_tmp_r1_nt, lr_tmp_r1_r2, lr_tmp_ns_r2
        );
      } else {
        tensor_product_blas_scaled_f(
          base_sym, hC1, hC2, T, symmetric_cross_scale, tens, scratch
        );
      }
    } else {
      tensor_product_asym_blas_scaled_f(
        base_asym, hC1, hC2, hC1t, hC2_asym, T, asym_cross_scale,
        tens, scratch, scratch_asym
      );
    }

    if (use_ppa) {
      tens -= epsilon_f * arma::log(T + kTinyF);
    }

    SinkhornBalancedResultF sk = sinkhorn_balanced_f(
      pf, qf, tens, epsilon_f, sinkhorn_max_iter, sinkhorn_tol_f, u_ws, v_ws
    );
    T = std::move(sk.plan);
    u_ws = std::move(sk.u);
    v_ws = std::move(sk.v);

    if ((it + 1) % check_every == 0) {
      err = static_cast<double>(arma::norm(T - T_prev_check, "fro"));
      err_trace.push_back(err);
      if (err < tol_f) {
        ++it;
        break;
      }
    }
  }

  arma::mat Td = arma::conv_to<arma::mat>::from(T);
  arma::mat constCd, hC1d, hC2d;
  init_matrices_square(C1, C2, p, q, constCd, hC1d, hC2d);
  const double fgw_dist = (1.0 - alpha) * arma::accu(M % Td) + alpha * gw_loss(constCd, hC1d, hC2d, Td);

  return Rcpp::List::create(
    Rcpp::Named("plan") = Td,
    Rcpp::Named("fgw_dist") = fgw_dist,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = err,
    Rcpp::Named("err_trace") = err_trace
  );
}

struct FgwEntropicCoreResult {
  arma::mat plan;
  double fgw_dist;
  int iterations;
  double error;
  bool used_lowrank;
  bool used_c1_cache;
  bool used_c2_cache;
  bool used_square_cache;
};

struct SrfgwEntropicCoreResult {
  arma::mat plan;
  arma::vec q;
  double lin_loss;
  double quad_loss;
  double srfgw_dist;
  double srgw_dist;
  int iterations;
  double error;
};

inline void normalize_rows_to_p(
    arma::mat& G,
    const arma::vec& p,
    arma::uword nt) {
  arma::vec rs = arma::sum(G, 1);
  for (arma::uword i = 0; i < G.n_rows; ++i) {
    if (std::isfinite(rs[i]) && rs[i] > 0.0) {
      G.row(i) *= (p[i] / rs[i]);
    } else {
      G.row(i).fill(p[i] / static_cast<double>(nt));
    }
  }
}

inline SrfgwEntropicCoreResult entropic_semirelaxed_fgw_square_core(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    double epsilon,
    double alpha,
    bool symmetric,
    int max_iter,
    double tol,
    int check_every,
    const arma::mat& init_plan) {
  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  if (check_every <= 0) {
    check_every = 1;
  }

  const arma::mat C1_sqr = C1 % C1;
  const arma::mat C2_sqr = C2 % C2;
  const arma::mat fC2t = C2_sqr.t();
  const arma::vec left = C1_sqr * p;

  arma::mat constC(ns, nt);
  for (arma::uword j = 0; j < nt; ++j) {
    constC.col(j) = left;
  }
  const arma::mat hC1 = C1;
  const arma::mat hC2 = 2.0 * C2;

  arma::mat constCt;
  arma::mat hC1t;
  arma::mat hC2t;
  arma::mat fC2;
  if (!symmetric) {
    const arma::mat C1_t = C1.t();
    const arma::mat C2_t = C2.t();
    const arma::mat C1_t_sqr = C1_t % C1_t;
    const arma::vec left_t = C1_t_sqr * p;
    constCt.set_size(ns, nt);
    for (arma::uword j = 0; j < nt; ++j) {
      constCt.col(j) = left_t;
    }
    hC1t = C1_t;
    hC2t = 2.0 * C2_t;
    fC2 = C2_sqr;
  }

  arma::mat G;
  if (init_plan.n_rows == ns && init_plan.n_cols == nt) {
    G = init_plan;
    G.for_each([](double& v) {
      if (!std::isfinite(v) || v < 0.0) {
        v = 0.0;
      }
    });
    normalize_rows_to_p(G, p, nt);
  } else {
    G = p * (arma::ones<arma::rowvec>(nt) / static_cast<double>(nt));
  }

  arma::mat G_prev = G;
  arma::mat scratch(ns, nt);
  arma::mat scratch2;
  if (!symmetric) {
    scratch2.set_size(ns, nt);
  }
  arma::mat grad_quad(ns, nt);
  arma::mat grad1;
  arma::mat grad2;
  if (!symmetric) {
    grad1.set_size(ns, nt);
    grad2.set_size(ns, nt);
  }
  arma::mat K(ns, nt);
  arma::vec rs(ns);

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  const double inv_eps = 1.0 / epsilon;
  const double alpha_lin = 1.0 - alpha;
  const bool alpha_full = (alpha >= 1.0);
  const bool alpha_zero = (alpha <= 0.0);
  const double quad_scale = symmetric ? (2.0 * alpha) : alpha;
  for (int k = 0; k < max_iter; ++k) {
    const bool do_check = ((k + 1) % check_every == 0) || (k == 0);
    if (do_check) {
      G_prev = G;
    }

    const arma::rowvec qG = arma::sum(G, 0);
    if (symmetric) {
      const arma::rowvec marg = qG * fC2t;
      grad_quad = constC;
      grad_quad.each_row() += marg;
      dgemm_nn(hC1, G, scratch);
      dgemm_nt_accum(scratch, hC2, grad_quad, -1.0, 1.0);
    } else {
      const arma::rowvec marg1 = qG * fC2t;
      grad1 = constC;
      grad1.each_row() += marg1;
      dgemm_nn(hC1, G, scratch);
      dgemm_nt_accum(scratch, hC2, grad1, -1.0, 1.0);
      grad1 *= 2.0;

      const arma::rowvec marg2 = qG * fC2;
      grad2 = constCt;
      grad2.each_row() += marg2;
      dgemm_nn(hC1t, G, scratch2);
      dgemm_nt_accum(scratch2, hC2t, grad2, -1.0, 1.0);
      grad2 *= 2.0;

      grad_quad = 0.5 * (grad1 + grad2);
    }

    if (alpha_full) {
      const double inv_eps_grad = inv_eps * (symmetric ? 2.0 : 1.0);
      build_kernel_from_cost_mul(grad_quad, G, inv_eps_grad, K);
    } else if (alpha_zero) {
      build_kernel_from_cost_mul(M, G, inv_eps, K);
    } else {
      build_kernel_from_fused_cost_mul(grad_quad, M, G, inv_eps, quad_scale, alpha_lin, K);
    }
    rs = arma::sum(K, 1);
    for (arma::uword i = 0; i < ns; ++i) {
      const double rsi = rs[i];
      if (std::isfinite(rsi) && rsi > 0.0) {
        G.row(i) = (p[i] / rsi) * K.row(i);
      } else {
        G.row(i).fill(p[i] / static_cast<double>(nt));
      }
    }

    if (do_check) {
      err = arma::norm(G - G_prev, "fro");
      if (err <= tol) {
        it = k + 1;
        break;
      }
    }
    it = k + 1;
  }

  const arma::rowvec qG = arma::sum(G, 0);
  const arma::rowvec marg = qG * fC2t;
  arma::mat tens = constC;
  tens.each_row() += marg;
  dgemm_nn(hC1, G, scratch);
  dgemm_nt_accum(scratch, hC2, tens, -1.0, 1.0);
  const double quad_raw = arma::accu(tens % G);
  const double lin_loss = (alpha >= 1.0 || M.n_elem == 0) ? 0.0 : (1.0 - alpha) * arma::accu(M % G);
  const double quad_loss = alpha * quad_raw;

  SrfgwEntropicCoreResult out;
  out.plan = std::move(G);
  out.q = qG.t();
  out.lin_loss = lin_loss;
  out.quad_loss = quad_loss;
  out.srfgw_dist = lin_loss + quad_loss;
  out.srgw_dist = quad_raw;
  out.iterations = it;
  out.error = err;
  return out;
}

inline SrfgwEntropicCoreResult entropic_semirelaxed_fgw_square_core_mixed(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    double epsilon,
    double alpha,
    bool symmetric,
    int max_iter,
    double tol,
    int check_every,
    const arma::mat& init_plan) {
  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  if (check_every <= 0) {
    check_every = 1;
  }

  const float eps_f = static_cast<float>(epsilon);
  const float alpha_f = static_cast<float>(alpha);
  const float alpha_lin_f = 1.0f - alpha_f;
  const float inv_eps_f = 1.0f / eps_f;
  const bool alpha_full = (alpha_f >= 1.0f);
  const bool alpha_zero = (alpha_f <= 0.0f);
  const float quad_scale_f = symmetric ? (2.0f * alpha_f) : alpha_f;

  arma::fmat Mf;
  if (!alpha_full) {
    Mf = arma::conv_to<arma::fmat>::from(M);
  }
  const arma::fmat C1f = arma::conv_to<arma::fmat>::from(C1);
  const arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
  const arma::fvec pf = arma::conv_to<arma::fvec>::from(p);

  const arma::fmat C1_sqr = C1f % C1f;
  const arma::fmat C2_sqr = C2f % C2f;
  const arma::fmat fC2t = C2_sqr.t();
  const arma::fvec left = C1_sqr * pf;

  arma::fmat constC(ns, nt);
  for (arma::uword j = 0; j < nt; ++j) {
    constC.col(j) = left;
  }
  const arma::fmat hC1 = C1f;
  const arma::fmat hC2 = 2.0f * C2f;

  arma::fmat constCt;
  arma::fmat hC1t;
  arma::fmat hC2t;
  arma::fmat fC2;
  if (!symmetric) {
    const arma::fmat C1_t = C1f.t();
    const arma::fmat C2_t = C2f.t();
    const arma::fmat C1_t_sqr = C1_t % C1_t;
    const arma::fvec left_t = C1_t_sqr * pf;
    constCt.set_size(ns, nt);
    for (arma::uword j = 0; j < nt; ++j) {
      constCt.col(j) = left_t;
    }
    hC1t = C1_t;
    hC2t = 2.0f * C2_t;
    fC2 = C2_sqr;
  }

  arma::fmat Gf;
  if (init_plan.n_rows == ns && init_plan.n_cols == nt) {
    Gf = arma::conv_to<arma::fmat>::from(init_plan);
    Gf.for_each([](float& v) {
      if (!std::isfinite(v) || v < 0.0f) {
        v = 0.0f;
      }
    });
    arma::fvec rs0 = arma::sum(Gf, 1);
    for (arma::uword i = 0; i < ns; ++i) {
      const float rsi = rs0[i];
      if (std::isfinite(rsi) && rsi > 0.0f) {
        Gf.row(i) *= (pf[i] / rsi);
      } else {
        Gf.row(i).fill(pf[i] / static_cast<float>(nt));
      }
    }
  } else {
    Gf = pf * (arma::ones<arma::frowvec>(nt) / static_cast<float>(nt));
  }

  arma::fmat G_prev = Gf;
  arma::fmat scratch(ns, nt);
  arma::fmat scratch2;
  if (!symmetric) {
    scratch2.set_size(ns, nt);
  }
  arma::fmat grad_quad(ns, nt);
  arma::fmat grad1;
  arma::fmat grad2;
  if (!symmetric) {
    grad1.set_size(ns, nt);
    grad2.set_size(ns, nt);
  }
  arma::fmat K(ns, nt);
  arma::fvec rs(ns);

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (int k = 0; k < max_iter; ++k) {
    const bool do_check = ((k + 1) % check_every) == 0 || (k == 0);
    if (do_check) {
      G_prev = Gf;
    }

    const arma::frowvec qG = arma::sum(Gf, 0);
    if (symmetric) {
      const arma::frowvec marg = qG * fC2t;
      grad_quad = constC;
      grad_quad.each_row() += marg;
      sgemm_nn(hC1, Gf, scratch);
      sgemm_nt_accum(scratch, hC2, grad_quad, -1.0f, 1.0f);
    } else {
      const arma::frowvec marg1 = qG * fC2t;
      grad1 = constC;
      grad1.each_row() += marg1;
      sgemm_nn(hC1, Gf, scratch);
      sgemm_nt_accum(scratch, hC2, grad1, -1.0f, 1.0f);
      grad1 *= 2.0f;

      const arma::frowvec marg2 = qG * fC2;
      grad2 = constCt;
      grad2.each_row() += marg2;
      sgemm_nn(hC1t, Gf, scratch2);
      sgemm_nt_accum(scratch2, hC2t, grad2, -1.0f, 1.0f);
      grad2 *= 2.0f;

      grad_quad = 0.5f * (grad1 + grad2);
    }

    if (alpha_full) {
      const float inv_eps_grad = inv_eps_f * (symmetric ? 2.0f : 1.0f);
      build_kernel_from_cost_mul_f(grad_quad, Gf, inv_eps_grad, K);
    } else if (alpha_zero) {
      build_kernel_from_cost_mul_f(Mf, Gf, inv_eps_f, K);
    } else {
      build_kernel_from_fused_cost_mul_f(grad_quad, Mf, Gf, inv_eps_f, quad_scale_f, alpha_lin_f, K);
    }
    rs = arma::sum(K, 1);
    for (arma::uword i = 0; i < ns; ++i) {
      const float rsi = rs[i];
      if (std::isfinite(rsi) && rsi > 0.0f) {
        Gf.row(i) = (pf[i] / rsi) * K.row(i);
      } else {
        Gf.row(i).fill(pf[i] / static_cast<float>(nt));
      }
    }

    if (do_check) {
      err = static_cast<double>(arma::norm(Gf - G_prev, "fro"));
      if (err <= tol) {
        it = k + 1;
        break;
      }
    }
    it = k + 1;
  }

  arma::mat G = arma::conv_to<arma::mat>::from(Gf);
  normalize_rows_to_p(G, p, nt);
  const arma::mat C1_sqr_d = C1 % C1;
  const arma::mat C2_sqr_d = C2 % C2;
  const arma::mat fC2t_d = C2_sqr_d.t();
  const arma::vec left_d = C1_sqr_d * p;
  arma::mat constC_d(ns, nt);
  for (arma::uword j = 0; j < nt; ++j) {
    constC_d.col(j) = left_d;
  }
  const arma::mat hC1_d = C1;
  const arma::mat hC2_d = 2.0 * C2;

  arma::mat tens = constC_d;
  const arma::rowvec qG = arma::sum(G, 0);
  tens.each_row() += qG * fC2t_d;
  arma::mat scratch_d(ns, nt);
  dgemm_nn(hC1_d, G, scratch_d);
  dgemm_nt_accum(scratch_d, hC2_d, tens, -1.0, 1.0);
  const double quad_raw = arma::accu(tens % G);
  const double lin_loss = (alpha >= 1.0 || M.n_elem == 0) ? 0.0 : (1.0 - alpha) * arma::accu(M % G);
  const double quad_loss = alpha * quad_raw;

  SrfgwEntropicCoreResult out;
  out.plan = G;
  out.q = qG.t();
  out.lin_loss = lin_loss;
  out.quad_loss = quad_loss;
  out.srfgw_dist = lin_loss + quad_loss;
  out.srgw_dist = quad_raw;
  out.iterations = it;
  out.error = err;
  return out;
}

inline FgwEntropicCoreResult fgw_entropic_square_double_core(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    bool use_log_sinkhorn,
    int check_every,
    int approx_rank,
    const arma::mat& init_plan,
    const arma::mat* c1_A1_cache_d,
    const arma::mat* c1_B1t_cache_d,
    const LowRankC2CacheD* c2_cache_d,
    const SquareC2CacheD* c2_square_cache_d) {
  if (check_every <= 0) {
    check_every = 1;
  }

  arma::mat constC, hC1, hC2_local;
  const arma::mat* hC2 = nullptr;
  const bool used_square_cache = symmetric &&
    (c2_square_cache_d != nullptr) &&
    square_cache_compatible(*c2_square_cache_d, C2, q);
  if (used_square_cache) {
    init_matrices_square_from_cache_d(C1, p, *c2_square_cache_d, constC, hC1);
    hC2 = &(c2_square_cache_d->hC2);
  } else {
    init_matrices_square(C1, C2, p, q, constC, hC1, hC2_local);
    hC2 = &hC2_local;
  }

  arma::mat constCt, hC1t, hC2_asym;
  if (!symmetric) {
    init_matrices_square(C1.t(), C2.t(), p, q, constCt, hC1t, hC2_asym);
  }
  const arma::mat M_lin = (1.0 - alpha) * M;
  bool use_lowrank = false;
  bool used_c1_cache = false;
  bool used_c2_cache = false;
  arma::mat lr_A1_local, lr_B1t_local, lr_A2, lr_A2t, lr_B2t;
  const arma::mat* lr_A1 = nullptr;
  const arma::mat* lr_B1t = nullptr;
  if (symmetric && approx_rank > 0) {
    const bool can_use_c2_cache = (c2_cache_d != nullptr) && lowrank_cache_compatible(*c2_cache_d, C2);
    const bool can_use_c1_cache = (c1_A1_cache_d != nullptr) &&
      (c1_B1t_cache_d != nullptr) &&
      lowrank_c1_cache_compatible(*c1_A1_cache_d, *c1_B1t_cache_d, C1);
    if (can_use_c1_cache) {
      lr_A1 = c1_A1_cache_d;
      lr_B1t = c1_B1t_cache_d;
      use_lowrank = true;
      used_c1_cache = true;
    } else {
      use_lowrank = svd_lowrank_factors(C1, approx_rank, lr_A1_local, lr_B1t_local);
      if (use_lowrank) {
        lr_A1 = &lr_A1_local;
        lr_B1t = &lr_B1t_local;
      }
    }
    if (use_lowrank) {
      if (can_use_c2_cache) {
        lr_A2t = c2_cache_d->A2t_scaled;
        lr_B2t = c2_cache_d->B2t;
        used_c2_cache = true;
      } else {
        use_lowrank = svd_lowrank_factors(C2, approx_rank, lr_A2, lr_B2t);
        if (use_lowrank) {
          lr_A2t = lr_A2.t();
        }
      }
    }
  }
  const double symmetric_cross_scale = use_lowrank ? (-4.0 * alpha) : (-2.0 * alpha);
  const double asym_cross_scale = -alpha;
  arma::mat base_sym;
  arma::mat base_asym;
  if (symmetric) {
    base_sym = (2.0 * alpha) * constC;
    base_sym += M_lin;
  } else {
    base_asym = alpha * (constC + constCt);
    base_asym += M_lin;
  }
  arma::mat T;
  if (init_plan.n_rows == p.n_elem && init_plan.n_cols == q.n_elem) {
    T = init_plan;
    T.for_each([](double& v) {
      if (!std::isfinite(v) || v < 0.0) {
        v = 0.0;
      }
    });
    const double m = arma::accu(T);
    if (m > 0.0) {
      T /= m;
    } else {
      T = p * q.t();
    }
  } else {
    T = p * q.t();
  }
  arma::mat T_prev_check = T;
  arma::vec u_ws = arma::ones<arma::vec>(p.n_elem);
  arma::vec v_ws = arma::ones<arma::vec>(q.n_elem);
  arma::vec f_ws = arma::zeros<arma::vec>(p.n_elem);
  arma::vec g_ws = arma::zeros<arma::vec>(q.n_elem);
  arma::mat tens(C1.n_rows, C2.n_rows);
  arma::mat scratch(C1.n_rows, C2.n_rows);
  arma::mat scratch_asym;
  if (!symmetric) {
    scratch_asym.set_size(C1.n_rows, C2.n_rows);
  }
  arma::mat lr_tmp_r1_nt;
  arma::mat lr_tmp_r1_r2;
  arma::mat lr_tmp_ns_r2;
  if (use_lowrank) {
    lr_tmp_r1_nt.set_size(lr_B1t->n_rows, C2.n_rows);
    lr_tmp_r1_r2.set_size(lr_B1t->n_rows, lr_B2t.n_rows);
    lr_tmp_ns_r2.set_size(C1.n_rows, lr_B2t.n_rows);
  }
  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if ((it + 1) % check_every == 0) {
      T_prev_check = T;
    }

    if (symmetric) {
      if (use_lowrank) {
        tensor_product_lowrank_scaled(
          base_sym, *lr_A1, *lr_B1t, lr_A2t, lr_B2t, T, symmetric_cross_scale,
          tens, lr_tmp_r1_nt, lr_tmp_r1_r2, lr_tmp_ns_r2
        );
      } else {
        tensor_product_blas_scaled(
          base_sym, hC1, *hC2, T, symmetric_cross_scale, tens, scratch
        );
      }
    } else {
      tensor_product_asym_blas_scaled(
        base_asym, hC1, *hC2, hC1t, hC2_asym, T, asym_cross_scale,
        tens, scratch, scratch_asym
      );
    }

    if (use_ppa) {
      tens -= epsilon * arma::log(T + kTiny);
    }

    SinkhornBalancedResult sk;
    if (use_log_sinkhorn) {
      sk = sinkhorn_balanced_log(
        p, q, tens, epsilon, sinkhorn_max_iter, sinkhorn_tol, f_ws, g_ws
      );
    } else {
      sk = sinkhorn_balanced(
        p, q, tens, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
      );
    }
    T = std::move(sk.plan);
    if (use_log_sinkhorn) {
      f_ws = std::move(sk.f);
      g_ws = std::move(sk.g);
    } else {
      u_ws = std::move(sk.u);
      v_ws = std::move(sk.v);
    }

    if ((it + 1) % check_every == 0) {
      err = arma::norm(T - T_prev_check, "fro");
      if (err < tol) {
        ++it;
        break;
      }
    }
  }

  FgwEntropicCoreResult out;
  out.fgw_dist = (1.0 - alpha) * arma::accu(M % T) + alpha * gw_loss(constC, hC1, *hC2, T);
  out.plan = std::move(T);
  out.iterations = it;
  out.error = err;
  out.used_lowrank = use_lowrank;
  out.used_c1_cache = used_c1_cache;
  out.used_c2_cache = used_c2_cache;
  out.used_square_cache = used_square_cache;
  return out;
}

inline FgwEntropicCoreResult fgw_entropic_square_mixed_core(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    int check_every,
    int approx_rank,
    const arma::mat& init_plan,
    const arma::fmat* M_cache_f,
    const arma::fmat* C1_cache_f,
    const arma::fvec* p_cache_f,
    const arma::fmat* c1_A1_cache_f,
    const arma::fmat* c1_B1t_cache_f,
    const LowRankC2CacheF* c2_cache_f,
    const SquareC2CacheF* c2_square_cache_f) {
  if (check_every <= 0) {
    check_every = 1;
  }

  const float alpha_f = static_cast<float>(alpha);
  const float epsilon_f = static_cast<float>(epsilon);
  const float sinkhorn_tol_f = static_cast<float>(std::max(sinkhorn_tol, 1e-6));
  const float tol_f = static_cast<float>(std::max(tol, 1e-6));

  arma::fmat Mf_local;
  arma::fmat C1f_local;
  arma::fvec pf_local;
  const arma::fmat* Mf_ptr = nullptr;
  const arma::fmat* C1f_ptr = nullptr;
  const arma::fvec* pf_ptr = nullptr;
  if (M_cache_f != nullptr &&
      M_cache_f->n_rows == M.n_rows &&
      M_cache_f->n_cols == M.n_cols) {
    Mf_ptr = M_cache_f;
  } else {
    Mf_local = arma::conv_to<arma::fmat>::from(M);
    Mf_ptr = &Mf_local;
  }
  if (C1_cache_f != nullptr &&
      C1_cache_f->n_rows == C1.n_rows &&
      C1_cache_f->n_cols == C1.n_cols) {
    C1f_ptr = C1_cache_f;
  } else {
    C1f_local = arma::conv_to<arma::fmat>::from(C1);
    C1f_ptr = &C1f_local;
  }
  if (p_cache_f != nullptr && p_cache_f->n_elem == p.n_elem) {
    pf_ptr = p_cache_f;
  } else {
    pf_local = arma::conv_to<arma::fvec>::from(p);
    pf_ptr = &pf_local;
  }
  const arma::fmat& Mf = *Mf_ptr;
  const arma::fmat& C1f = *C1f_ptr;
  const arma::fvec& pf = *pf_ptr;

  const bool can_use_square_cache = symmetric &&
    (c2_square_cache_f != nullptr) &&
    square_cache_compatible_f(*c2_square_cache_f, C2, q);
  arma::fmat C2f_local;
  arma::fvec qf_local;
  const arma::fmat* C2f = nullptr;
  const arma::fvec* qf = nullptr;
  if (can_use_square_cache) {
    C2f = &(c2_square_cache_f->C2);
    qf = &(c2_square_cache_f->q);
  } else {
    C2f_local = arma::conv_to<arma::fmat>::from(C2);
    qf_local = arma::conv_to<arma::fvec>::from(q);
    C2f = &C2f_local;
    qf = &qf_local;
  }

  arma::fmat constC, hC1, hC2_local;
  const arma::fmat* hC2 = nullptr;
  if (can_use_square_cache) {
    init_matrices_square_from_cache_f(C1f, pf, *c2_square_cache_f, constC, hC1);
    hC2 = &(c2_square_cache_f->hC2);
  } else {
    init_matrices_square_f(C1f, *C2f, pf, *qf, constC, hC1, hC2_local);
    hC2 = &hC2_local;
  }
  arma::fmat constCt, hC1t, hC2_asym;
  if (!symmetric) {
    init_matrices_square_f(C1f.t(), C2f->t(), pf, *qf, constCt, hC1t, hC2_asym);
  }
  const arma::fmat M_lin = (1.0f - alpha_f) * Mf;
  bool use_lowrank = false;
  bool used_c1_cache = false;
  bool used_c2_cache = false;
  arma::fmat lr_A1_local, lr_B1t_local, lr_A2, lr_A2t, lr_B2t;
  const arma::fmat* lr_A1 = nullptr;
  const arma::fmat* lr_B1t = nullptr;
  if (symmetric && approx_rank > 0) {
    const bool can_use_c2_cache = (c2_cache_f != nullptr) && lowrank_cache_compatible_f(*c2_cache_f, *C2f);
    const bool can_use_c1_cache = (c1_A1_cache_f != nullptr) &&
      (c1_B1t_cache_f != nullptr) &&
      lowrank_c1_cache_compatible_f(*c1_A1_cache_f, *c1_B1t_cache_f, C1f);
    if (can_use_c1_cache) {
      lr_A1 = c1_A1_cache_f;
      lr_B1t = c1_B1t_cache_f;
      use_lowrank = true;
      used_c1_cache = true;
    } else {
      use_lowrank = svd_lowrank_factors_f(C1f, approx_rank, lr_A1_local, lr_B1t_local);
      if (use_lowrank) {
        lr_A1 = &lr_A1_local;
        lr_B1t = &lr_B1t_local;
      }
    }
    if (use_lowrank) {
      if (can_use_c2_cache) {
        lr_A2t = c2_cache_f->A2t_scaled;
        lr_B2t = c2_cache_f->B2t;
        used_c2_cache = true;
      } else {
        use_lowrank = svd_lowrank_factors_f(*C2f, approx_rank, lr_A2, lr_B2t);
        if (use_lowrank) {
          lr_A2t = lr_A2.t();
        }
      }
    }
  }
  const float symmetric_cross_scale = use_lowrank ? (-4.0f * alpha_f) : (-2.0f * alpha_f);
  const float asym_cross_scale = -alpha_f;
  arma::fmat base_sym;
  arma::fmat base_asym;
  if (symmetric) {
    base_sym = (2.0f * alpha_f) * constC;
    base_sym += M_lin;
  } else {
    base_asym = alpha_f * (constC + constCt);
    base_asym += M_lin;
  }

  arma::fmat T;
  if (init_plan.n_rows == p.n_elem && init_plan.n_cols == q.n_elem) {
    T = arma::conv_to<arma::fmat>::from(init_plan);
    T.for_each([](float& v) {
      if (!std::isfinite(v) || v < 0.0f) {
        v = 0.0f;
      }
    });
    const float m = arma::accu(T);
    if (m > 0.0f) {
      T /= m;
    } else {
      T = pf * qf->t();
    }
  } else {
    T = pf * qf->t();
  }
  arma::fmat T_prev_check = T;
  arma::fvec u_ws = arma::ones<arma::fvec>(pf.n_elem);
  arma::fvec v_ws = arma::ones<arma::fvec>(qf->n_elem);
  arma::fmat tens(C1f.n_rows, C2f->n_rows);
  arma::fmat scratch(C1f.n_rows, C2f->n_rows);
  arma::fmat scratch_asym;
  if (!symmetric) {
    scratch_asym.set_size(C1f.n_rows, C2f->n_rows);
  }
  arma::fmat lr_tmp_r1_nt;
  arma::fmat lr_tmp_r1_r2;
  arma::fmat lr_tmp_ns_r2;
  if (use_lowrank) {
    lr_tmp_r1_nt.set_size(lr_B1t->n_rows, C2f->n_rows);
    lr_tmp_r1_r2.set_size(lr_B1t->n_rows, lr_B2t.n_rows);
    lr_tmp_ns_r2.set_size(C1f.n_rows, lr_B2t.n_rows);
  }

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if ((it + 1) % check_every == 0) {
      T_prev_check = T;
    }

    if (symmetric) {
      if (use_lowrank) {
        tensor_product_lowrank_scaled_f(
          base_sym, *lr_A1, *lr_B1t, lr_A2t, lr_B2t, T, symmetric_cross_scale,
          tens, lr_tmp_r1_nt, lr_tmp_r1_r2, lr_tmp_ns_r2
        );
      } else {
        tensor_product_blas_scaled_f(
          base_sym, hC1, *hC2, T, symmetric_cross_scale, tens, scratch
        );
      }
    } else {
      tensor_product_asym_blas_scaled_f(
        base_asym, hC1, *hC2, hC1t, hC2_asym, T, asym_cross_scale,
        tens, scratch, scratch_asym
      );
    }

    if (use_ppa) {
      tens -= epsilon_f * arma::log(T + kTinyF);
    }

    SinkhornBalancedResultF sk = sinkhorn_balanced_f(
      pf, *qf, tens, epsilon_f, sinkhorn_max_iter, sinkhorn_tol_f, u_ws, v_ws
    );
    T = std::move(sk.plan);
    u_ws = std::move(sk.u);
    v_ws = std::move(sk.v);

    if ((it + 1) % check_every == 0) {
      err = static_cast<double>(arma::norm(T - T_prev_check, "fro"));
      if (err < tol_f) {
        ++it;
        break;
      }
    }
  }

  arma::mat Td = arma::conv_to<arma::mat>::from(T);
  arma::mat constCd, hC1d, hC2d;
  init_matrices_square(C1, C2, p, q, constCd, hC1d, hC2d);

  FgwEntropicCoreResult out;
  out.fgw_dist = (1.0 - alpha) * arma::accu(M % Td) + alpha * gw_loss(constCd, hC1d, hC2d, Td);
  out.plan = std::move(Td);
  out.iterations = it;
  out.error = err;
  out.used_lowrank = use_lowrank;
  out.used_c1_cache = used_c1_cache;
  out.used_c2_cache = used_c2_cache;
  out.used_square_cache = can_use_square_cache;
  return out;
}

inline FgwEntropicCoreResult fgw_entropic_square_core(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    bool use_log_sinkhorn,
    bool use_mixed_precision,
    int check_every,
    int approx_rank,
    const arma::mat& init_plan,
    const arma::mat* c1_A1_cache_d,
    const arma::mat* c1_B1t_cache_d,
    const arma::fmat* c1_A1_cache_f,
    const arma::fmat* c1_B1t_cache_f,
    const arma::fmat* M_cache_f,
    const arma::fmat* C1_cache_f,
    const arma::fvec* p_cache_f,
    const LowRankC2CacheD* c2_cache_d,
    const LowRankC2CacheF* c2_cache_f,
    const SquareC2CacheD* c2_square_cache_d,
    const SquareC2CacheF* c2_square_cache_f) {
  const bool use_double_mixed_accel =
    (!use_mixed_precision) &&
    (!use_log_sinkhorn) &&
    (!use_ppa) &&
    (C1.n_rows >= 32) &&
    (C2.n_rows >= 32);
  if ((use_mixed_precision || use_double_mixed_accel) && !use_log_sinkhorn) {
    return fgw_entropic_square_mixed_core(
      M, C1, C2, p, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, check_every, approx_rank, init_plan,
      M_cache_f, C1_cache_f, p_cache_f,
      c1_A1_cache_f, c1_B1t_cache_f, c2_cache_f, c2_square_cache_f
    );
  }
  return fgw_entropic_square_double_core(
    M, C1, C2, p, q,
    alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
    symmetric, use_ppa, use_log_sinkhorn, check_every, approx_rank, init_plan,
    c1_A1_cache_d, c1_B1t_cache_d, c2_cache_d, c2_square_cache_d
  );
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_fgw_entropic_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    bool use_log_sinkhorn,
    bool use_mixed_precision,
    int check_every,
    int approx_rank,
    const arma::mat& init_plan) {
  if (check_every <= 0) {
    check_every = 1;
  }
  if (approx_rank < 0) {
    approx_rank = 0;
  }
  const bool use_double_mixed_accel =
    (!use_mixed_precision) &&
    (!use_log_sinkhorn) &&
    (!use_ppa) &&
    (C1.n_rows >= 32) &&
    (C2.n_rows >= 32);
  if ((use_mixed_precision || use_double_mixed_accel) && !use_log_sinkhorn) {
    return fgw_entropic_square_mixed_impl(
      M, C1, C2, p, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, check_every, approx_rank, init_plan
    );
  }

  arma::mat constC, hC1, hC2;
  init_matrices_square(C1, C2, p, q, constC, hC1, hC2);

  arma::mat constCt, hC1t, hC2_asym;
  if (!symmetric) {
    init_matrices_square(C1.t(), C2.t(), p, q, constCt, hC1t, hC2_asym);
  }
  const arma::mat M_lin = (1.0 - alpha) * M;
  bool use_lowrank = false;
  arma::mat lr_A1, lr_B1t, lr_A2, lr_A2t, lr_B2t;
  if (symmetric && approx_rank > 0) {
    use_lowrank = svd_lowrank_factors(C1, approx_rank, lr_A1, lr_B1t) &&
      svd_lowrank_factors(C2, approx_rank, lr_A2, lr_B2t);
    if (use_lowrank) {
      lr_A2t = lr_A2.t();
    }
  }
  const double symmetric_cross_scale = use_lowrank ? (-4.0 * alpha) : (-2.0 * alpha);
  const double asym_cross_scale = -alpha;
  arma::mat base_sym;
  arma::mat base_asym;
  if (symmetric) {
    base_sym = (2.0 * alpha) * constC;
    base_sym += M_lin;
  } else {
    base_asym = alpha * (constC + constCt);
    base_asym += M_lin;
  }

  arma::mat T;
  if (init_plan.n_rows == p.n_elem && init_plan.n_cols == q.n_elem) {
    T = init_plan;
    T.for_each([](double& v) {
      if (!std::isfinite(v) || v < 0.0) {
        v = 0.0;
      }
    });
    const double m = arma::accu(T);
    if (m > 0.0) {
      T /= m;
    } else {
      T = p * q.t();
    }
  } else {
    T = p * q.t();
  }
  arma::mat T_prev_check = T;
  arma::vec u_ws = arma::ones<arma::vec>(p.n_elem);
  arma::vec v_ws = arma::ones<arma::vec>(q.n_elem);
  arma::vec f_ws = arma::zeros<arma::vec>(p.n_elem);
  arma::vec g_ws = arma::zeros<arma::vec>(q.n_elem);
  Rcpp::NumericVector err_trace;
  arma::mat tens(C1.n_rows, C2.n_rows);
  arma::mat scratch(C1.n_rows, C2.n_rows);
  arma::mat scratch_asym;
  if (!symmetric) {
    scratch_asym.set_size(C1.n_rows, C2.n_rows);
  }
  arma::mat lr_tmp_r1_nt;
  arma::mat lr_tmp_r1_r2;
  arma::mat lr_tmp_ns_r2;
  if (use_lowrank) {
    lr_tmp_r1_nt.set_size(lr_B1t.n_rows, C2.n_rows);
    lr_tmp_r1_r2.set_size(lr_B1t.n_rows, lr_B2t.n_rows);
    lr_tmp_ns_r2.set_size(C1.n_rows, lr_B2t.n_rows);
  }

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    if ((it + 1) % check_every == 0) {
      T_prev_check = T;
    }

    if (symmetric) {
      if (use_lowrank) {
        tensor_product_lowrank_scaled(
          base_sym, lr_A1, lr_B1t, lr_A2t, lr_B2t, T, symmetric_cross_scale,
          tens, lr_tmp_r1_nt, lr_tmp_r1_r2, lr_tmp_ns_r2
        );
      } else {
        tensor_product_blas_scaled(
          base_sym, hC1, hC2, T, symmetric_cross_scale, tens, scratch
        );
      }
    } else {
      tensor_product_asym_blas_scaled(
        base_asym, hC1, hC2, hC1t, hC2_asym, T, asym_cross_scale,
        tens, scratch, scratch_asym
      );
    }

    if (use_ppa) {
      tens -= epsilon * arma::log(T + kTiny);
    }

    SinkhornBalancedResult sk;
    if (use_log_sinkhorn) {
      sk = sinkhorn_balanced_log(
        p, q, tens, epsilon, sinkhorn_max_iter, sinkhorn_tol, f_ws, g_ws
      );
    } else {
      sk = sinkhorn_balanced(
        p, q, tens, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
      );
    }
    T = std::move(sk.plan);
    if (use_log_sinkhorn) {
      f_ws = std::move(sk.f);
      g_ws = std::move(sk.g);
    } else {
      u_ws = std::move(sk.u);
      v_ws = std::move(sk.v);
    }

    if ((it + 1) % check_every == 0) {
      err = arma::norm(T - T_prev_check, "fro");
      err_trace.push_back(err);
      if (err < tol) {
        ++it;
        break;
      }
    }
  }

  const double fgw_dist = (1.0 - alpha) * arma::accu(M % T) + alpha * gw_loss(constC, hC1, hC2, T);
  return Rcpp::List::create(
    Rcpp::Named("plan") = T,
    Rcpp::Named("fgw_dist") = fgw_dist,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = err,
    Rcpp::Named("err_trace") = err_trace
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_entropic_semirelaxed_fgw_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    double epsilon,
    double alpha,
    bool symmetric,
    bool use_mixed_precision,
    int max_iter,
    double tol,
    int check_every,
    const arma::mat& init_plan) {
  if (C1.n_rows != C1.n_cols) {
    Rcpp::stop("`C1` must be square.");
  }
  if (C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C2` must be square.");
  }
  const bool skip_M = (alpha >= 1.0) && (M.n_elem == 0);
  if (!skip_M && (M.n_rows != C1.n_rows || M.n_cols != C2.n_rows)) {
    Rcpp::stop("`M` must have shape nrow(C1) x nrow(C2).");
  }
  if (p.n_elem != C1.n_rows) {
    Rcpp::stop("`p` must have length nrow(C1).");
  }
  if (!std::isfinite(epsilon) || epsilon <= 0.0) {
    Rcpp::stop("`epsilon` must be positive.");
  }
  if (!std::isfinite(alpha) || alpha < 0.0 || alpha > 1.0) {
    Rcpp::stop("`alpha` must be in [0, 1].");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }
  if (check_every <= 0) {
    check_every = 1;
  }

  arma::vec p_norm = p;
  bool p_invalid = false;
  for (arma::uword i = 0; i < p_norm.n_elem; ++i) {
    const double v = p_norm[i];
    if (!std::isfinite(v) || v < 0.0) {
      p_invalid = true;
      break;
    }
  }
  if (p_invalid) {
    Rcpp::stop("`p` must be finite and nonnegative.");
  }
  const double p_sum = arma::accu(p_norm);
  if (!std::isfinite(p_sum) || p_sum <= 0.0) {
    Rcpp::stop("`p` must have positive total mass.");
  }
  p_norm /= p_sum;

  const bool use_double_mixed_accel =
    (!use_mixed_precision) &&
    (C1.n_rows >= 32) &&
    (C2.n_rows >= 32);
  SrfgwEntropicCoreResult out = (use_mixed_precision || use_double_mixed_accel) ?
    entropic_semirelaxed_fgw_square_core_mixed(
      M, C1, C2, p_norm, epsilon, alpha, symmetric, max_iter, tol, check_every, init_plan
    ) :
    entropic_semirelaxed_fgw_square_core(
      M, C1, C2, p_norm, epsilon, alpha, symmetric, max_iter, tol, check_every, init_plan
    );
  return Rcpp::List::create(
    Rcpp::Named("plan") = out.plan,
    Rcpp::Named("q") = out.q,
    Rcpp::Named("lin_loss") = out.lin_loss,
    Rcpp::Named("quad_loss") = out.quad_loss,
    Rcpp::Named("srfgw_dist") = out.srfgw_dist,
    Rcpp::Named("srgw_dist") = out.srgw_dist,
    Rcpp::Named("iterations") = out.iterations,
    Rcpp::Named("error") = out.error
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_fgw_entropic_square_batch(
    const Rcpp::List& M_list,
    const Rcpp::List& C1_list,
    const Rcpp::List& p_list,
    const arma::mat& C2,
    const arma::vec& q,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    bool use_log_sinkhorn,
    bool use_mixed_precision,
    int check_every,
    const Rcpp::List& init_plan_list,
    const Rcpp::List& c1_A_scaled_list,
    const Rcpp::List& c1_Bt_list,
    int approx_rank,
    int n_threads) {
  const int n_jobs = static_cast<int>(M_list.size());
  if (n_jobs <= 0) {
    Rcpp::stop("`M_list` must be non-empty.");
  }
  if (C1_list.size() != M_list.size() || p_list.size() != M_list.size()) {
    Rcpp::stop("`M_list`, `C1_list`, and `p_list` must have the same length.");
  }
  if (init_plan_list.size() > 0 && init_plan_list.size() != M_list.size()) {
    Rcpp::stop("`init_plan_list` must be empty or the same length as `M_list`.");
  }
  if ((c1_A_scaled_list.size() > 0 && c1_A_scaled_list.size() != M_list.size()) ||
      (c1_Bt_list.size() > 0 && c1_Bt_list.size() != M_list.size())) {
    Rcpp::stop("`c1_A_scaled_list` and `c1_Bt_list` must be empty or the same length as `M_list`.");
  }
  if ((c1_A_scaled_list.size() > 0) != (c1_Bt_list.size() > 0)) {
    Rcpp::stop("Provide both `c1_A_scaled_list` and `c1_Bt_list` together, or neither.");
  }
  if (q.n_elem != C2.n_rows || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C2` must be square and compatible with `q`.");
  }
  if (check_every <= 0) {
    check_every = 1;
  }
  if (approx_rank < 0) {
    approx_rank = 0;
  }
  if (n_threads <= 0) {
    n_threads = 1;
  }

  std::vector<Rcpp::NumericMatrix> M_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> C1_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericVector> p_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> init_refs(static_cast<std::size_t>(n_jobs));
  std::vector<unsigned char> has_init(static_cast<std::size_t>(n_jobs), 0);
  std::vector<Rcpp::NumericMatrix> c1_A_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> c1_Bt_refs(static_cast<std::size_t>(n_jobs));
  std::vector<unsigned char> has_c1_cache(static_cast<std::size_t>(n_jobs), 0);
  for (int i = 0; i < n_jobs; ++i) {
    M_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(M_list[i]);
    C1_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(C1_list[i]);
    p_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericVector>(p_list[i]);
    if (init_plan_list.size() == M_list.size()) {
      SEXP init_i = init_plan_list[i];
      if (init_i != R_NilValue) {
        init_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(init_i);
        has_init[static_cast<std::size_t>(i)] = 1;
      }
    }
    if (c1_A_scaled_list.size() == M_list.size()) {
      SEXP c1a_i = c1_A_scaled_list[i];
      SEXP c1b_i = c1_Bt_list[i];
      if (c1a_i != R_NilValue && c1b_i != R_NilValue) {
        c1_A_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(c1a_i);
        c1_Bt_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(c1b_i);
        has_c1_cache[static_cast<std::size_t>(i)] = 1;
      }
    }
    const Rcpp::NumericMatrix& Mi = M_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericMatrix& C1i = C1_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericVector& pi = p_refs[static_cast<std::size_t>(i)];
    if (C1i.nrow() != C1i.ncol()) {
      Rcpp::stop("Every `C1_list[[i]]` must be square.");
    }
    if (static_cast<int>(pi.size()) != C1i.nrow()) {
      Rcpp::stop("Each `p_list[[i]]` must match `nrow(C1_list[[i]])`.");
    }
    if (Mi.nrow() != C1i.nrow() || Mi.ncol() != static_cast<int>(C2.n_rows)) {
      Rcpp::stop("Each `M_list[[i]]` must have shape nrow(C1_list[[i]]) x nrow(C2).");
    }
    if (has_c1_cache[static_cast<std::size_t>(i)] == 1) {
      const Rcpp::NumericMatrix& c1A = c1_A_refs[static_cast<std::size_t>(i)];
      const Rcpp::NumericMatrix& c1B = c1_Bt_refs[static_cast<std::size_t>(i)];
      if (c1A.nrow() != C1i.nrow() || c1B.ncol() != C1i.ncol() || c1A.ncol() != c1B.nrow() || c1A.ncol() <= 0) {
        Rcpp::stop("Each cached C1 low-rank pair must satisfy A(n x r), Bt(r x n).");
      }
    }
  }

  std::vector<arma::mat> plans(static_cast<std::size_t>(n_jobs));
  std::vector<double> dists(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<double> errs(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<int> iters(static_cast<std::size_t>(n_jobs), 0);
  std::vector<double> kernel_ms(static_cast<std::size_t>(n_jobs), 0.0);
  std::vector<int> lowrank_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> c1_cache_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> c2_cache_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> square_cache_used(static_cast<std::size_t>(n_jobs), 0);
  SquareC2CacheD c2_square_cache_d;
  SquareC2CacheF c2_square_cache_f;
  LowRankC2CacheD c2_cache_d;
  LowRankC2CacheF c2_cache_f;
  const bool use_mixed_kernel = use_mixed_precision && !use_log_sinkhorn;
  if (symmetric) {
    if (use_mixed_kernel) {
      const arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
      const arma::fvec qf = arma::conv_to<arma::fvec>::from(q);
      c2_square_cache_f = build_square_c2_cache_f(C2f, qf);
    } else {
      c2_square_cache_d = build_square_c2_cache_d(C2, q);
    }
  }
  const bool want_lowrank = symmetric && (approx_rank > 0);
  if (want_lowrank) {
    if (use_mixed_kernel) {
      if (c2_square_cache_f.valid) {
        c2_cache_f = build_lowrank_c2_cache_f(c2_square_cache_f.C2, approx_rank);
      } else {
        const arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
        c2_cache_f = build_lowrank_c2_cache_f(C2f, approx_rank);
      }
    } else {
      c2_cache_d = build_lowrank_c2_cache_d(C2, approx_rank);
    }
  }
  const LowRankC2CacheD* c2_cache_d_ptr = c2_cache_d.valid ? &c2_cache_d : nullptr;
  const LowRankC2CacheF* c2_cache_f_ptr = c2_cache_f.valid ? &c2_cache_f : nullptr;
  const SquareC2CacheD* c2_square_cache_d_ptr = c2_square_cache_d.valid ? &c2_square_cache_d : nullptr;
  const SquareC2CacheF* c2_square_cache_f_ptr = c2_square_cache_f.valid ? &c2_square_cache_f : nullptr;
  std::vector<arma::fmat> c1_A_cache_f;
  std::vector<arma::fmat> c1_Bt_cache_f;
  std::vector<arma::fmat> M_kernel_cache_f;
  std::vector<arma::fmat> C1_kernel_cache_f;
  std::vector<arma::fvec> p_kernel_cache_f;
  if (use_mixed_kernel) {
    M_kernel_cache_f.resize(static_cast<std::size_t>(n_jobs));
    C1_kernel_cache_f.resize(static_cast<std::size_t>(n_jobs));
    p_kernel_cache_f.resize(static_cast<std::size_t>(n_jobs));
    for (int i = 0; i < n_jobs; ++i) {
      const std::size_t idx = static_cast<std::size_t>(i);
      const Rcpp::NumericMatrix& Mref = M_refs[idx];
      const Rcpp::NumericMatrix& C1ref = C1_refs[idx];
      const Rcpp::NumericVector& pref = p_refs[idx];
      const arma::mat Mi(const_cast<double*>(Mref.begin()), static_cast<arma::uword>(Mref.nrow()), static_cast<arma::uword>(Mref.ncol()), false, true);
      const arma::mat C1i(const_cast<double*>(C1ref.begin()), static_cast<arma::uword>(C1ref.nrow()), static_cast<arma::uword>(C1ref.ncol()), false, true);
      const arma::vec pi(const_cast<double*>(pref.begin()), static_cast<arma::uword>(pref.size()), false, true);
      M_kernel_cache_f[idx] = arma::conv_to<arma::fmat>::from(Mi);
      C1_kernel_cache_f[idx] = arma::conv_to<arma::fmat>::from(C1i);
      p_kernel_cache_f[idx] = arma::conv_to<arma::fvec>::from(pi);
    }
  }
  if (use_mixed_kernel && c1_A_scaled_list.size() == M_list.size()) {
    c1_A_cache_f.resize(static_cast<std::size_t>(n_jobs));
    c1_Bt_cache_f.resize(static_cast<std::size_t>(n_jobs));
    for (int i = 0; i < n_jobs; ++i) {
      const std::size_t idx = static_cast<std::size_t>(i);
      if (has_c1_cache[idx] == 0) {
        continue;
      }
      const Rcpp::NumericMatrix& c1A = c1_A_refs[idx];
      const Rcpp::NumericMatrix& c1B = c1_Bt_refs[idx];
      const arma::mat c1A_d(const_cast<double*>(c1A.begin()), static_cast<arma::uword>(c1A.nrow()), static_cast<arma::uword>(c1A.ncol()), false, true);
      const arma::mat c1B_d(const_cast<double*>(c1B.begin()), static_cast<arma::uword>(c1B.nrow()), static_cast<arma::uword>(c1B.ncol()), false, true);
      c1_A_cache_f[idx] = arma::conv_to<arma::fmat>::from(c1A_d);
      c1_Bt_cache_f[idx] = arma::conv_to<arma::fmat>::from(c1B_d);
    }
  }

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  const int used_threads = std::min(std::max(1, n_threads), max_threads);
  omp_set_schedule((n_jobs >= (2 * used_threads)) ? omp_sched_static : omp_sched_dynamic, 1);
#pragma omp parallel for schedule(runtime) num_threads(used_threads) if(n_jobs > 1 && used_threads > 1)
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& Mref = M_refs[idx];
    const Rcpp::NumericMatrix& C1ref = C1_refs[idx];
    const Rcpp::NumericVector& pref = p_refs[idx];
    const arma::mat Mi(const_cast<double*>(Mref.begin()), static_cast<arma::uword>(Mref.nrow()), static_cast<arma::uword>(Mref.ncol()), false, true);
    const arma::mat C1i(const_cast<double*>(C1ref.begin()), static_cast<arma::uword>(C1ref.nrow()), static_cast<arma::uword>(C1ref.ncol()), false, true);
    const arma::vec pi(const_cast<double*>(pref.begin()), static_cast<arma::uword>(pref.size()), false, true);
    arma::mat init_i;
    arma::mat c1A_i;
    arma::mat c1Bt_i;
    const arma::mat* c1A_ptr = nullptr;
    const arma::mat* c1Bt_ptr = nullptr;
    const arma::fmat* c1A_ptr_f = nullptr;
    const arma::fmat* c1Bt_ptr_f = nullptr;
    const arma::fmat* M_ptr_f = nullptr;
    const arma::fmat* C1_ptr_f = nullptr;
    const arma::fvec* p_ptr_f = nullptr;
    if (use_mixed_kernel) {
      M_ptr_f = &(M_kernel_cache_f[idx]);
      C1_ptr_f = &(C1_kernel_cache_f[idx]);
      p_ptr_f = &(p_kernel_cache_f[idx]);
    }
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    }
    if (has_c1_cache[idx]) {
      if (use_mixed_kernel) {
        c1A_ptr_f = &(c1_A_cache_f[idx]);
        c1Bt_ptr_f = &(c1_Bt_cache_f[idx]);
      } else {
        const Rcpp::NumericMatrix& c1A_ref = c1_A_refs[idx];
        const Rcpp::NumericMatrix& c1Bt_ref = c1_Bt_refs[idx];
        c1A_i = arma::mat(const_cast<double*>(c1A_ref.begin()), static_cast<arma::uword>(c1A_ref.nrow()), static_cast<arma::uword>(c1A_ref.ncol()), false, true);
        c1Bt_i = arma::mat(const_cast<double*>(c1Bt_ref.begin()), static_cast<arma::uword>(c1Bt_ref.nrow()), static_cast<arma::uword>(c1Bt_ref.ncol()), false, true);
        c1A_ptr = &c1A_i;
        c1Bt_ptr = &c1Bt_i;
      }
    }
    const auto t0 = std::chrono::steady_clock::now();
    FgwEntropicCoreResult out = fgw_entropic_square_core(
      Mi, C1i, C2, pi, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, use_log_sinkhorn, use_mixed_precision, check_every, approx_rank,
      init_i, c1A_ptr, c1Bt_ptr, c1A_ptr_f, c1Bt_ptr_f,
      M_ptr_f, C1_ptr_f, p_ptr_f,
      c2_cache_d_ptr, c2_cache_f_ptr, c2_square_cache_d_ptr, c2_square_cache_f_ptr
    );
    const auto t1 = std::chrono::steady_clock::now();
    plans[idx] = std::move(out.plan);
    dists[idx] = out.fgw_dist;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
    lowrank_used[idx] = out.used_lowrank ? 1 : 0;
    c1_cache_used[idx] = out.used_c1_cache ? 1 : 0;
    c2_cache_used[idx] = out.used_c2_cache ? 1 : 0;
    square_cache_used[idx] = out.used_square_cache ? 1 : 0;
  }
#else
  const int used_threads = 1;
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& Mref = M_refs[idx];
    const Rcpp::NumericMatrix& C1ref = C1_refs[idx];
    const Rcpp::NumericVector& pref = p_refs[idx];
    const arma::mat Mi(const_cast<double*>(Mref.begin()), static_cast<arma::uword>(Mref.nrow()), static_cast<arma::uword>(Mref.ncol()), false, true);
    const arma::mat C1i(const_cast<double*>(C1ref.begin()), static_cast<arma::uword>(C1ref.nrow()), static_cast<arma::uword>(C1ref.ncol()), false, true);
    const arma::vec pi(const_cast<double*>(pref.begin()), static_cast<arma::uword>(pref.size()), false, true);
    arma::mat init_i;
    arma::mat c1A_i;
    arma::mat c1Bt_i;
    const arma::mat* c1A_ptr = nullptr;
    const arma::mat* c1Bt_ptr = nullptr;
    const arma::fmat* c1A_ptr_f = nullptr;
    const arma::fmat* c1Bt_ptr_f = nullptr;
    const arma::fmat* M_ptr_f = nullptr;
    const arma::fmat* C1_ptr_f = nullptr;
    const arma::fvec* p_ptr_f = nullptr;
    if (use_mixed_kernel) {
      M_ptr_f = &(M_kernel_cache_f[idx]);
      C1_ptr_f = &(C1_kernel_cache_f[idx]);
      p_ptr_f = &(p_kernel_cache_f[idx]);
    }
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    }
    if (has_c1_cache[idx]) {
      if (use_mixed_kernel) {
        c1A_ptr_f = &(c1_A_cache_f[idx]);
        c1Bt_ptr_f = &(c1_Bt_cache_f[idx]);
      } else {
        const Rcpp::NumericMatrix& c1A_ref = c1_A_refs[idx];
        const Rcpp::NumericMatrix& c1Bt_ref = c1_Bt_refs[idx];
        c1A_i = arma::mat(const_cast<double*>(c1A_ref.begin()), static_cast<arma::uword>(c1A_ref.nrow()), static_cast<arma::uword>(c1A_ref.ncol()), false, true);
        c1Bt_i = arma::mat(const_cast<double*>(c1Bt_ref.begin()), static_cast<arma::uword>(c1Bt_ref.nrow()), static_cast<arma::uword>(c1Bt_ref.ncol()), false, true);
        c1A_ptr = &c1A_i;
        c1Bt_ptr = &c1Bt_i;
      }
    }
    const auto t0 = std::chrono::steady_clock::now();
    FgwEntropicCoreResult out = fgw_entropic_square_core(
      Mi, C1i, C2, pi, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, use_log_sinkhorn, use_mixed_precision, check_every, approx_rank,
      init_i, c1A_ptr, c1Bt_ptr, c1A_ptr_f, c1Bt_ptr_f,
      M_ptr_f, C1_ptr_f, p_ptr_f,
      c2_cache_d_ptr, c2_cache_f_ptr, c2_square_cache_d_ptr, c2_square_cache_f_ptr
    );
    const auto t1 = std::chrono::steady_clock::now();
    plans[idx] = std::move(out.plan);
    dists[idx] = out.fgw_dist;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
    lowrank_used[idx] = out.used_lowrank ? 1 : 0;
    c1_cache_used[idx] = out.used_c1_cache ? 1 : 0;
    c2_cache_used[idx] = out.used_c2_cache ? 1 : 0;
    square_cache_used[idx] = out.used_square_cache ? 1 : 0;
  }
#endif

  int n_lowrank = 0;
  int n_c1_cache = 0;
  int n_c2_cache = 0;
  int n_square_cache = 0;
  for (int i = 0; i < n_jobs; ++i) {
    n_lowrank += lowrank_used[static_cast<std::size_t>(i)];
    n_c1_cache += c1_cache_used[static_cast<std::size_t>(i)];
    n_c2_cache += c2_cache_used[static_cast<std::size_t>(i)];
    n_square_cache += square_cache_used[static_cast<std::size_t>(i)];
  }

  Rcpp::List out_plans(n_jobs);
  for (int i = 0; i < n_jobs; ++i) {
    out_plans[i] = plans[static_cast<std::size_t>(i)];
  }

  return Rcpp::List::create(
    Rcpp::Named("plans") = out_plans,
    Rcpp::Named("fgw_dist") = dists,
    Rcpp::Named("iterations") = iters,
    Rcpp::Named("error") = errs,
    Rcpp::Named("used_threads") = used_threads,
    Rcpp::Named("kernel_ms") = kernel_ms,
    Rcpp::Named("lowrank_used") = lowrank_used,
    Rcpp::Named("c1_cache_used") = c1_cache_used,
    Rcpp::Named("c2_cache_used") = c2_cache_used,
    Rcpp::Named("square_cache_used") = square_cache_used,
    Rcpp::Named("n_lowrank") = n_lowrank,
    Rcpp::Named("n_c1_cache") = n_c1_cache,
    Rcpp::Named("n_c2_cache") = n_c2_cache,
    Rcpp::Named("n_square_cache") = n_square_cache
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_fgw_entropic_square_batch_features(
    const Rcpp::List& F1_list,
    const Rcpp::List& C1_list,
    const Rcpp::List& p_list,
    const arma::mat& F2,
    const arma::mat& C2,
    const arma::vec& q,
    bool use_euclidean,
    bool rescale01,
    double alpha,
    double epsilon,
    int max_iter,
    double tol,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    bool use_ppa,
    bool use_log_sinkhorn,
    bool use_mixed_precision,
    int check_every,
    const Rcpp::List& init_plan_list,
    const Rcpp::List& c1_A_scaled_list,
    const Rcpp::List& c1_Bt_list,
    int approx_rank,
    int n_threads) {
  const int n_jobs = static_cast<int>(C1_list.size());
  if (n_jobs <= 0) {
    Rcpp::stop("`C1_list` must be non-empty.");
  }
  if (F1_list.size() != C1_list.size() || p_list.size() != C1_list.size()) {
    Rcpp::stop("`F1_list`, `C1_list`, and `p_list` must have the same length.");
  }
  if (init_plan_list.size() > 0 && init_plan_list.size() != C1_list.size()) {
    Rcpp::stop("`init_plan_list` must be empty or the same length as `C1_list`.");
  }
  if ((c1_A_scaled_list.size() > 0 && c1_A_scaled_list.size() != C1_list.size()) ||
      (c1_Bt_list.size() > 0 && c1_Bt_list.size() != C1_list.size())) {
    Rcpp::stop("`c1_A_scaled_list` and `c1_Bt_list` must be empty or the same length as `C1_list`.");
  }
  if ((c1_A_scaled_list.size() > 0) != (c1_Bt_list.size() > 0)) {
    Rcpp::stop("Provide both `c1_A_scaled_list` and `c1_Bt_list` together, or neither.");
  }
  if (q.n_elem != C2.n_rows || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C2` must be square and compatible with `q`.");
  }
  if (F2.n_rows != C2.n_rows) {
    Rcpp::stop("`F2` must have nrow = nrow(C2).");
  }
  if (check_every <= 0) {
    check_every = 1;
  }
  if (approx_rank < 0) {
    approx_rank = 0;
  }
  if (n_threads <= 0) {
    n_threads = 1;
  }

  std::vector<Rcpp::NumericMatrix> F1_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> C1_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericVector> p_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> init_refs(static_cast<std::size_t>(n_jobs));
  std::vector<unsigned char> has_init(static_cast<std::size_t>(n_jobs), 0);
  std::vector<Rcpp::NumericMatrix> c1_A_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> c1_Bt_refs(static_cast<std::size_t>(n_jobs));
  std::vector<unsigned char> has_c1_cache(static_cast<std::size_t>(n_jobs), 0);
  for (int i = 0; i < n_jobs; ++i) {
    F1_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(F1_list[i]);
    C1_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(C1_list[i]);
    p_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericVector>(p_list[i]);
    if (init_plan_list.size() == C1_list.size()) {
      SEXP init_i = init_plan_list[i];
      if (init_i != R_NilValue) {
        init_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(init_i);
        has_init[static_cast<std::size_t>(i)] = 1;
      }
    }
    if (c1_A_scaled_list.size() == C1_list.size()) {
      SEXP c1a_i = c1_A_scaled_list[i];
      SEXP c1b_i = c1_Bt_list[i];
      if (c1a_i != R_NilValue && c1b_i != R_NilValue) {
        c1_A_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(c1a_i);
        c1_Bt_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(c1b_i);
        has_c1_cache[static_cast<std::size_t>(i)] = 1;
      }
    }
    const Rcpp::NumericMatrix& F1i = F1_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericMatrix& C1i = C1_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericVector& pi = p_refs[static_cast<std::size_t>(i)];
    if (C1i.nrow() != C1i.ncol()) {
      Rcpp::stop("Every `C1_list[[i]]` must be square.");
    }
    if (F1i.nrow() != C1i.nrow()) {
      Rcpp::stop("Each `F1_list[[i]]` must have nrow = nrow(C1_list[[i]]).");
    }
    if (F1i.ncol() != static_cast<int>(F2.n_cols)) {
      Rcpp::stop("Each `F1_list[[i]]` must have ncol = ncol(F2).");
    }
    if (static_cast<int>(pi.size()) != C1i.nrow()) {
      Rcpp::stop("Each `p_list[[i]]` must match `nrow(C1_list[[i]])`.");
    }
    if (has_c1_cache[static_cast<std::size_t>(i)] == 1) {
      const Rcpp::NumericMatrix& c1A = c1_A_refs[static_cast<std::size_t>(i)];
      const Rcpp::NumericMatrix& c1B = c1_Bt_refs[static_cast<std::size_t>(i)];
      if (c1A.nrow() != C1i.nrow() || c1B.ncol() != C1i.ncol() || c1A.ncol() != c1B.nrow() || c1A.ncol() <= 0) {
        Rcpp::stop("Each cached C1 low-rank pair must satisfy A(n x r), Bt(r x n).");
      }
    }
  }

  std::vector<arma::mat> plans(static_cast<std::size_t>(n_jobs));
  std::vector<double> dists(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<double> errs(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<int> iters(static_cast<std::size_t>(n_jobs), 0);
  std::vector<double> kernel_ms(static_cast<std::size_t>(n_jobs), 0.0);
  std::vector<double> feature_ms(static_cast<std::size_t>(n_jobs), 0.0);
  std::vector<double> solve_ms(static_cast<std::size_t>(n_jobs), 0.0);
  std::vector<int> lowrank_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> c1_cache_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> c2_cache_used(static_cast<std::size_t>(n_jobs), 0);
  std::vector<int> square_cache_used(static_cast<std::size_t>(n_jobs), 0);
  SquareC2CacheD c2_square_cache_d;
  SquareC2CacheF c2_square_cache_f;
  LowRankC2CacheD c2_cache_d;
  LowRankC2CacheF c2_cache_f;
  const bool use_mixed_kernel = use_mixed_precision && !use_log_sinkhorn;
  if (symmetric) {
    if (use_mixed_kernel) {
      const arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
      const arma::fvec qf = arma::conv_to<arma::fvec>::from(q);
      c2_square_cache_f = build_square_c2_cache_f(C2f, qf);
    } else {
      c2_square_cache_d = build_square_c2_cache_d(C2, q);
    }
  }
  const bool want_lowrank = symmetric && (approx_rank > 0);
  if (want_lowrank) {
    if (use_mixed_kernel) {
      if (c2_square_cache_f.valid) {
        c2_cache_f = build_lowrank_c2_cache_f(c2_square_cache_f.C2, approx_rank);
      } else {
        const arma::fmat C2f = arma::conv_to<arma::fmat>::from(C2);
        c2_cache_f = build_lowrank_c2_cache_f(C2f, approx_rank);
      }
    } else {
      c2_cache_d = build_lowrank_c2_cache_d(C2, approx_rank);
    }
  }
  const LowRankC2CacheD* c2_cache_d_ptr = c2_cache_d.valid ? &c2_cache_d : nullptr;
  const LowRankC2CacheF* c2_cache_f_ptr = c2_cache_f.valid ? &c2_cache_f : nullptr;
  const SquareC2CacheD* c2_square_cache_d_ptr = c2_square_cache_d.valid ? &c2_square_cache_d : nullptr;
  const SquareC2CacheF* c2_square_cache_f_ptr = c2_square_cache_f.valid ? &c2_square_cache_f : nullptr;
  std::vector<arma::fmat> c1_A_cache_f;
  std::vector<arma::fmat> c1_Bt_cache_f;
  if (use_mixed_kernel && c1_A_scaled_list.size() == C1_list.size()) {
    c1_A_cache_f.resize(static_cast<std::size_t>(n_jobs));
    c1_Bt_cache_f.resize(static_cast<std::size_t>(n_jobs));
    for (int i = 0; i < n_jobs; ++i) {
      const std::size_t idx = static_cast<std::size_t>(i);
      if (has_c1_cache[idx] == 0) {
        continue;
      }
      const Rcpp::NumericMatrix& c1A = c1_A_refs[idx];
      const Rcpp::NumericMatrix& c1B = c1_Bt_refs[idx];
      const arma::mat c1A_d(const_cast<double*>(c1A.begin()), static_cast<arma::uword>(c1A.nrow()), static_cast<arma::uword>(c1A.ncol()), false, true);
      const arma::mat c1B_d(const_cast<double*>(c1B.begin()), static_cast<arma::uword>(c1B.nrow()), static_cast<arma::uword>(c1B.ncol()), false, true);
      c1_A_cache_f[idx] = arma::conv_to<arma::fmat>::from(c1A_d);
      c1_Bt_cache_f[idx] = arma::conv_to<arma::fmat>::from(c1B_d);
    }
  }

  const arma::vec y2 = arma::sum(arma::square(F2), 1);

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  const int used_threads = std::min(std::max(1, n_threads), max_threads);
  omp_set_schedule((n_jobs >= (2 * used_threads)) ? omp_sched_static : omp_sched_dynamic, 1);
#pragma omp parallel for schedule(runtime) num_threads(used_threads) if(n_jobs > 1 && used_threads > 1)
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& F1ref = F1_refs[idx];
    const Rcpp::NumericMatrix& C1ref = C1_refs[idx];
    const Rcpp::NumericVector& pref = p_refs[idx];
    const arma::mat F1i(const_cast<double*>(F1ref.begin()), static_cast<arma::uword>(F1ref.nrow()), static_cast<arma::uword>(F1ref.ncol()), false, true);
    const arma::mat C1i(const_cast<double*>(C1ref.begin()), static_cast<arma::uword>(C1ref.nrow()), static_cast<arma::uword>(C1ref.ncol()), false, true);
    const arma::vec pi(const_cast<double*>(pref.begin()), static_cast<arma::uword>(pref.size()), false, true);
    const auto t0 = std::chrono::steady_clock::now();
    arma::mat Mi = cross_feature_cost_matrix(F1i, F2, y2, use_euclidean, rescale01);
    arma::mat init_i;
    arma::mat c1A_i;
    arma::mat c1Bt_i;
    const arma::mat* c1A_ptr = nullptr;
    const arma::mat* c1Bt_ptr = nullptr;
    const arma::fmat* c1A_ptr_f = nullptr;
    const arma::fmat* c1Bt_ptr_f = nullptr;
    const arma::fmat* M_ptr_f = nullptr;
    const arma::fmat* C1_ptr_f = nullptr;
    const arma::fvec* p_ptr_f = nullptr;
    arma::fmat Mi_f;
    arma::fmat C1i_f;
    arma::fvec pi_f;
    if (use_mixed_kernel) {
      Mi_f = arma::conv_to<arma::fmat>::from(Mi);
      C1i_f = arma::conv_to<arma::fmat>::from(C1i);
      pi_f = arma::conv_to<arma::fvec>::from(pi);
      M_ptr_f = &Mi_f;
      C1_ptr_f = &C1i_f;
      p_ptr_f = &pi_f;
    }
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    }
    if (has_c1_cache[idx]) {
      if (use_mixed_kernel) {
        c1A_ptr_f = &(c1_A_cache_f[idx]);
        c1Bt_ptr_f = &(c1_Bt_cache_f[idx]);
      } else {
        const Rcpp::NumericMatrix& c1A_ref = c1_A_refs[idx];
        const Rcpp::NumericMatrix& c1Bt_ref = c1_Bt_refs[idx];
        c1A_i = arma::mat(const_cast<double*>(c1A_ref.begin()), static_cast<arma::uword>(c1A_ref.nrow()), static_cast<arma::uword>(c1A_ref.ncol()), false, true);
        c1Bt_i = arma::mat(const_cast<double*>(c1Bt_ref.begin()), static_cast<arma::uword>(c1Bt_ref.nrow()), static_cast<arma::uword>(c1Bt_ref.ncol()), false, true);
        c1A_ptr = &c1A_i;
        c1Bt_ptr = &c1Bt_i;
      }
    }
    const auto tf = std::chrono::steady_clock::now();
    FgwEntropicCoreResult out = fgw_entropic_square_core(
      Mi, C1i, C2, pi, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, use_log_sinkhorn, use_mixed_precision, check_every, approx_rank,
      init_i, c1A_ptr, c1Bt_ptr, c1A_ptr_f, c1Bt_ptr_f,
      M_ptr_f, C1_ptr_f, p_ptr_f,
      c2_cache_d_ptr, c2_cache_f_ptr, c2_square_cache_d_ptr, c2_square_cache_f_ptr
    );
    const auto t1 = std::chrono::steady_clock::now();
    plans[idx] = std::move(out.plan);
    dists[idx] = out.fgw_dist;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
    feature_ms[idx] = std::chrono::duration<double, std::milli>(tf - t0).count();
    solve_ms[idx] = std::chrono::duration<double, std::milli>(t1 - tf).count();
    lowrank_used[idx] = out.used_lowrank ? 1 : 0;
    c1_cache_used[idx] = out.used_c1_cache ? 1 : 0;
    c2_cache_used[idx] = out.used_c2_cache ? 1 : 0;
    square_cache_used[idx] = out.used_square_cache ? 1 : 0;
  }
#else
  const int used_threads = 1;
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& F1ref = F1_refs[idx];
    const Rcpp::NumericMatrix& C1ref = C1_refs[idx];
    const Rcpp::NumericVector& pref = p_refs[idx];
    const arma::mat F1i(const_cast<double*>(F1ref.begin()), static_cast<arma::uword>(F1ref.nrow()), static_cast<arma::uword>(F1ref.ncol()), false, true);
    const arma::mat C1i(const_cast<double*>(C1ref.begin()), static_cast<arma::uword>(C1ref.nrow()), static_cast<arma::uword>(C1ref.ncol()), false, true);
    const arma::vec pi(const_cast<double*>(pref.begin()), static_cast<arma::uword>(pref.size()), false, true);
    const auto t0 = std::chrono::steady_clock::now();
    arma::mat Mi = cross_feature_cost_matrix(F1i, F2, y2, use_euclidean, rescale01);
    arma::mat init_i;
    arma::mat c1A_i;
    arma::mat c1Bt_i;
    const arma::mat* c1A_ptr = nullptr;
    const arma::mat* c1Bt_ptr = nullptr;
    const arma::fmat* c1A_ptr_f = nullptr;
    const arma::fmat* c1Bt_ptr_f = nullptr;
    const arma::fmat* M_ptr_f = nullptr;
    const arma::fmat* C1_ptr_f = nullptr;
    const arma::fvec* p_ptr_f = nullptr;
    arma::fmat Mi_f;
    arma::fmat C1i_f;
    arma::fvec pi_f;
    if (use_mixed_kernel) {
      Mi_f = arma::conv_to<arma::fmat>::from(Mi);
      C1i_f = arma::conv_to<arma::fmat>::from(C1i);
      pi_f = arma::conv_to<arma::fvec>::from(pi);
      M_ptr_f = &Mi_f;
      C1_ptr_f = &C1i_f;
      p_ptr_f = &pi_f;
    }
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    }
    if (has_c1_cache[idx]) {
      if (use_mixed_kernel) {
        c1A_ptr_f = &(c1_A_cache_f[idx]);
        c1Bt_ptr_f = &(c1_Bt_cache_f[idx]);
      } else {
        const Rcpp::NumericMatrix& c1A_ref = c1_A_refs[idx];
        const Rcpp::NumericMatrix& c1Bt_ref = c1_Bt_refs[idx];
        c1A_i = arma::mat(const_cast<double*>(c1A_ref.begin()), static_cast<arma::uword>(c1A_ref.nrow()), static_cast<arma::uword>(c1A_ref.ncol()), false, true);
        c1Bt_i = arma::mat(const_cast<double*>(c1Bt_ref.begin()), static_cast<arma::uword>(c1Bt_ref.nrow()), static_cast<arma::uword>(c1Bt_ref.ncol()), false, true);
        c1A_ptr = &c1A_i;
        c1Bt_ptr = &c1Bt_i;
      }
    }
    const auto tf = std::chrono::steady_clock::now();
    FgwEntropicCoreResult out = fgw_entropic_square_core(
      Mi, C1i, C2, pi, q,
      alpha, epsilon, max_iter, tol, sinkhorn_max_iter, sinkhorn_tol,
      symmetric, use_ppa, use_log_sinkhorn, use_mixed_precision, check_every, approx_rank,
      init_i, c1A_ptr, c1Bt_ptr, c1A_ptr_f, c1Bt_ptr_f,
      M_ptr_f, C1_ptr_f, p_ptr_f,
      c2_cache_d_ptr, c2_cache_f_ptr, c2_square_cache_d_ptr, c2_square_cache_f_ptr
    );
    const auto t1 = std::chrono::steady_clock::now();
    plans[idx] = std::move(out.plan);
    dists[idx] = out.fgw_dist;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
    feature_ms[idx] = std::chrono::duration<double, std::milli>(tf - t0).count();
    solve_ms[idx] = std::chrono::duration<double, std::milli>(t1 - tf).count();
    lowrank_used[idx] = out.used_lowrank ? 1 : 0;
    c1_cache_used[idx] = out.used_c1_cache ? 1 : 0;
    c2_cache_used[idx] = out.used_c2_cache ? 1 : 0;
    square_cache_used[idx] = out.used_square_cache ? 1 : 0;
  }
#endif

  int n_lowrank = 0;
  int n_c1_cache = 0;
  int n_c2_cache = 0;
  int n_square_cache = 0;
  for (int i = 0; i < n_jobs; ++i) {
    n_lowrank += lowrank_used[static_cast<std::size_t>(i)];
    n_c1_cache += c1_cache_used[static_cast<std::size_t>(i)];
    n_c2_cache += c2_cache_used[static_cast<std::size_t>(i)];
    n_square_cache += square_cache_used[static_cast<std::size_t>(i)];
  }

  Rcpp::List out_plans(n_jobs);
  for (int i = 0; i < n_jobs; ++i) {
    out_plans[i] = plans[static_cast<std::size_t>(i)];
  }

  return Rcpp::List::create(
    Rcpp::Named("plans") = out_plans,
    Rcpp::Named("fgw_dist") = dists,
    Rcpp::Named("iterations") = iters,
    Rcpp::Named("error") = errs,
    Rcpp::Named("used_threads") = used_threads,
    Rcpp::Named("kernel_ms") = kernel_ms,
    Rcpp::Named("feature_ms") = feature_ms,
    Rcpp::Named("solve_ms") = solve_ms,
    Rcpp::Named("lowrank_used") = lowrank_used,
    Rcpp::Named("c1_cache_used") = c1_cache_used,
    Rcpp::Named("c2_cache_used") = c2_cache_used,
    Rcpp::Named("square_cache_used") = square_cache_used,
    Rcpp::Named("n_lowrank") = n_lowrank,
    Rcpp::Named("n_c1_cache") = n_c1_cache,
    Rcpp::Named("n_c2_cache") = n_c2_cache,
    Rcpp::Named("n_square_cache") = n_square_cache
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_feature_cost_batch(
    const Rcpp::List& F1_list,
    const arma::mat& F2,
    bool use_euclidean,
    bool rescale01,
    int n_threads) {
  const int n_jobs = static_cast<int>(F1_list.size());
  if (n_jobs <= 0) {
    return Rcpp::List();
  }
  if (n_threads <= 0) {
    n_threads = 1;
  }

  std::vector<Rcpp::NumericMatrix> F1_refs(static_cast<std::size_t>(n_jobs));
  for (int i = 0; i < n_jobs; ++i) {
    F1_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(F1_list[i]);
    const Rcpp::NumericMatrix& F1i = F1_refs[static_cast<std::size_t>(i)];
    if (F1i.ncol() != static_cast<int>(F2.n_cols)) {
      Rcpp::stop("Each `F1_list[[i]]` must have ncol = ncol(`F2`).");
    }
  }

  const arma::vec y2 = arma::sum(arma::square(F2), 1);
  std::vector<arma::mat> out(static_cast<std::size_t>(n_jobs));

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  const int used_threads = std::min(std::max(1, n_threads), max_threads);
  omp_set_schedule((n_jobs >= (2 * used_threads)) ? omp_sched_static : omp_sched_dynamic, 1);
#pragma omp parallel for schedule(runtime) num_threads(used_threads) if(n_jobs > 1 && used_threads > 1)
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& F1ref = F1_refs[idx];
    const arma::mat F1(
      const_cast<double*>(F1ref.begin()),
      static_cast<arma::uword>(F1ref.nrow()),
      static_cast<arma::uword>(F1ref.ncol()),
      false,
      true
    );
    out[idx] = cross_feature_cost_matrix(F1, F2, y2, use_euclidean, rescale01);
  }
#else
  const int used_threads = 1;
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& F1ref = F1_refs[idx];
    const arma::mat F1(
      const_cast<double*>(F1ref.begin()),
      static_cast<arma::uword>(F1ref.nrow()),
      static_cast<arma::uword>(F1ref.ncol()),
      false,
      true
    );
    out[idx] = cross_feature_cost_matrix(F1, F2, y2, use_euclidean, rescale01);
  }
#endif

  Rcpp::List out_list(n_jobs);
  for (int i = 0; i < n_jobs; ++i) {
    out_list[i] = std::move(out[static_cast<std::size_t>(i)]);
  }
  out_list.attr("used_threads") = used_threads;
  return out_list;
}

// [[Rcpp::export]]
Rcpp::List cpp_fgw_exact_cg_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double alpha,
    bool symmetric,
    int max_iter,
    double tol_rel,
    double tol_abs,
    int lp_max_iter,
    double lp_tol,
    const arma::mat& init_plan) {
  arma::mat constC, hC1, hC2;
  init_matrices_square(C1, C2, p, q, constC, hC1, hC2);
  arma::mat constCt, hC1t, hC2t;
  if (!symmetric) {
    init_matrices_square(C1.t(), C2.t(), p, q, constCt, hC1t, hC2t);
  }

  const arma::mat M_lin = (1.0 - alpha) * M;
  arma::mat G;
  if (init_plan.n_rows == 0 || init_plan.n_cols == 0) {
    G = p * q.t();
  } else {
    if (init_plan.n_rows != p.n_elem || init_plan.n_cols != q.n_elem) {
      Rcpp::stop("`init_plan` must have shape length(p) x length(q).");
    }
    G = init_plan;
  }
  arma::mat scratch(C1.n_rows, C2.n_rows);
  arma::mat Acur(C1.n_rows, C2.n_rows);
  dgemm_nn(hC1, G, scratch);
  dgemm_nt(scratch, hC2, Acur);
  arma::mat tens = constC - Acur;

  arma::mat Acurt;
  arma::mat tenst;
  if (!symmetric) {
    Acurt.set_size(C1.n_rows, C2.n_rows);
    dgemm_nn(hC1t, G, scratch);
    dgemm_nt(scratch, hC2t, Acurt);
    tenst = constCt - Acurt;
  }

  const double gw_term = symmetric ? arma::accu(tens % G) : 0.5 * (arma::accu(tens % G) + arma::accu(tenst % G));
  double cost_G = arma::accu(M_lin % G) + alpha * gw_term;
  std::vector<double> loss_trace;
  loss_trace.reserve(static_cast<std::size_t>(max_iter) + 1u);
  loss_trace.push_back(cost_G);

  double abs_delta = std::numeric_limits<double>::quiet_NaN();
  double rel_delta = std::numeric_limits<double>::quiet_NaN();
  int it = 0;
  bool lp_ok = true;
  int inner_iterations = 0;

  arma::mat gG(C1.n_rows, C2.n_rows);
  arma::mat Mi(C1.n_rows, C2.n_rows);
  arma::mat deltaG(C1.n_rows, C2.n_rows);
  arma::mat dot(C1.n_rows, C2.n_rows);
  arma::mat dot_t;
  if (!symmetric) {
    dot_t.set_size(C1.n_rows, C2.n_rows);
  }

  for (int k = 0; k < max_iter; ++k) {
    if (symmetric) {
      gG = 2.0 * tens;
    } else {
      gG = tens + tenst;
    }
    Mi = M_lin + alpha * gG;

    const TransportSimplexResult dir = transport_simplex_solve(Mi, p, q, lp_max_iter, lp_tol);
    lp_ok = lp_ok && dir.converged;
    inner_iterations += dir.iterations;
    arma::mat Gc = dir.plan;
    deltaG = Gc - G;

    dgemm_nn(hC1, deltaG, scratch);
    dgemm_nt(scratch, hC2, dot);
    const double a_ls = -alpha * arma::accu(dot % deltaG);
    double b_ls;
    if (symmetric) {
      b_ls = arma::accu(M_lin % deltaG) - 2.0 * alpha * arma::accu(dot % G);
    } else {
      b_ls = arma::accu(M_lin % deltaG) - alpha * (arma::accu(dot % G) + arma::accu(Acur % deltaG));
    }

    double step = solve_1d_linesearch_quad(a_ls, b_ls);
    step = std::min(1.0, std::max(0.0, step));
    const double new_cost = cost_G + a_ls * step * step + b_ls * step;

    G += step * deltaG;
    Acur += step * dot;
    tens = constC - Acur;
    if (!symmetric) {
      dgemm_nn(hC1t, deltaG, scratch);
      dgemm_nt(scratch, hC2t, dot_t);
      Acurt += step * dot_t;
      tenst = constCt - Acurt;
    }

    abs_delta = std::abs(new_cost - cost_G);
    rel_delta = (new_cost != 0.0) ? (abs_delta / std::abs(new_cost)) : std::numeric_limits<double>::quiet_NaN();
    cost_G = new_cost;
    loss_trace.push_back(cost_G);
    it = k + 1;

    if ((std::isfinite(rel_delta) && rel_delta < tol_rel) || (std::isfinite(abs_delta) && abs_delta < tol_abs)) {
      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("plan") = G,
    Rcpp::Named("fgw_dist") = cost_G,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = abs_delta,
    Rcpp::Named("rel_error") = rel_delta,
    Rcpp::Named("loss_trace") = loss_trace,
    Rcpp::Named("lp_ok") = lp_ok,
    Rcpp::Named("inner_iterations") = inner_iterations
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_fugw_kl_square(
    const arma::mat& Cx,
    const arma::mat& Cy,
    const arma::vec& wx,
    const arma::vec& wy,
    const arma::vec& reg_marginals,
    double epsilon,
    double alpha,
    const arma::mat& M,
    const arma::mat& init_pi,
    int max_iter,
    double tol,
    int max_iter_ot,
    double tol_ot,
    bool rescale_plan,
    int check_every,
    bool use_mixed_precision) {
  if (check_every <= 0) {
    check_every = 1;
  }

  const double rho_x = reg_marginals(0);
  const double rho_y = reg_marginals(1);
  const arma::mat M_samp = (alpha / 2.0) * M;
  const arma::mat M_feat = (alpha / 2.0) * M;
  const arma::mat Cx_sqr = Cx % Cx;
  const arma::mat Cy_sqr = Cy % Cy;
  const arma::mat Cy_t = Cy.t();
  const arma::mat wxy = wx * wy.t();
  const arma::vec log_wx = arma::log(wx + kTiny);
  const arma::vec log_wy = arma::log(wy + kTiny);

  Rcpp::NumericVector err_trace;
  std::vector<int> inner_iters_feat;
  std::vector<int> inner_iters_samp;
  std::vector<int> inner_warm_feat;
  std::vector<int> inner_warm_samp;
  std::vector<int> inner_fallback_feat;
  std::vector<int> inner_fallback_samp;
  inner_iters_feat.reserve(static_cast<std::size_t>(std::max(0, max_iter)));
  inner_iters_samp.reserve(static_cast<std::size_t>(std::max(0, max_iter)));
  inner_warm_feat.reserve(static_cast<std::size_t>(std::max(0, max_iter)));
  inner_warm_samp.reserve(static_cast<std::size_t>(std::max(0, max_iter)));
  inner_fallback_feat.reserve(static_cast<std::size_t>(std::max(0, max_iter)));
  inner_fallback_samp.reserve(static_cast<std::size_t>(std::max(0, max_iter)));

  if (!use_mixed_precision) {
    arma::mat pi_samp = init_pi;
    arma::mat pi_feat = init_pi;
    SinkhornUnbalancedWorkspace ws_samp;
    SinkhornUnbalancedWorkspace ws_feat;
    arma::mat uot_cost(Cx.n_rows, Cy.n_rows);
    arma::mat scratch(Cx.n_rows, Cy.n_rows);
    arma::vec pi1(wx.n_elem);
    arma::vec pi2(wy.n_elem);
    arma::vec A(wx.n_elem);
    arma::vec B(wy.n_elem);

    double err = std::numeric_limits<double>::infinity();
    double outer_err_prev = std::numeric_limits<double>::infinity();
    int it = 0;
    for (; it < max_iter; ++it) {
      const bool do_check = ((it + 1) % check_every == 0);
      const double tol_ot_eff = fugw_effective_inner_tol(tol_ot, outer_err_prev, tol);

      {
        const double mass = arma::accu(pi_samp);
        uot_cost_matrix_kl_joint_inplace(
          Cx_sqr, Cy_sqr, Cx, Cy_t, M_feat, pi_samp, wx, wy, log_wx, log_wy, rho_x, rho_y, epsilon,
          uot_cost, scratch, pi1, pi2, A, B
        );
        SinkhornUnbalancedResult su = sinkhorn_unbalanced_kl(
          uot_cost,
          wx,
          wy,
          wxy,
          rho_x * mass,
          rho_y * mass,
          epsilon * mass,
          max_iter_ot,
          tol_ot_eff,
          pi_feat,
          ws_feat,
          fugw_enable_warm_start()
        );
        inner_iters_feat.push_back(su.iters);
        inner_warm_feat.push_back(su.warm_started ? 1 : 0);
        inner_fallback_feat.push_back(su.warm_fallback ? 1 : 0);
        pi_feat = std::move(su.plan);
        if (rescale_plan) {
          const double mass_feat = arma::accu(pi_feat);
          if (mass_feat > 0.0) {
            pi_feat *= std::sqrt(mass / mass_feat);
          }
        }
      }

      {
        const double mass = arma::accu(pi_feat);
        uot_cost_matrix_kl_joint_inplace(
          Cx_sqr, Cy_sqr, Cx, Cy_t, M_samp, pi_feat, wx, wy, log_wx, log_wy, rho_x, rho_y, epsilon,
          uot_cost, scratch, pi1, pi2, A, B
        );
        SinkhornUnbalancedResult su = sinkhorn_unbalanced_kl(
          uot_cost,
          wx,
          wy,
          wxy,
          rho_x * mass,
          rho_y * mass,
          epsilon * mass,
          max_iter_ot,
          tol_ot_eff,
          pi_samp,
          ws_samp,
          fugw_enable_warm_start()
        );
        inner_iters_samp.push_back(su.iters);
        inner_warm_samp.push_back(su.warm_started ? 1 : 0);
        inner_fallback_samp.push_back(su.warm_fallback ? 1 : 0);
        arma::mat pi_next = std::move(su.plan);
        if (rescale_plan) {
          const double mass_samp = arma::accu(pi_next);
          if (mass_samp > 0.0) {
            pi_next *= std::sqrt(mass / mass_samp);
          }
        }
        if (do_check) {
          err = arma::accu(arma::abs(pi_next - pi_samp));
        }
        pi_samp = std::move(pi_next);
      }

      if (do_check) {
        err_trace.push_back(err);
        outer_err_prev = err;
        if (err < tol) {
          ++it;
          break;
        }
      }
    }

    const auto costs = fused_unbalanced_cost_square_joint(
      Cx, Cy, Cx_sqr, Cy_sqr, M_samp, M_feat, wx, wy, wxy, pi_samp, pi_feat, rho_x, rho_y, epsilon
    );
    const int inner_total = std::accumulate(inner_iters_feat.begin(), inner_iters_feat.end(), 0) +
      std::accumulate(inner_iters_samp.begin(), inner_iters_samp.end(), 0);

    return Rcpp::List::create(
      Rcpp::Named("pi_samp") = pi_samp,
      Rcpp::Named("pi_feat") = pi_feat,
      Rcpp::Named("linear_cost") = costs.first,
      Rcpp::Named("fugw_cost") = costs.second,
      Rcpp::Named("iterations") = it,
      Rcpp::Named("error") = err,
      Rcpp::Named("err_trace") = err_trace,
      Rcpp::Named("inner_iters_feat") = Rcpp::IntegerVector(inner_iters_feat.begin(), inner_iters_feat.end()),
      Rcpp::Named("inner_iters_samp") = Rcpp::IntegerVector(inner_iters_samp.begin(), inner_iters_samp.end()),
      Rcpp::Named("inner_warm_feat") = Rcpp::LogicalVector(inner_warm_feat.begin(), inner_warm_feat.end()),
      Rcpp::Named("inner_warm_samp") = Rcpp::LogicalVector(inner_warm_samp.begin(), inner_warm_samp.end()),
      Rcpp::Named("inner_warm_fallback_feat") = Rcpp::LogicalVector(inner_fallback_feat.begin(), inner_fallback_feat.end()),
      Rcpp::Named("inner_warm_fallback_samp") = Rcpp::LogicalVector(inner_fallback_samp.begin(), inner_fallback_samp.end()),
      Rcpp::Named("inner_iters_total") = inner_total
    );
  }

  const float rho_x_f = static_cast<float>(rho_x);
  const float rho_y_f = static_cast<float>(rho_y);
  const float epsilon_f = static_cast<float>(epsilon);
  const arma::fmat Cxf = arma::conv_to<arma::fmat>::from(Cx);
  const arma::fmat Cyf = arma::conv_to<arma::fmat>::from(Cy);
  const arma::fvec wxf = arma::conv_to<arma::fvec>::from(wx);
  const arma::fvec wyf = arma::conv_to<arma::fvec>::from(wy);
  const arma::fmat Mf = arma::conv_to<arma::fmat>::from(M);
  const arma::fmat M_samp_f = (static_cast<float>(alpha) / 2.0f) * Mf;
  const arma::fmat M_feat_f = (static_cast<float>(alpha) / 2.0f) * Mf;
  const arma::fmat Cx_sqr_f = Cxf % Cxf;
  const arma::fmat Cy_sqr_f = Cyf % Cyf;
  const arma::fmat Cy_t_f = Cyf.t();
  const arma::fmat wxy_f = wxf * wyf.t();
  const arma::fvec log_wxf = arma::log(wxf + kTinyF);
  const arma::fvec log_wyf = arma::log(wyf + kTinyF);

  arma::fmat pi_samp_f = arma::conv_to<arma::fmat>::from(init_pi);
  arma::fmat pi_feat_f = arma::conv_to<arma::fmat>::from(init_pi);
  SinkhornUnbalancedWorkspaceF ws_samp_f;
  SinkhornUnbalancedWorkspaceF ws_feat_f;
  arma::fmat uot_cost_f(Cx.n_rows, Cy.n_rows);
  arma::fmat scratch_f(Cx.n_rows, Cy.n_rows);
  arma::fvec pi1_f(wx.n_elem);
  arma::fvec pi2_f(wy.n_elem);
  arma::fvec A_f(wx.n_elem);
  arma::fvec B_f(wy.n_elem);

  double err = std::numeric_limits<double>::infinity();
  double outer_err_prev = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    const bool do_check = ((it + 1) % check_every == 0);
    const float tol_ot_eff_f = static_cast<float>(fugw_effective_inner_tol(tol_ot, outer_err_prev, tol));

    {
      const float mass = arma::accu(pi_samp_f);
      uot_cost_matrix_kl_joint_inplace_f(
        Cx_sqr_f, Cy_sqr_f, Cxf, Cy_t_f, M_feat_f, pi_samp_f, wxf, wyf, log_wxf, log_wyf, rho_x_f, rho_y_f, epsilon_f,
        uot_cost_f, scratch_f, pi1_f, pi2_f, A_f, B_f
      );
      SinkhornUnbalancedResultF su = sinkhorn_unbalanced_kl_f(
        uot_cost_f,
        wxf,
        wyf,
        wxy_f,
        rho_x_f * mass,
        rho_y_f * mass,
        epsilon_f * mass,
        max_iter_ot,
        tol_ot_eff_f,
        pi_feat_f,
        ws_feat_f,
        fugw_enable_warm_start()
      );
      inner_iters_feat.push_back(su.iters);
      inner_warm_feat.push_back(su.warm_started ? 1 : 0);
      inner_fallback_feat.push_back(su.warm_fallback ? 1 : 0);
      pi_feat_f = std::move(su.plan);
      if (rescale_plan) {
        const float mass_feat = arma::accu(pi_feat_f);
        if (mass_feat > 0.0f) {
          pi_feat_f *= std::sqrt(mass / mass_feat);
        }
      }
    }

    {
      const float mass = arma::accu(pi_feat_f);
      uot_cost_matrix_kl_joint_inplace_f(
        Cx_sqr_f, Cy_sqr_f, Cxf, Cy_t_f, M_samp_f, pi_feat_f, wxf, wyf, log_wxf, log_wyf, rho_x_f, rho_y_f, epsilon_f,
        uot_cost_f, scratch_f, pi1_f, pi2_f, A_f, B_f
      );
      SinkhornUnbalancedResultF su = sinkhorn_unbalanced_kl_f(
        uot_cost_f,
        wxf,
        wyf,
        wxy_f,
        rho_x_f * mass,
        rho_y_f * mass,
        epsilon_f * mass,
        max_iter_ot,
        tol_ot_eff_f,
        pi_samp_f,
        ws_samp_f,
        fugw_enable_warm_start()
      );
      inner_iters_samp.push_back(su.iters);
      inner_warm_samp.push_back(su.warm_started ? 1 : 0);
      inner_fallback_samp.push_back(su.warm_fallback ? 1 : 0);
      arma::fmat pi_next = std::move(su.plan);
      if (rescale_plan) {
        const float mass_samp = arma::accu(pi_next);
        if (mass_samp > 0.0f) {
          pi_next *= std::sqrt(mass / mass_samp);
        }
      }
      if (do_check) {
        err = static_cast<double>(arma::accu(arma::abs(pi_next - pi_samp_f)));
      }
      pi_samp_f = std::move(pi_next);
    }

    if (do_check) {
      err_trace.push_back(err);
      outer_err_prev = err;
      if (err < tol) {
        ++it;
        break;
      }
    }
  }

  const arma::mat pi_samp = arma::conv_to<arma::mat>::from(pi_samp_f);
  const arma::mat pi_feat = arma::conv_to<arma::mat>::from(pi_feat_f);
  const auto costs = fused_unbalanced_cost_square_joint(
    Cx, Cy, Cx_sqr, Cy_sqr, M_samp, M_feat, wx, wy, wxy, pi_samp, pi_feat, rho_x, rho_y, epsilon
  );
  const int inner_total = std::accumulate(inner_iters_feat.begin(), inner_iters_feat.end(), 0) +
    std::accumulate(inner_iters_samp.begin(), inner_iters_samp.end(), 0);

  return Rcpp::List::create(
    Rcpp::Named("pi_samp") = pi_samp,
    Rcpp::Named("pi_feat") = pi_feat,
    Rcpp::Named("linear_cost") = costs.first,
    Rcpp::Named("fugw_cost") = costs.second,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = err,
    Rcpp::Named("err_trace") = err_trace,
    Rcpp::Named("inner_iters_feat") = Rcpp::IntegerVector(inner_iters_feat.begin(), inner_iters_feat.end()),
    Rcpp::Named("inner_iters_samp") = Rcpp::IntegerVector(inner_iters_samp.begin(), inner_iters_samp.end()),
    Rcpp::Named("inner_warm_feat") = Rcpp::LogicalVector(inner_warm_feat.begin(), inner_warm_feat.end()),
    Rcpp::Named("inner_warm_samp") = Rcpp::LogicalVector(inner_warm_samp.begin(), inner_warm_samp.end()),
    Rcpp::Named("inner_warm_fallback_feat") = Rcpp::LogicalVector(inner_fallback_feat.begin(), inner_fallback_feat.end()),
    Rcpp::Named("inner_warm_fallback_samp") = Rcpp::LogicalVector(inner_fallback_samp.begin(), inner_fallback_samp.end()),
    Rcpp::Named("inner_iters_total") = inner_total
  );
}

namespace {

struct FugwKlCoreResult {
  arma::mat pi_samp;
  arma::mat pi_feat;
  double linear_cost;
  double fugw_cost;
  int iterations;
  double error;
};

inline bool fugw_template_cache_compatible(
    const arma::mat& Cy_sqr_cache,
    const arma::mat& Cy_t_cache,
    const arma::mat& Cy) {
  return Cy_sqr_cache.n_rows == Cy.n_rows &&
    Cy_sqr_cache.n_cols == Cy.n_cols &&
    Cy_t_cache.n_rows == Cy.n_cols &&
    Cy_t_cache.n_cols == Cy.n_rows;
}

inline FugwKlCoreResult fugw_kl_square_core(
    const arma::mat& Cx,
    const arma::mat& Cy,
    const arma::vec& wx,
    const arma::vec& wy,
    const arma::vec& reg_marginals,
    double epsilon,
    double alpha,
    const arma::mat& M,
    const arma::mat& init_pi,
    int max_iter,
    double tol,
    int max_iter_ot,
    double tol_ot,
    bool rescale_plan,
    int check_every,
    const arma::mat* Cy_sqr_cache,
    const arma::mat* Cy_t_cache) {
  if (check_every <= 0) {
    check_every = 1;
  }

  const double rho_x = reg_marginals(0);
  const double rho_y = reg_marginals(1);
  const arma::mat M_samp = (alpha / 2.0) * M;
  const arma::mat M_feat = (alpha / 2.0) * M;
  const arma::mat Cx_sqr = Cx % Cx;
  arma::mat Cy_sqr_local;
  arma::mat Cy_t_local;
  const arma::mat* Cy_sqr_ptr = nullptr;
  const arma::mat* Cy_t_ptr = nullptr;
  if (Cy_sqr_cache != nullptr && Cy_t_cache != nullptr &&
      fugw_template_cache_compatible(*Cy_sqr_cache, *Cy_t_cache, Cy)) {
    Cy_sqr_ptr = Cy_sqr_cache;
    Cy_t_ptr = Cy_t_cache;
  } else {
    Cy_sqr_local = Cy % Cy;
    Cy_t_local = Cy.t();
    Cy_sqr_ptr = &Cy_sqr_local;
    Cy_t_ptr = &Cy_t_local;
  }
  const arma::mat wxy = wx * wy.t();
  const arma::vec log_wx = arma::log(wx + kTiny);
  const arma::vec log_wy = arma::log(wy + kTiny);

  arma::mat pi_samp = init_pi;
  arma::mat pi_feat = init_pi;
  SinkhornUnbalancedWorkspace ws_samp;
  SinkhornUnbalancedWorkspace ws_feat;
  arma::mat uot_cost(Cx.n_rows, Cy.n_rows);
  arma::mat scratch(Cx.n_rows, Cy.n_rows);
  arma::vec pi1(wx.n_elem);
  arma::vec pi2(wy.n_elem);
  arma::vec A(wx.n_elem);
  arma::vec B(wy.n_elem);

  double err = std::numeric_limits<double>::infinity();
  double outer_err_prev = std::numeric_limits<double>::infinity();
  int it = 0;
  for (; it < max_iter; ++it) {
    const bool do_check = ((it + 1) % check_every == 0);
    const double tol_ot_eff = fugw_effective_inner_tol(tol_ot, outer_err_prev, tol);

    {
      const double mass = arma::accu(pi_samp);
      uot_cost_matrix_kl_joint_inplace(
        Cx_sqr, *Cy_sqr_ptr, Cx, *Cy_t_ptr, M_feat, pi_samp, wx, wy, log_wx, log_wy, rho_x, rho_y, epsilon,
        uot_cost, scratch, pi1, pi2, A, B
      );
      SinkhornUnbalancedResult su = sinkhorn_unbalanced_kl(
        uot_cost,
        wx,
        wy,
        wxy,
        rho_x * mass,
        rho_y * mass,
        epsilon * mass,
        max_iter_ot,
        tol_ot_eff,
        pi_feat,
        ws_feat,
        fugw_enable_warm_start()
      );
      pi_feat = std::move(su.plan);
      if (rescale_plan) {
        const double mass_feat = arma::accu(pi_feat);
        if (mass_feat > 0.0) {
          pi_feat *= std::sqrt(mass / mass_feat);
        }
      }
    }

    {
      const double mass = arma::accu(pi_feat);
      uot_cost_matrix_kl_joint_inplace(
        Cx_sqr, *Cy_sqr_ptr, Cx, *Cy_t_ptr, M_samp, pi_feat, wx, wy, log_wx, log_wy, rho_x, rho_y, epsilon,
        uot_cost, scratch, pi1, pi2, A, B
      );
      SinkhornUnbalancedResult su = sinkhorn_unbalanced_kl(
        uot_cost,
        wx,
        wy,
        wxy,
        rho_x * mass,
        rho_y * mass,
        epsilon * mass,
        max_iter_ot,
        tol_ot_eff,
        pi_samp,
        ws_samp,
        fugw_enable_warm_start()
      );
      arma::mat pi_next = std::move(su.plan);
      if (rescale_plan) {
        const double mass_samp = arma::accu(pi_next);
        if (mass_samp > 0.0) {
          pi_next *= std::sqrt(mass / mass_samp);
        }
      }
      if (do_check) {
        err = arma::accu(arma::abs(pi_next - pi_samp));
      }
      pi_samp = std::move(pi_next);
    }

    if (do_check) {
      outer_err_prev = err;
      if (err < tol) {
        ++it;
        break;
      }
    }
  }

  const auto costs = fused_unbalanced_cost_square_joint(
    Cx, Cy, Cx_sqr, *Cy_sqr_ptr, M_samp, M_feat, wx, wy, wxy, pi_samp, pi_feat, rho_x, rho_y, epsilon
  );

  FugwKlCoreResult out;
  out.pi_samp = std::move(pi_samp);
  out.pi_feat = std::move(pi_feat);
  out.linear_cost = costs.first;
  out.fugw_cost = costs.second;
  out.iterations = it;
  out.error = err;
  return out;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_fugw_kl_square_batch(
    const Rcpp::List& Cx_list,
    const Rcpp::List& wx_list,
    const Rcpp::List& M_list,
    const arma::mat& Cy,
    const arma::vec& wy,
    const arma::vec& reg_marginals,
    double epsilon,
    double alpha,
    int max_iter,
    double tol,
    int max_iter_ot,
    double tol_ot,
    bool rescale_plan,
    int check_every,
    const Rcpp::List& init_pi_list,
    int n_threads) {
  const int n_jobs = static_cast<int>(Cx_list.size());
  if (n_jobs <= 0) {
    Rcpp::stop("`Cx_list` must be non-empty.");
  }
  if (wx_list.size() != Cx_list.size() || M_list.size() != Cx_list.size()) {
    Rcpp::stop("`Cx_list`, `wx_list`, and `M_list` must have the same length.");
  }
  if (init_pi_list.size() > 0 && init_pi_list.size() != Cx_list.size()) {
    Rcpp::stop("`init_pi_list` must be empty or the same length as `Cx_list`.");
  }
  if (Cy.n_rows != Cy.n_cols || wy.n_elem != Cy.n_rows) {
    Rcpp::stop("`Cy` must be square and compatible with `wy`.");
  }
  if (reg_marginals.n_elem != 2) {
    Rcpp::stop("`reg_marginals` must have length 2.");
  }
  if (check_every <= 0) {
    check_every = 1;
  }
  if (n_threads <= 0) {
    n_threads = 1;
  }

  std::vector<Rcpp::NumericMatrix> Cx_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericVector> wx_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> M_refs(static_cast<std::size_t>(n_jobs));
  std::vector<Rcpp::NumericMatrix> init_refs(static_cast<std::size_t>(n_jobs));
  std::vector<unsigned char> has_init(static_cast<std::size_t>(n_jobs), 0);
  for (int i = 0; i < n_jobs; ++i) {
    Cx_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(Cx_list[i]);
    wx_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericVector>(wx_list[i]);
    M_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(M_list[i]);
    if (init_pi_list.size() == Cx_list.size()) {
      SEXP init_i = init_pi_list[i];
      if (init_i != R_NilValue) {
        init_refs[static_cast<std::size_t>(i)] = Rcpp::as<Rcpp::NumericMatrix>(init_i);
        has_init[static_cast<std::size_t>(i)] = 1;
      }
    }
    const Rcpp::NumericMatrix& Cxi = Cx_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericVector& wxi = wx_refs[static_cast<std::size_t>(i)];
    const Rcpp::NumericMatrix& Mi = M_refs[static_cast<std::size_t>(i)];
    if (Cxi.nrow() != Cxi.ncol()) {
      Rcpp::stop("Every `Cx_list[[i]]` must be square.");
    }
    if (static_cast<int>(wxi.size()) != Cxi.nrow()) {
      Rcpp::stop("Each `wx_list[[i]]` must match `nrow(Cx_list[[i]])`.");
    }
    if (Mi.nrow() != Cxi.nrow() || Mi.ncol() != static_cast<int>(Cy.n_rows)) {
      Rcpp::stop("Each `M_list[[i]]` must have shape nrow(Cx_list[[i]]) x nrow(Cy).");
    }
  }

  const arma::mat Cy_sqr = Cy % Cy;
  const arma::mat Cy_t = Cy.t();

  std::vector<arma::mat> pi_samp(static_cast<std::size_t>(n_jobs));
  std::vector<double> fugw_cost(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<double> linear_cost(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<double> errs(static_cast<std::size_t>(n_jobs), NA_REAL);
  std::vector<int> iters(static_cast<std::size_t>(n_jobs), 0);
  std::vector<double> kernel_ms(static_cast<std::size_t>(n_jobs), 0.0);

#ifdef _OPENMP
  const int max_threads = omp_get_max_threads();
  const int used_threads = std::min(std::max(1, n_threads), max_threads);
  omp_set_schedule((n_jobs >= (2 * used_threads)) ? omp_sched_static : omp_sched_dynamic, 1);
#pragma omp parallel for schedule(runtime) num_threads(used_threads) if(n_jobs > 1 && used_threads > 1)
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& Cxref = Cx_refs[idx];
    const Rcpp::NumericVector& wxref = wx_refs[idx];
    const Rcpp::NumericMatrix& Mref = M_refs[idx];
    const arma::mat Cx(const_cast<double*>(Cxref.begin()), static_cast<arma::uword>(Cxref.nrow()), static_cast<arma::uword>(Cxref.ncol()), false, true);
    const arma::vec wx(const_cast<double*>(wxref.begin()), static_cast<arma::uword>(wxref.size()), false, true);
    const arma::mat M(const_cast<double*>(Mref.begin()), static_cast<arma::uword>(Mref.nrow()), static_cast<arma::uword>(Mref.ncol()), false, true);
    arma::mat init_i;
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    } else {
      init_i = wx * wy.t();
    }
    const auto t0 = std::chrono::steady_clock::now();
    FugwKlCoreResult out = fugw_kl_square_core(
      Cx, Cy, wx, wy, reg_marginals, epsilon, alpha, M, init_i,
      max_iter, tol, max_iter_ot, tol_ot, rescale_plan, check_every,
      &Cy_sqr, &Cy_t
    );
    const auto t1 = std::chrono::steady_clock::now();
    pi_samp[idx] = std::move(out.pi_samp);
    fugw_cost[idx] = out.fugw_cost;
    linear_cost[idx] = out.linear_cost;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
#else
  const int used_threads = 1;
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    const Rcpp::NumericMatrix& Cxref = Cx_refs[idx];
    const Rcpp::NumericVector& wxref = wx_refs[idx];
    const Rcpp::NumericMatrix& Mref = M_refs[idx];
    const arma::mat Cx(const_cast<double*>(Cxref.begin()), static_cast<arma::uword>(Cxref.nrow()), static_cast<arma::uword>(Cxref.ncol()), false, true);
    const arma::vec wx(const_cast<double*>(wxref.begin()), static_cast<arma::uword>(wxref.size()), false, true);
    const arma::mat M(const_cast<double*>(Mref.begin()), static_cast<arma::uword>(Mref.nrow()), static_cast<arma::uword>(Mref.ncol()), false, true);
    arma::mat init_i;
    if (has_init[idx]) {
      const Rcpp::NumericMatrix& init_ref = init_refs[idx];
      init_i = arma::mat(const_cast<double*>(init_ref.begin()), static_cast<arma::uword>(init_ref.nrow()), static_cast<arma::uword>(init_ref.ncol()), false, true);
    } else {
      init_i = wx * wy.t();
    }
    const auto t0 = std::chrono::steady_clock::now();
    FugwKlCoreResult out = fugw_kl_square_core(
      Cx, Cy, wx, wy, reg_marginals, epsilon, alpha, M, init_i,
      max_iter, tol, max_iter_ot, tol_ot, rescale_plan, check_every,
      &Cy_sqr, &Cy_t
    );
    const auto t1 = std::chrono::steady_clock::now();
    pi_samp[idx] = std::move(out.pi_samp);
    fugw_cost[idx] = out.fugw_cost;
    linear_cost[idx] = out.linear_cost;
    errs[idx] = out.error;
    iters[idx] = out.iterations;
    kernel_ms[idx] = std::chrono::duration<double, std::milli>(t1 - t0).count();
  }
#endif

  Rcpp::List out_plans(n_jobs);
  for (int i = 0; i < n_jobs; ++i) {
    out_plans[i] = pi_samp[static_cast<std::size_t>(i)];
  }

  return Rcpp::List::create(
    Rcpp::Named("pi_samp") = out_plans,
    Rcpp::Named("fugw_cost") = fugw_cost,
    Rcpp::Named("linear_cost") = linear_cost,
    Rcpp::Named("iterations") = iters,
    Rcpp::Named("error") = errs,
    Rcpp::Named("used_threads") = used_threads,
    Rcpp::Named("kernel_ms") = kernel_ms
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_sampled_gromov_wasserstein_entropic_square(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    int nb_p,
    int nb_q,
    double epsilon,
    int max_iter,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    bool symmetric,
    const arma::mat& init_plan,
    bool verbose,
    bool use_mixed_precision) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  if (C1.n_rows != p.n_elem || C2.n_rows != q.n_elem) {
    Rcpp::stop("Weight vectors must match matrix sizes.");
  }
  if (!(epsilon > 0.0)) {
    Rcpp::stop("`epsilon` must be positive for the entropic sampled GW kernel.");
  }
  if (nb_p < 1 || nb_q < 1) {
    Rcpp::stop("`nb_p` and `nb_q` must be >= 1.");
  }
  if (init_plan.n_rows != C1.n_rows || init_plan.n_cols != C2.n_rows) {
    Rcpp::stop("`init_plan` has incompatible shape.");
  }

  Rcpp::RNGScope scope;
  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const int nbp_eff = std::min<int>(nb_p, static_cast<int>(ns));
  const int nbq_eff = std::min<int>(nb_q, static_cast<int>(nt));
  const double exp_floor = std::exp(-200.0);

  arma::mat T = init_plan;
  arma::mat Lik(ns, nt, arma::fill::zeros);
  arma::mat Lik_eff(ns, nt, arma::fill::zeros);
  const arma::mat C1_sq = C1 % C1;
  const arma::mat C1_t = C1.t();
  const arma::mat C1_sq_t = C1_sq.t();
  const arma::mat C2_t = C2.t();
  arma::vec row_prob(nt, arma::fill::zeros);
  arma::rowvec mu2(nt, arma::fill::zeros);
  arma::rowvec mu2_sq(nt, arma::fill::zeros);
  arma::vec mu2_t(nt, arma::fill::zeros);
  arma::vec mu2_t_sq(nt, arma::fill::zeros);
  arma::vec u_ws;
  arma::vec v_ws;
  arma::fvec u_ws_f;
  arma::fvec v_ws_f;
  const arma::fvec p_f = arma::conv_to<arma::fvec>::from(p);
  const arma::fvec q_f = arma::conv_to<arma::fvec>::from(q);
  const float epsilon_f = static_cast<float>(epsilon);
  const float sinkhorn_tol_f = static_cast<float>(std::max(sinkhorn_tol, 1e-6));

  int continue_small = 0;
  int it_last = 0;
  double change = std::numeric_limits<double>::infinity();

  for (int it = 0; it < max_iter; ++it) {
    it_last = it + 1;
    if (((it + 1) % 5) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const arma::uvec idx0 = weighted_sample_no_replace(p, nbp_eff);
    Lik.zeros();

    for (arma::uword r = 0; r < idx0.n_elem; ++r) {
      const arma::uword i = idx0[r];

      row_prob = T.row(i).t();
      double srow = arma::accu(row_prob);
      if (!(srow > 0.0) || !std::isfinite(srow)) {
        row_prob = q;
      } else {
        row_prob /= srow;
      }

      const int nnz = static_cast<int>(arma::accu(row_prob > 0.0));
      const bool replace_q = (nnz < nbq_eff);
      const arma::uvec idx1 = replace_q
        ? weighted_sample_replace(row_prob, nbq_eff)
        : weighted_sample_no_replace(row_prob, nbq_eff);

      mu2.zeros();
      mu2_sq.zeros();
      for (arma::uword t = 0; t < idx1.n_elem; ++t) {
        const arma::uword j = idx1[t];
        const double* row_j = C2_t.colptr(j);
        double* mu_ptr = mu2.memptr();
        double* mu_sq_ptr = mu2_sq.memptr();
        for (arma::uword b = 0; b < nt; ++b) {
          const double v = row_j[b];
          mu_ptr[b] += v;
          mu_sq_ptr[b] += v * v;
        }
      }
      const double kq_inv = 1.0 / static_cast<double>(std::max<arma::uword>(1, idx1.n_elem));
      mu2 *= kq_inv;
      mu2_sq *= kq_inv;

      bool use_asym = false;
      if (!symmetric) {
        use_asym = (R::runif(0.0, 1.0) > 0.5);
      }

      const double* c1_ptr = nullptr;
      const double* c1_sq_ptr = nullptr;
      if (use_asym) {
        mu2_t.zeros();
        mu2_t_sq.zeros();
        for (arma::uword t = 0; t < idx1.n_elem; ++t) {
          const arma::uword j = idx1[t];
          const double* col_j = C2.colptr(j);
          double* mu_t_ptr = mu2_t.memptr();
          double* mu_t_sq_ptr = mu2_t_sq.memptr();
          for (arma::uword b = 0; b < nt; ++b) {
            const double v = col_j[b];
            mu_t_ptr[b] += v;
            mu_t_sq_ptr[b] += v * v;
          }
        }
        mu2_t *= kq_inv;
        mu2_t_sq *= kq_inv;
        c1_ptr = C1.colptr(i);
        c1_sq_ptr = C1_sq.colptr(i);
      } else {
        c1_ptr = C1_t.colptr(i);
        c1_sq_ptr = C1_sq_t.colptr(i);
      }

      for (arma::uword b = 0; b < nt; ++b) {
        const double mb = use_asym ? mu2_t[b] : mu2[b];
        const double mb_sq = use_asym ? mu2_t_sq[b] : mu2_sq[b];
        double* lk_col = Lik.colptr(b);
        for (arma::uword a = 0; a < ns; ++a) {
          lk_col[a] += c1_sq_ptr[a] + mb_sq - 2.0 * c1_ptr[a] * mb;
        }
      }
    }

    const double max_lik = Lik.max();
    if (std::isfinite(max_lik) && max_lik > 0.0) {
      Lik /= max_lik;
    }

    const double* Tp = T.memptr();
    const double* Lp = Lik.memptr();
    double* Lep = Lik_eff.memptr();
    for (arma::uword idx = 0; idx < T.n_elem; ++idx) {
      const double tij = Tp[idx];
      if (tij <= exp_floor || !std::isfinite(tij)) {
        Lep[idx] = std::numeric_limits<double>::infinity();
      } else {
        Lep[idx] = Lp[idx] - epsilon * std::log(tij);
      }
    }

    arma::mat new_T;
    if (use_mixed_precision) {
      const arma::fmat Lik_eff_f = arma::conv_to<arma::fmat>::from(Lik_eff);
      SinkhornBalancedResultF sk_f = sinkhorn_balanced_f(
        p_f, q_f, Lik_eff_f, epsilon_f, sinkhorn_max_iter, sinkhorn_tol_f, u_ws_f, v_ws_f
      );
      new_T = arma::conv_to<arma::mat>::from(sk_f.plan);
      u_ws_f = std::move(sk_f.u);
      v_ws_f = std::move(sk_f.v);
      bool finite_plan = true;
      const double* new_ptr = new_T.memptr();
      for (arma::uword idx = 0; idx < new_T.n_elem; ++idx) {
        if (!std::isfinite(new_ptr[idx])) {
          finite_plan = false;
          break;
        }
      }
      if (!finite_plan) {
        SinkhornBalancedResult sk = sinkhorn_balanced(
          p, q, Lik_eff, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
        );
        new_T = std::move(sk.plan);
        u_ws = std::move(sk.u);
        v_ws = std::move(sk.v);
      }
    } else {
      SinkhornBalancedResult sk = sinkhorn_balanced(
        p, q, Lik_eff, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
      );
      new_T = std::move(sk.plan);
      u_ws = std::move(sk.u);
      v_ws = std::move(sk.v);
    }

    const arma::mat diff = T - new_T;
    change = arma::accu(diff % diff) / static_cast<double>(diff.n_elem);
    if (!std::isfinite(change)) {
      break;
    }

    if (change <= 1e-19) {
      ++continue_small;
      if (continue_small > 100) {
        T = std::move(new_T);
        break;
      }
    } else {
      continue_small = 0;
    }

    if (verbose && (((it + 1) % 10) == 0 || it == 0)) {
      Rcpp::Rcout << "iter=" << (it + 1) << " change=" << change << "\n";
    }

    T = std::move(new_T);
  }

  return Rcpp::List::create(
    Rcpp::Named("plan") = T,
    Rcpp::Named("iterations") = it_last,
    Rcpp::Named("change") = change
  );
}

inline void coord_metric_row(
    const arma::mat& X,
    const arma::vec& x_norm2,
    arma::uword idx,
    bool metric_euclidean,
    arma::vec& xi_buf,
    arma::vec& dot_buf,
    arma::vec& d_metric,
    arma::vec& d_metric_sq) {
  const arma::uword n = X.n_rows;
  const arma::uword d = X.n_cols;
  if (xi_buf.n_elem != d) {
    xi_buf.set_size(d);
  }
  for (arma::uword c = 0; c < d; ++c) {
    xi_buf[c] = X.colptr(c)[idx];
  }

  if (dot_buf.n_elem != n) {
    dot_buf.set_size(n);
  }
  const bool small_dim_manual = (d <= 16);
  if (small_dim_manual) {
    double* dot_ptr = dot_buf.memptr();
    const double* xi_ptr = xi_buf.memptr();
    for (arma::uword a = 0; a < n; ++a) {
      double acc = 0.0;
      for (arma::uword c = 0; c < d; ++c) {
        acc += X.colptr(c)[a] * xi_ptr[c];
      }
      dot_ptr[a] = acc;
    }
  } else {
    dgemv_n(X, xi_buf, dot_buf);
  }

  const double self = x_norm2[idx];
  for (arma::uword a = 0; a < n; ++a) {
    double d2 = x_norm2[a] + self - 2.0 * dot_buf[a];
    if (!std::isfinite(d2) || d2 < 0.0) {
      d2 = 0.0;
    }
    if (metric_euclidean) {
      const double d = std::sqrt(d2);
      d_metric[a] = d;
      d_metric_sq[a] = d * d;
    } else {
      d_metric[a] = d2;
      d_metric_sq[a] = d2 * d2;
    }
  }
}

inline void coord_metric_rows_mean_batch(
    const arma::mat& X,
    const arma::vec& x_norm2,
    const arma::uvec& idx_rows,
    bool metric_euclidean,
    arma::mat& X_sel,
    arma::mat& dot_block,
    arma::vec& sel_norm2,
    arma::vec& mu,
    arma::vec& mu_sq) {
  const arma::uword n = X.n_rows;
  const arma::uword d = X.n_cols;
  const arma::uword k = idx_rows.n_elem;
  if (k == 0) {
    mu.zeros(n);
    mu_sq.zeros(n);
    return;
  }

  if (X_sel.n_rows != d || X_sel.n_cols < k) {
    X_sel.set_size(d, k);
  }
  if (dot_block.n_rows != n || dot_block.n_cols < k) {
    dot_block.set_size(n, k);
  }
  if (sel_norm2.n_elem < k) {
    sel_norm2.set_size(k);
  }
  if (mu.n_elem != n) {
    mu.set_size(n);
  }
  if (mu_sq.n_elem != n) {
    mu_sq.set_size(n);
  }
  mu.zeros();
  mu_sq.zeros();

  for (arma::uword t = 0; t < k; ++t) {
    const arma::uword j = idx_rows[t];
    sel_norm2[t] = x_norm2[j];
    for (arma::uword c = 0; c < d; ++c) {
      X_sel(c, t) = X.colptr(c)[j];
    }
  }

  arma::mat X_sel_view = X_sel.cols(0, k - 1);
  dgemm_nn(X, X_sel_view, dot_block);

  for (arma::uword t = 0; t < k; ++t) {
    const double self = sel_norm2[t];
    const double* dot_col = dot_block.colptr(t);
    for (arma::uword a = 0; a < n; ++a) {
      double d2 = x_norm2[a] + self - 2.0 * dot_col[a];
      if (!std::isfinite(d2) || d2 < 0.0) {
        d2 = 0.0;
      }
      if (metric_euclidean) {
        const double dv = std::sqrt(d2);
        mu[a] += dv;
        mu_sq[a] += dv * dv;
      } else {
        mu[a] += d2;
        mu_sq[a] += d2 * d2;
      }
    }
  }
  const double inv_k = 1.0 / static_cast<double>(k);
  mu *= inv_k;
  mu_sq *= inv_k;
}

// [[Rcpp::export]]
Rcpp::List cpp_sampled_gromov_wasserstein_coords_entropic_square(
    const arma::mat& X1,
    const arma::mat& X2,
    const arma::vec& p,
    const arma::vec& q,
    int nb_p,
    int nb_q,
    double epsilon,
    int max_iter,
    int sinkhorn_max_iter,
    double sinkhorn_tol,
    int metric_code,
    const arma::mat& init_plan,
    bool verbose,
    bool use_mixed_precision,
    bool deterministic_sampling) {
  if (X1.n_cols != X2.n_cols) {
    Rcpp::stop("`X1` and `X2` must have the same number of columns.");
  }
  if (X1.n_rows < 2 || X2.n_rows < 2) {
    Rcpp::stop("`X1` and `X2` must each have at least 2 rows.");
  }
  if (X1.n_rows != p.n_elem || X2.n_rows != q.n_elem) {
    Rcpp::stop("Weight vectors must match input sizes.");
  }
  if (!(epsilon > 0.0)) {
    Rcpp::stop("`epsilon` must be positive for the entropic sampled GW kernel.");
  }
  if (nb_p < 1 || nb_q < 1) {
    Rcpp::stop("`nb_p` and `nb_q` must be >= 1.");
  }
  if (metric_code != 0 && metric_code != 1) {
    Rcpp::stop("`metric_code` must be 0 (`euclidean`) or 1 (`sqeuclidean`).");
  }
  if (init_plan.n_rows != X1.n_rows || init_plan.n_cols != X2.n_rows) {
    Rcpp::stop("`init_plan` has incompatible shape.");
  }

  const bool metric_euclidean = (metric_code == 0);
  Rcpp::RNGScope scope;
  const arma::uword ns = X1.n_rows;
  const arma::uword nt = X2.n_rows;
  const int nbp_eff = std::min<int>(nb_p, static_cast<int>(ns));
  const int nbq_eff = std::min<int>(nb_q, static_cast<int>(nt));
  const double exp_floor = std::exp(-200.0);

  arma::vec x1_norm2 = arma::sum(X1 % X1, 1);
  arma::vec x2_norm2 = arma::sum(X2 % X2, 1);

  arma::vec xi1(X1.n_cols, arma::fill::zeros);
  arma::vec xi2(X2.n_cols, arma::fill::zeros);
  arma::vec dot1(ns, arma::fill::zeros);
  arma::vec dot2(nt, arma::fill::zeros);
  arma::vec d1(ns, arma::fill::zeros);
  arma::vec d1_sq(ns, arma::fill::zeros);
  arma::vec mu2(nt, arma::fill::zeros);
  arma::vec mu2_sq(nt, arma::fill::zeros);
  arma::mat X2_sel(X2.n_cols, static_cast<arma::uword>(std::max(1, nbq_eff)), arma::fill::zeros);
  arma::mat dot_block(nt, static_cast<arma::uword>(std::max(1, nbq_eff)), arma::fill::zeros);
  arma::vec sel_norm2(static_cast<arma::uword>(std::max(1, nbq_eff)), arma::fill::zeros);

  arma::mat T = init_plan;
  arma::mat Lik(ns, nt, arma::fill::zeros);
  arma::mat Lik_eff(ns, nt, arma::fill::zeros);
  arma::vec row_prob(nt, arma::fill::zeros);
  arma::vec u_ws;
  arma::vec v_ws;
  arma::fvec u_ws_f;
  arma::fvec v_ws_f;
  const arma::fvec p_f = arma::conv_to<arma::fvec>::from(p);
  const arma::fvec q_f = arma::conv_to<arma::fvec>::from(q);
  const float epsilon_f = static_cast<float>(epsilon);
  const float sinkhorn_tol_f = static_cast<float>(std::max(sinkhorn_tol, 1e-6));

  int continue_small = 0;
  int it_last = 0;
  double change = std::numeric_limits<double>::infinity();

  for (int it = 0; it < max_iter; ++it) {
    it_last = it + 1;
    if (((it + 1) % 5) == 0) {
      Rcpp::checkUserInterrupt();
    }

    const arma::uvec idx0 = deterministic_sampling
      ? deterministic_topk_indices(p, nbp_eff)
      : weighted_sample_no_replace(p, nbp_eff);
    Lik.zeros();

    for (arma::uword r = 0; r < idx0.n_elem; ++r) {
      const arma::uword i = idx0[r];

      row_prob = T.row(i).t();
      double srow = arma::accu(row_prob);
      if (!(srow > 0.0) || !std::isfinite(srow)) {
        row_prob = q;
      } else {
        row_prob /= srow;
      }

      arma::uvec idx1;
      if (deterministic_sampling) {
        idx1 = deterministic_topk_indices(row_prob, nbq_eff);
      } else {
        const int nnz = static_cast<int>(arma::accu(row_prob > 0.0));
        const bool replace_q = (nnz < nbq_eff);
        idx1 = replace_q
          ? weighted_sample_replace(row_prob, nbq_eff)
          : weighted_sample_no_replace(row_prob, nbq_eff);
      }

      coord_metric_row(X1, x1_norm2, i, metric_euclidean, xi1, dot1, d1, d1_sq);

      if (idx1.n_elem > 1) {
        coord_metric_rows_mean_batch(
          X2, x2_norm2, idx1, metric_euclidean,
          X2_sel, dot_block, sel_norm2, mu2, mu2_sq
        );
      } else {
        const arma::uword j = idx1[0];
        coord_metric_row(X2, x2_norm2, j, metric_euclidean, xi2, dot2, mu2, mu2_sq);
      }

      const double* d1_ptr = d1.memptr();
      const double* d1_sq_ptr = d1_sq.memptr();
      const double* mu_ptr = mu2.memptr();
      const double* mu_sq_ptr = mu2_sq.memptr();

      for (arma::uword b = 0; b < nt; ++b) {
        const double mb = mu_ptr[b];
        const double mb_sq = mu_sq_ptr[b];
        double* lk_col = Lik.colptr(b);
        for (arma::uword a = 0; a < ns; ++a) {
          lk_col[a] += d1_sq_ptr[a] + mb_sq - 2.0 * d1_ptr[a] * mb;
        }
      }
    }

    const double max_lik = Lik.max();
    if (std::isfinite(max_lik) && max_lik > 0.0) {
      Lik /= max_lik;
    }

    const double* Tp = T.memptr();
    const double* Lp = Lik.memptr();
    double* Lep = Lik_eff.memptr();
    for (arma::uword idx = 0; idx < T.n_elem; ++idx) {
      const double tij = Tp[idx];
      if (tij <= exp_floor || !std::isfinite(tij)) {
        Lep[idx] = std::numeric_limits<double>::infinity();
      } else {
        Lep[idx] = Lp[idx] - epsilon * std::log(tij);
      }
    }

    arma::mat new_T;
    if (use_mixed_precision) {
      const arma::fmat Lik_eff_f = arma::conv_to<arma::fmat>::from(Lik_eff);
      SinkhornBalancedResultF sk_f = sinkhorn_balanced_f(
        p_f, q_f, Lik_eff_f, epsilon_f, sinkhorn_max_iter, sinkhorn_tol_f, u_ws_f, v_ws_f
      );
      new_T = arma::conv_to<arma::mat>::from(sk_f.plan);
      u_ws_f = std::move(sk_f.u);
      v_ws_f = std::move(sk_f.v);
      bool finite_plan = true;
      const double* new_ptr = new_T.memptr();
      for (arma::uword idx = 0; idx < new_T.n_elem; ++idx) {
        if (!std::isfinite(new_ptr[idx])) {
          finite_plan = false;
          break;
        }
      }
      if (!finite_plan) {
        SinkhornBalancedResult sk = sinkhorn_balanced(
          p, q, Lik_eff, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
        );
        new_T = std::move(sk.plan);
        u_ws = std::move(sk.u);
        v_ws = std::move(sk.v);
      }
    } else {
      SinkhornBalancedResult sk = sinkhorn_balanced(
        p, q, Lik_eff, epsilon, sinkhorn_max_iter, sinkhorn_tol, u_ws, v_ws
      );
      new_T = std::move(sk.plan);
      u_ws = std::move(sk.u);
      v_ws = std::move(sk.v);
    }

    const arma::mat diff = T - new_T;
    change = arma::accu(diff % diff) / static_cast<double>(diff.n_elem);
    if (!std::isfinite(change)) {
      break;
    }

    if (change <= 1e-19) {
      ++continue_small;
      if (continue_small > 100) {
        T = std::move(new_T);
        break;
      }
    } else {
      continue_small = 0;
    }

    if (verbose && (((it + 1) % 10) == 0 || it == 0)) {
      Rcpp::Rcout << "iter=" << (it + 1) << " change=" << change << "\n";
    }

    T = std::move(new_T);
  }

  return Rcpp::List::create(
    Rcpp::Named("plan") = T,
    Rcpp::Named("iterations") = it_last,
    Rcpp::Named("change") = change
  );
}

// [[Rcpp::export]]
arma::mat cpp_srfgw_row_min_direction(const arma::mat& Mi, const arma::vec& p) {
  if (Mi.n_rows != p.n_elem) {
    Rcpp::stop("`Mi` rows must match length of `p`.");
  }
  const arma::uword ns = Mi.n_rows;
  const arma::uword nt = Mi.n_cols;
  arma::mat out(ns, nt, arma::fill::zeros);
  const double* mi_ptr = Mi.memptr();
  double* out_ptr = out.memptr();

  for (arma::uword i = 0; i < ns; ++i) {
    double row_min = std::numeric_limits<double>::infinity();
    for (arma::uword j = 0; j < nt; ++j) {
      const double v = mi_ptr[i + ns * j];
      if (v < row_min) {
        row_min = v;
      }
    }

    int cnt = 0;
    for (arma::uword j = 0; j < nt; ++j) {
      if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
        ++cnt;
      }
    }
    if (cnt <= 0) {
      cnt = 1;
    }
    const double val = p[i] / static_cast<double>(cnt);
    for (arma::uword j = 0; j < nt; ++j) {
      if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
        out_ptr[i + ns * j] = val;
      }
    }
  }
  return out;
}

struct GwSquareTermsResult {
  double loss;
  arma::mat grad;
};

inline GwSquareTermsResult gw_square_terms_exact(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    bool symmetric) {
  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const arma::vec pG = arma::sum(G, 1);
  const arma::vec qG = arma::sum(G, 0).t();

  const arma::mat C1_sq = C1 % C1;
  const arma::mat C2_sq = C2 % C2;
  const arma::vec left = C1_sq * pG;
  const arma::vec right = C2_sq.t() * qG;

  arma::mat constC(ns, nt, arma::fill::zeros);
  for (arma::uword j = 0; j < nt; ++j) {
    constC.col(j) = left + right[j];
  }

  arma::mat scratch;
  arma::mat Acur;
  dgemm_nn(C1, G, scratch);
  dgemm_nt(scratch, 2.0 * C2, Acur);
  arma::mat tens = constC - Acur;

  double loss = arma::accu(tens % G);
  arma::mat grad = 2.0 * tens;

  if (!symmetric) {
    const arma::mat C1t = C1.t();
    const arma::mat C2t = C2.t();
    const arma::mat C1t_sq = C1t % C1t;
    const arma::vec left_t = C1t_sq * pG;

    arma::mat constCt(ns, nt, arma::fill::zeros);
    for (arma::uword j = 0; j < nt; ++j) {
      constCt.col(j) = left_t + right[j];
    }

    arma::mat Acurt;
    dgemm_nn(C1t, G, scratch);
    dgemm_nn(scratch, 2.0 * C2, Acurt);
    arma::mat tenst = constCt - Acurt;

    loss = 0.5 * (loss + arma::accu(tenst % G));
    grad = 0.5 * (grad + 2.0 * tenst);
  }

  GwSquareTermsResult out;
  out.loss = loss;
  out.grad = std::move(grad);
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_gw_square_terms_square(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    bool symmetric = true) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  if (G.n_rows != C1.n_rows || G.n_cols != C2.n_rows) {
    Rcpp::stop("`G` has incompatible shape.");
  }

  GwSquareTermsResult out = gw_square_terms_exact(C1, C2, G, symmetric);
  return Rcpp::List::create(
    Rcpp::Named("loss") = out.loss,
    Rcpp::Named("grad") = out.grad
  );
}

struct EntropicPartialOtResult {
  arma::mat plan;
  int iterations;
  double error;
  std::vector<double> err_trace;
  bool finite;
};

inline EntropicPartialOtResult entropic_partial_wasserstein_core(
    const arma::vec& a,
    const arma::vec& b,
    const arma::mat& M,
    double reg,
    double m,
    int numItermax,
    double stopThr,
    bool verbose,
    bool log) {
  EntropicPartialOtResult out;
  out.iterations = 0;
  out.error = std::numeric_limits<double>::infinity();
  out.finite = true;

  const arma::uword ns = M.n_rows;
  const arma::uword nt = M.n_cols;
  const double floor_denom = 1e-300;
  const arma::vec dx(ns, arma::fill::ones);
  const arma::vec dy(nt, arma::fill::ones);

  arma::mat K = arma::exp(-M / reg);
  double Ksum = arma::accu(K);
  if (!std::isfinite(Ksum) || Ksum <= floor_denom) {
    Ksum = floor_denom;
  }
  K *= (m / Ksum);

  arma::mat q1(ns, nt, arma::fill::ones);
  arma::mat q2(ns, nt, arma::fill::ones);
  arma::mat q3(ns, nt, arma::fill::ones);

  arma::mat Kprev(ns, nt, arma::fill::zeros);
  arma::mat K1(ns, nt, arma::fill::zeros);
  arma::mat K2(ns, nt, arma::fill::zeros);
  arma::mat K_save(ns, nt, arma::fill::zeros);
  arma::mat denom(ns, nt, arma::fill::zeros);
  arma::vec row_sum(ns, arma::fill::zeros);
  arma::vec col_sum(nt, arma::fill::zeros);
  arma::vec row_scale(ns, arma::fill::zeros);
  arma::vec col_scale(nt, arma::fill::zeros);

  double err = std::numeric_limits<double>::infinity();
  int cpt = 0;
  std::vector<double> err_trace;
  if (log) {
    err_trace.reserve(static_cast<std::size_t>(std::max(2, numItermax / 10 + 2)));
  }

  while (err > stopThr && cpt < numItermax) {
    Kprev = K;

    K %= q1;
    row_sum = arma::sum(K, 1);
    for (arma::uword i = 0; i < ns; ++i) {
      double denom_row = row_sum[i];
      if (!std::isfinite(denom_row) || denom_row <= floor_denom) {
        denom_row = floor_denom;
      }
      const double ratio = a[i] / denom_row;
      row_scale[i] = (ratio < dx[i]) ? ratio : dx[i];
    }
    K1 = K;
    K1.each_col() %= row_scale;
    denom = K1;
    double* denom_ptr = denom.memptr();
    const arma::uword denom_n = denom.n_elem;
    for (arma::uword idx = 0; idx < denom_n; ++idx) {
      if (!std::isfinite(denom_ptr[idx]) || denom_ptr[idx] <= floor_denom) {
        denom_ptr[idx] = floor_denom;
      }
    }
    q1 %= (Kprev / denom);

    K_save = K1;
    K1 %= q2;
    col_sum = arma::sum(K1, 0).t();
    for (arma::uword j = 0; j < nt; ++j) {
      double denom_col = col_sum[j];
      if (!std::isfinite(denom_col) || denom_col <= floor_denom) {
        denom_col = floor_denom;
      }
      const double ratio = b[j] / denom_col;
      col_scale[j] = (ratio < dy[j]) ? ratio : dy[j];
    }
    K2 = K1;
    K2.each_row() %= col_scale.t();
    denom = K2;
    denom_ptr = denom.memptr();
    for (arma::uword idx = 0; idx < denom_n; ++idx) {
      if (!std::isfinite(denom_ptr[idx]) || denom_ptr[idx] <= floor_denom) {
        denom_ptr[idx] = floor_denom;
      }
    }
    q2 %= (K_save / denom);

    K_save = K2;
    K2 %= q3;
    Ksum = arma::accu(K2);
    if (!std::isfinite(Ksum) || Ksum <= floor_denom) {
      Ksum = floor_denom;
    }
    K = K2 * (m / Ksum);
    denom = K;
    denom_ptr = denom.memptr();
    for (arma::uword idx = 0; idx < denom_n; ++idx) {
      if (!std::isfinite(denom_ptr[idx]) || denom_ptr[idx] <= floor_denom) {
        denom_ptr[idx] = floor_denom;
      }
    }
    q3 %= (K_save / denom);

    bool finite_plan = true;
    const double* k_ptr = K.memptr();
    const arma::uword k_n = K.n_elem;
    for (arma::uword idx = 0; idx < k_n; ++idx) {
      if (!std::isfinite(k_ptr[idx])) {
        finite_plan = false;
        break;
      }
    }
    if (!finite_plan) {
      out.finite = false;
      K = Kprev;
      Rcpp::warning("Numerical instability in entropic partial OT at iteration %d.", cpt);
      break;
    }

    if ((cpt % 10) == 0) {
      arma::mat diff = Kprev - K;
      err = std::sqrt(arma::accu(diff % diff));
      if (log) {
        err_trace.push_back(err);
      }
      if (verbose && ((cpt % 200) == 0 || cpt == 0)) {
        Rcpp::Rcout << "it=" << cpt << " err=" << err << "\n";
      }
    }

    cpt += 1;
  }

  out.plan = std::move(K);
  out.iterations = cpt;
  out.error = err;
  out.err_trace = std::move(err_trace);
  return out;
}

// [[Rcpp::export]]
SEXP cpp_entropic_partial_wasserstein(
    const arma::vec& a,
    const arma::vec& b,
    const arma::mat& M,
    double reg,
    double m,
    int numItermax = 1000,
    double stopThr = 1e-100,
    bool verbose = false,
    bool log = false) {
  if (!std::isfinite(reg) || reg <= 0) {
    Rcpp::stop("`reg` must be positive.");
  }
  if (M.n_rows != a.n_elem || M.n_cols != b.n_elem) {
    Rcpp::stop("`a` and `b` must match the dimensions of `M`.");
  }
  if (!std::isfinite(m) || m < 0) {
    Rcpp::stop("`m` must be finite and >= 0.");
  }
  if (numItermax < 1) {
    Rcpp::stop("`numItermax` must be >= 1.");
  }

  const EntropicPartialOtResult core = entropic_partial_wasserstein_core(
    a, b, M, reg, m, numItermax, stopThr, verbose, log
  );
  if (log) {
    Rcpp::List log_obj = Rcpp::List::create(
      Rcpp::Named("err") = core.err_trace,
      Rcpp::Named("partial_w_dist") = arma::accu(M % core.plan)
    );
    return Rcpp::wrap(Rcpp::List::create(
      Rcpp::Named("plan") = core.plan,
      Rcpp::Named("log") = log_obj
    ));
  }
  return Rcpp::wrap(core.plan);
}

struct SrfgwSquareTermsResult {
  double quad;
  arma::mat grad;
};

inline arma::mat srfgw_row_min_direction_core(const arma::mat& Mi, const arma::vec& p) {
  if (Mi.n_rows != p.n_elem) {
    Rcpp::stop("`Mi` rows must match length of `p`.");
  }
  const arma::uword ns = Mi.n_rows;
  const arma::uword nt = Mi.n_cols;
  arma::mat out(ns, nt, arma::fill::zeros);
  const double* mi_ptr = Mi.memptr();
  double* out_ptr = out.memptr();

  for (arma::uword i = 0; i < ns; ++i) {
    double row_min = std::numeric_limits<double>::infinity();
    for (arma::uword j = 0; j < nt; ++j) {
      const double v = mi_ptr[i + ns * j];
      if (v < row_min) {
        row_min = v;
      }
    }

    int cnt = 0;
    for (arma::uword j = 0; j < nt; ++j) {
      if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
        ++cnt;
      }
    }
    if (cnt <= 0) {
      cnt = 1;
    }
    const double val = p[i] / static_cast<double>(cnt);
    for (arma::uword j = 0; j < nt; ++j) {
      if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
        out_ptr[i + ns * j] = val;
      }
    }
  }
  return out;
}

inline SrfgwSquareTermsResult srfgw_square_terms_exact(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    const arma::vec& p,
    bool symmetric) {
  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const arma::vec qG = arma::sum(G, 0).t();
  const arma::vec ones_p = arma::ones<arma::vec>(ns);

  const arma::mat C1_sq = C1 % C1;
  const arma::mat C2_sq = C2 % C2;
  const arma::mat fC2t = C2_sq.t();

  arma::mat constC(ns, nt, arma::fill::zeros);
  const arma::vec left = C1_sq * p;
  const arma::vec right = C2_sq * qG;
  for (arma::uword j = 0; j < nt; ++j) {
    constC.col(j) = left + right[j];
  }

  arma::mat scratch;
  arma::mat Acur;
  dgemm_nn(C1, G, scratch);
  dgemm_nt(scratch, 2.0 * C2, Acur);
  arma::mat tens = constC - Acur;

  double quad = arma::accu(tens % G);
  arma::mat grad = 2.0 * tens;

  if (!symmetric) {
    const arma::mat C1t = C1.t();
    const arma::mat C2t = C2.t();
    const arma::mat C1t_sq = C1t % C1t;
    const arma::mat C2t_sq = C2t % C2t;

    arma::mat constCt(ns, nt, arma::fill::zeros);
    const arma::vec left_t = C1t_sq * p;
    const arma::vec right_t = C2t_sq * qG;
    for (arma::uword j = 0; j < nt; ++j) {
      constCt.col(j) = left_t + right_t[j];
    }

    arma::mat Acurt;
    dgemm_nn(C1t, G, scratch);
    dgemm_nt(scratch, 2.0 * C2t, Acurt);
    arma::mat tenst = constCt - Acurt;

    quad = 0.5 * (quad + arma::accu(tenst % G));
    grad = 0.5 * (grad + 2.0 * tenst);
  }

  SrfgwSquareTermsResult out;
  out.quad = quad;
  out.grad = std::move(grad);
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_semirelaxed_fgw_exact_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    double alpha,
    bool symmetric,
    const arma::mat& init_plan,
    int max_iter,
    double tol_rel,
    double tol_abs) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  const bool skip_M = (alpha >= 1.0) && (M.n_elem == 0);
  if (!skip_M && (M.n_rows != C1.n_rows || M.n_cols != C2.n_rows)) {
    Rcpp::stop("`M` must have shape nrow(C1) x nrow(C2).");
  }
  if (p.n_elem != C1.n_rows) {
    Rcpp::stop("`p` must match `nrow(C1)`.");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }

  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const double lin_w = 1.0 - alpha;
  const bool use_lin = !skip_M && lin_w != 0.0;

  arma::mat G;
  if (init_plan.n_elem == 0) {
    G = p * (arma::ones<arma::vec>(nt) / static_cast<double>(nt)).t();
  } else {
    if (init_plan.n_rows != ns || init_plan.n_cols != nt) {
      Rcpp::stop("`init_plan` has incompatible shape.");
    }
    G = init_plan;
  }

  const arma::mat hC2 = 2.0 * C2;
  const arma::mat C1_sq = C1 % C1;
  const arma::mat C2_sq = C2 % C2;
  const arma::mat fC2t = C2_sq.t();
  const arma::vec left = C1_sq * p;
  arma::mat C1t;
  arma::mat hC2t;
  arma::vec left_t;
  arma::mat Acurt;
  if (!symmetric) {
    C1t = C1.t();
    hC2t = 2.0 * C2.t();
    left_t = (C1t % C1t) * p;
  }

  arma::vec qG = arma::sum(G, 0).t();
  arma::mat scratch;
  arma::mat Acur;
  dgemm_nn(C1, G, scratch);
  dgemm_nt(scratch, hC2, Acur);
  arma::mat tens(ns, nt);
  {
    const arma::vec right = C2_sq * qG;
    for (arma::uword j = 0; j < nt; ++j) {
      tens.col(j) = left + right[j] - Acur.col(j);
    }
  }
  arma::mat grad = 2.0 * tens;
  double quad_raw = arma::accu(tens % G);
  if (!symmetric) {
    dgemm_nn(C1t, G, scratch);
    dgemm_nt(scratch, hC2t, Acurt);
    arma::mat tenst(ns, nt);
    const arma::vec right_t = fC2t * qG;
    for (arma::uword j = 0; j < nt; ++j) {
      tenst.col(j) = left_t + right_t[j] - Acurt.col(j);
    }
    quad_raw = 0.5 * (quad_raw + arma::accu(tenst % G));
    grad = 0.5 * (grad + 2.0 * tenst);
  }

  double lin_loss = use_lin ? lin_w * arma::accu(M % G) : 0.0;
  double cost = lin_loss + alpha * quad_raw;

  std::vector<double> loss_trace;
  loss_trace.reserve(static_cast<std::size_t>(max_iter) + 1u);
  loss_trace.push_back(cost);

  double rel_delta = std::numeric_limits<double>::infinity();
  double abs_delta = std::numeric_limits<double>::infinity();
  int it = 0;

  for (int k = 0; k < max_iter; ++k) {
    arma::mat Mi = alpha * grad;
    if (use_lin) {
      Mi += lin_w * M;
    }
    const arma::mat Gc = srfgw_row_min_direction_core(Mi, p);
    const arma::mat delta = Gc - G;
    const arma::vec qdelta = arma::sum(delta, 0).t();

    arma::mat dot;
    dgemm_nn(C1, delta, scratch);
    dgemm_nt(scratch, hC2, dot);

    const arma::vec right_qG = C2_sq * qG;
    const arma::vec right_qdelta = C2_sq * qdelta;
    const double a_ls = alpha * (arma::dot(right_qdelta, qdelta) - arma::accu(dot % delta));
    double b_ls = alpha * (
      arma::dot(right_qdelta, qG) - arma::accu(dot % G) +
        arma::dot(right_qG, qdelta) - arma::accu(Acur % delta)
    );
    if (use_lin) {
      b_ls += lin_w * arma::accu(M % delta);
    }

    double step = solve_1d_linesearch_quad(a_ls, b_ls);
    step = std::min(1.0, std::max(0.0, step));
    const double new_cost = cost + a_ls * step * step + b_ls * step;

    G += step * delta;
    qG += step * qdelta;
    Acur += step * dot;
    {
      const arma::vec right = C2_sq * qG;
      for (arma::uword j = 0; j < nt; ++j) {
        tens.col(j) = left + right[j] - Acur.col(j);
      }
    }
    grad = 2.0 * tens;
    quad_raw = arma::accu(tens % G);
    if (!symmetric) {
      arma::mat dott;
      dgemm_nn(C1t, delta, scratch);
      dgemm_nt(scratch, hC2t, dott);
      Acurt += step * dott;
      arma::mat tenst(ns, nt);
      const arma::vec right_t = fC2t * qG;
      for (arma::uword j = 0; j < nt; ++j) {
        tenst.col(j) = left_t + right_t[j] - Acurt.col(j);
      }
      quad_raw = 0.5 * (quad_raw + arma::accu(tenst % G));
      grad = 0.5 * (grad + 2.0 * tenst);
    }

    abs_delta = std::abs(new_cost - cost);
    rel_delta = abs_delta / (std::abs(new_cost) + 1e-15);
    cost = new_cost;
    loss_trace.push_back(cost);
    it = k + 1;

    if ((std::isfinite(rel_delta) && rel_delta <= tol_rel) || (std::isfinite(abs_delta) && abs_delta <= tol_abs)) {
      break;
    }
  }

  const arma::vec q = qG;
  lin_loss = use_lin ? lin_w * arma::accu(M % G) : 0.0;
  const double quad_loss = alpha * quad_raw;

  return Rcpp::List::create(
    Rcpp::Named("plan") = G,
    Rcpp::Named("q") = q,
    Rcpp::Named("lin_loss") = lin_loss,
    Rcpp::Named("quad_loss") = quad_loss,
    Rcpp::Named("srfgw_dist") = lin_loss + quad_loss,
    Rcpp::Named("srgw_dist") = quad_raw,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = rel_delta,
    Rcpp::Named("abs_error") = abs_delta,
    Rcpp::Named("loss_trace") = loss_trace,
    Rcpp::Named("symmetric") = symmetric
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_semirelaxed_fgw_cg_square_fast(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    double alpha,
    const arma::mat& init_plan,
    int max_iter,
    double tol_rel,
    double tol_abs,
    bool verbose,
    bool use_mixed_precision) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  const bool skip_M = (alpha >= 1.0) && (M.n_elem == 0);
  if (!skip_M && (M.n_rows != C1.n_rows || M.n_cols != C2.n_rows)) {
    Rcpp::stop("`M` must have shape nrow(C1) x nrow(C2).");
  }
  if (p.n_elem != C1.n_rows) {
    Rcpp::stop("`p` must match nrow(C1).");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }

  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const bool use_lin = !skip_M && (alpha < 1.0);
  if (use_mixed_precision) {
    const float alpha_f = static_cast<float>(alpha);
    const arma::fmat C1_f = arma::conv_to<arma::fmat>::from(C1);
    const arma::fmat C2_f = arma::conv_to<arma::fmat>::from(C2);
    const arma::fvec p_f = arma::conv_to<arma::fvec>::from(p);
    arma::fmat M_lin_f;
    if (use_lin) {
      M_lin_f = (1.0f - alpha_f) * arma::conv_to<arma::fmat>::from(M);
    }

    arma::fmat G_f;
    if (init_plan.n_elem == 0) {
      G_f = p_f * (arma::ones<arma::fvec>(nt) / static_cast<float>(nt)).t();
    } else {
      if (init_plan.n_rows != ns || init_plan.n_cols != nt) {
        Rcpp::stop("`init_plan` has incompatible shape.");
      }
      G_f = arma::conv_to<arma::fmat>::from(init_plan);
    }

    const arma::fmat hC1_f = C1_f;
    const arma::fmat hC2_f = 2.0f * C2_f;
    const arma::fmat C1_sq_f = C1_f % C1_f;
    const arma::fmat C2_sq_f = C2_f % C2_f;
    const arma::fmat fC2t_f = C2_sq_f.t();
    const arma::fvec left_f = C1_sq_f * p_f;

    arma::fvec q_f = arma::sum(G_f, 0).t();
    arma::fvec right_q_f = fC2t_f * q_f;

    arma::fmat Acur_f;
    arma::fmat scratch_f;
    sgemm_nn(hC1_f, G_f, scratch_f);
    sgemm_nt(scratch_f, hC2_f, Acur_f);

    arma::fmat tens_f(ns, nt, arma::fill::zeros);
    for (arma::uword j = 0; j < nt; ++j) {
      tens_f.col(j) = left_f + right_q_f[j] - Acur_f.col(j);
    }
    arma::fmat grad_f = 2.0f * tens_f;

    double quad_raw = static_cast<double>(arma::accu(tens_f % G_f));
    double lin_loss = use_lin ? static_cast<double>(arma::accu(M_lin_f % G_f)) : 0.0;
    double cost = lin_loss + alpha * quad_raw;

    std::vector<double> loss_trace;
    loss_trace.reserve(static_cast<std::size_t>(max_iter) + 1u);
    loss_trace.push_back(cost);

    arma::fmat Mi_f(ns, nt, arma::fill::zeros);
    arma::fmat delta_f(ns, nt, arma::fill::zeros);
    arma::fvec qdelta_f(nt, arma::fill::zeros);
    arma::fvec right_delta_f(nt, arma::fill::zeros);
    arma::fmat dot_f(ns, nt, arma::fill::zeros);

    double rel_delta = std::numeric_limits<double>::infinity();
    double abs_delta = std::numeric_limits<double>::infinity();
    const double tol_rel_eff = std::max(tol_rel, 1e-6);
    const double tol_abs_eff = std::max(tol_abs, 1e-8);
    int it = 0;

    for (int k = 0; k < max_iter; ++k) {
      if (use_lin) {
        Mi_f = M_lin_f;
        Mi_f += alpha_f * grad_f;
      } else {
        Mi_f = alpha_f * grad_f;
      }

      delta_f = -G_f;
      qdelta_f = -q_f;
      const float* mi_ptr = Mi_f.memptr();
      float* d_ptr = delta_f.memptr();
      float* qd_ptr = qdelta_f.memptr();

      for (arma::uword i = 0; i < ns; ++i) {
        float row_min = std::numeric_limits<float>::infinity();
        for (arma::uword j = 0; j < nt; ++j) {
          const float v = mi_ptr[i + ns * j];
          if (v < row_min) {
            row_min = v;
          }
        }

        int cnt = 0;
        const float tie_tol = 1e-8f * (1.0f + std::fabs(row_min));
        for (arma::uword j = 0; j < nt; ++j) {
          if (mi_ptr[i + ns * j] <= row_min + tie_tol) {
            ++cnt;
          }
        }
        if (cnt <= 0) {
          cnt = 1;
        }
        const float val = p_f[i] / static_cast<float>(cnt);
        for (arma::uword j = 0; j < nt; ++j) {
          if (mi_ptr[i + ns * j] <= row_min + tie_tol) {
            d_ptr[i + ns * j] += val;
            qd_ptr[j] += val;
          }
        }
      }

      sgemm_nn(hC1_f, delta_f, scratch_f);
      sgemm_nt(scratch_f, hC2_f, dot_f);
      right_delta_f = fC2t_f * qdelta_f;

      const double sum_dot_delta = static_cast<double>(arma::accu(dot_f % delta_f));
      const double sum_dot_G = static_cast<double>(arma::accu(dot_f % G_f));
      const double sum_Acur_delta = static_cast<double>(arma::accu(Acur_f % delta_f));

      const double a_ls = alpha * (static_cast<double>(arma::dot(right_delta_f, qdelta_f)) - sum_dot_delta);
      double b_ls = alpha * (
        static_cast<double>(arma::dot(right_delta_f, q_f)) - sum_dot_G +
          static_cast<double>(arma::dot(right_q_f, qdelta_f)) - sum_Acur_delta
      );
      if (use_lin) {
        b_ls += static_cast<double>(arma::accu(M_lin_f % delta_f));
      }

      double step = solve_1d_linesearch_quad(a_ls, b_ls);
      step = std::min(1.0, std::max(0.0, step));
      const float step_f = static_cast<float>(step);
      const double new_cost = cost + a_ls * step * step + b_ls * step;

      G_f += step_f * delta_f;
      Acur_f += step_f * dot_f;
      q_f += step_f * qdelta_f;
      right_q_f += step_f * right_delta_f;

      for (arma::uword j = 0; j < nt; ++j) {
        tens_f.col(j) = left_f + right_q_f[j] - Acur_f.col(j);
      }
      grad_f = 2.0f * tens_f;

      abs_delta = std::abs(new_cost - cost);
      rel_delta = abs_delta / (std::abs(new_cost) + 1e-15);
      cost = new_cost;
      loss_trace.push_back(cost);
      it = k + 1;

      if (verbose && ((k + 1) % 25 == 0 || k == 0)) {
        Rcpp::Rcout << "iter=" << (k + 1) << " cost=" << cost << " rel=" << rel_delta << " abs=" << abs_delta << "\n";
      }

      if ((std::isfinite(rel_delta) && rel_delta <= tol_rel_eff) || (std::isfinite(abs_delta) && abs_delta <= tol_abs_eff)) {
        break;
      }
    }

    const arma::mat G = arma::conv_to<arma::mat>::from(G_f);
    const arma::vec q = arma::conv_to<arma::vec>::from(q_f);
    const SrfgwSquareTermsResult state_d = srfgw_square_terms_exact(C1, C2, G, p, true);
    quad_raw = state_d.quad;
    lin_loss = use_lin ? (1.0 - alpha) * arma::accu(M % G) : 0.0;
    const double quad_loss = alpha * quad_raw;
    const double srfgw_dist = lin_loss + quad_loss;

    return Rcpp::List::create(
      Rcpp::Named("plan") = G,
      Rcpp::Named("q") = q,
      Rcpp::Named("lin_loss") = lin_loss,
      Rcpp::Named("quad_loss") = quad_loss,
      Rcpp::Named("srfgw_dist") = srfgw_dist,
      Rcpp::Named("srgw_dist") = quad_raw,
      Rcpp::Named("iterations") = it,
      Rcpp::Named("error") = rel_delta,
      Rcpp::Named("abs_error") = abs_delta,
      Rcpp::Named("loss_trace") = loss_trace,
      Rcpp::Named("symmetric") = true
    );
  }

  const double lin_w = 1.0 - alpha;
  arma::mat M_lin;
  if (use_lin) {
    M_lin = lin_w * M;
  }

  arma::mat G;
  if (init_plan.n_elem == 0) {
    G = p * (arma::ones<arma::vec>(nt) / static_cast<double>(nt)).t();
  } else {
    if (init_plan.n_rows != ns || init_plan.n_cols != nt) {
      Rcpp::stop("`init_plan` has incompatible shape.");
    }
    G = init_plan;
  }

  const arma::mat hC1 = C1;
  const arma::mat hC2 = 2.0 * C2;
  const arma::mat C1_sq = C1 % C1;
  const arma::mat C2_sq = C2 % C2;
  const arma::mat fC2t = C2_sq.t();
  const arma::vec left = C1_sq * p;

  arma::vec q = arma::sum(G, 0).t();
  arma::vec right_q = fC2t * q;

  arma::mat Acur;
  arma::mat scratch;
  dgemm_nn(hC1, G, scratch);
  dgemm_nt(scratch, hC2, Acur);

  arma::mat tens(ns, nt, arma::fill::zeros);
  for (arma::uword j = 0; j < nt; ++j) {
    tens.col(j) = left + right_q[j] - Acur.col(j);
  }
  arma::mat grad = 2.0 * tens;

  double quad_raw = arma::accu(tens % G);
  double lin_loss = use_lin ? arma::accu(M_lin % G) : 0.0;
  double cost = lin_loss + alpha * quad_raw;

  std::vector<double> loss_trace;
  loss_trace.reserve(static_cast<std::size_t>(max_iter) + 1u);
  loss_trace.push_back(cost);

  arma::mat Mi(ns, nt, arma::fill::zeros);
  arma::mat delta(ns, nt, arma::fill::zeros);
  arma::vec qdelta(nt, arma::fill::zeros);
  arma::vec right_delta(nt, arma::fill::zeros);
  arma::mat dot(ns, nt, arma::fill::zeros);

  double rel_delta = std::numeric_limits<double>::infinity();
  double abs_delta = std::numeric_limits<double>::infinity();
  int it = 0;

  for (int k = 0; k < max_iter; ++k) {
    if (use_lin) {
      Mi = M_lin;
      Mi += alpha * grad;
    } else {
      Mi = alpha * grad;
    }

    delta = -G;
    qdelta = -q;
    const double* mi_ptr = Mi.memptr();
    double* d_ptr = delta.memptr();
    double* qd_ptr = qdelta.memptr();

    for (arma::uword i = 0; i < ns; ++i) {
      double row_min = std::numeric_limits<double>::infinity();
      for (arma::uword j = 0; j < nt; ++j) {
        const double v = mi_ptr[i + ns * j];
        if (v < row_min) {
          row_min = v;
        }
      }

      int cnt = 0;
      for (arma::uword j = 0; j < nt; ++j) {
        if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
          ++cnt;
        }
      }
      if (cnt <= 0) {
        cnt = 1;
      }
      const double val = p[i] / static_cast<double>(cnt);
      for (arma::uword j = 0; j < nt; ++j) {
        if (mi_ptr[i + ns * j] <= row_min + 1e-15) {
          d_ptr[i + ns * j] += val;
          qd_ptr[j] += val;
        }
      }
    }

    dgemm_nn(hC1, delta, scratch);
    dgemm_nt(scratch, hC2, dot);
    right_delta = fC2t * qdelta;

    const double sum_dot_delta = arma::accu(dot % delta);
    const double sum_dot_G = arma::accu(dot % G);
    const double sum_Acur_delta = arma::accu(Acur % delta);

    const double a_ls = alpha * (arma::dot(right_delta, qdelta) - sum_dot_delta);
    double b_ls = alpha * (
      arma::dot(right_delta, q) - sum_dot_G +
        arma::dot(right_q, qdelta) - sum_Acur_delta
    );
    if (use_lin) {
      b_ls += arma::accu(M_lin % delta);
    }

    double step = solve_1d_linesearch_quad(a_ls, b_ls);
    step = std::min(1.0, std::max(0.0, step));
    const double new_cost = cost + a_ls * step * step + b_ls * step;

    G += step * delta;
    Acur += step * dot;
    q += step * qdelta;
    right_q += step * right_delta;

    for (arma::uword j = 0; j < nt; ++j) {
      tens.col(j) = left + right_q[j] - Acur.col(j);
    }
    grad = 2.0 * tens;
    quad_raw = arma::accu(tens % G);
    lin_loss = use_lin ? arma::accu(M_lin % G) : 0.0;

    abs_delta = std::abs(new_cost - cost);
    rel_delta = abs_delta / (std::abs(new_cost) + 1e-15);
    cost = new_cost;
    loss_trace.push_back(cost);
    it = k + 1;

    if (verbose && ((k + 1) % 25 == 0 || k == 0)) {
      Rcpp::Rcout << "iter=" << (k + 1) << " cost=" << cost << " rel=" << rel_delta << " abs=" << abs_delta << "\n";
    }

    if ((std::isfinite(rel_delta) && rel_delta <= tol_rel) || (std::isfinite(abs_delta) && abs_delta <= tol_abs)) {
      break;
    }
  }

  const double quad_loss = alpha * quad_raw;
  const double srfgw_dist = lin_loss + quad_loss;

  return Rcpp::List::create(
    Rcpp::Named("plan") = G,
    Rcpp::Named("q") = q,
    Rcpp::Named("lin_loss") = lin_loss,
    Rcpp::Named("quad_loss") = quad_loss,
    Rcpp::Named("srfgw_dist") = srfgw_dist,
    Rcpp::Named("srgw_dist") = quad_raw,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = rel_delta,
    Rcpp::Named("abs_error") = abs_delta,
    Rcpp::Named("loss_trace") = loss_trace,
    Rcpp::Named("symmetric") = true
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_ot_sinkhorn(
    const arma::mat& M,
    const arma::vec& p,
    const arma::vec& q,
    double epsilon,
    int max_iter,
    double tol,
    bool use_log) {
  arma::vec u;
  arma::vec v;
  const SinkhornBalancedResult res = use_log
    ? sinkhorn_balanced_log(p, q, M, epsilon, max_iter, tol, u, v)
    : sinkhorn_balanced(p, q, M, epsilon, max_iter, tol, u, v);
  const double ot_dist = arma::accu(M % res.plan);
  return Rcpp::List::create(
    Rcpp::Named("plan") = res.plan,
    Rcpp::Named("ot_dist") = ot_dist,
    Rcpp::Named("iterations") = res.iters,
    Rcpp::Named("error") = res.err
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_ot_emd(
    const arma::mat& M,
    const arma::vec& p,
    const arma::vec& q,
    int max_iter,
    double tol) {
  const TransportSimplexResult res = transport_simplex_solve(M, p, q, max_iter, tol);
  const double ot_dist = arma::accu(M % res.plan);
  return Rcpp::List::create(
    Rcpp::Named("plan") = res.plan,
    Rcpp::Named("ot_dist") = ot_dist,
    Rcpp::Named("iterations") = res.iterations,
    Rcpp::Named("error") = res.converged ? 0.0 : R_PosInf,
    Rcpp::Named("lp_ok") = res.converged
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_ot_sinkhorn_unbalanced(
    const arma::mat& M,
    const arma::vec& a,
    const arma::vec& b,
    double epsilon,
    double rho1,
    double rho2,
    int max_iter,
    double tol,
    const arma::mat& init_plan) {
  SinkhornUnbalancedWorkspace ws;
  const arma::mat c = a * b.t();
  const SinkhornUnbalancedResult res = sinkhorn_unbalanced_kl(
    M, a, b, c, rho1, rho2, epsilon, max_iter, tol, init_plan, ws, true
  );
  const double ot_dist = arma::accu(M % res.plan);
  return Rcpp::List::create(
    Rcpp::Named("plan") = res.plan,
    Rcpp::Named("ot_dist") = ot_dist,
    Rcpp::Named("iterations") = res.iters,
    Rcpp::Named("error") = res.err,
    Rcpp::Named("mass") = arma::accu(res.plan)
  );
}

inline double partial_penalty_cost(const arma::mat& cost) {
  const double ma = cost.max();
  const double aa = arma::abs(cost).max();
  const double aa_use = (std::isfinite(aa) && aa > 0.0) ? aa : 1.0;
  const double ma_use = std::isfinite(ma) ? ma : 0.0;
  return ma_use + aa_use + 1.0;
}

inline arma::mat partial_transport_direction(
    const arma::mat& cost,
    const arma::vec& a,
    const arma::vec& b,
    double m,
    int nb_dummies,
    int lp_max_iter,
    double lp_tol) {
  if (nb_dummies < 1) {
    Rcpp::stop("`nb_dummies` must be >= 1.");
  }
  const arma::uword ns = a.n_elem;
  const arma::uword nt = b.n_elem;
  if (cost.n_rows != ns || cost.n_cols != nt) {
    Rcpp::stop("`cost` has incompatible shape for partial LP direction step.");
  }
  const double sum_a = arma::accu(a);
  const double sum_b = arma::accu(b);
  const double dummy_a = (sum_b - m) / static_cast<double>(nb_dummies);
  const double dummy_b = (sum_a - m) / static_cast<double>(nb_dummies);
  if (dummy_a <= 1e-15 && dummy_b <= 1e-15) {
    return transport_simplex_solve(cost, a, b, lp_max_iter, lp_tol).plan;
  }

  const arma::uword ns_ext = ns + static_cast<arma::uword>(nb_dummies);
  const arma::uword nt_ext = nt + static_cast<arma::uword>(nb_dummies);
  arma::vec a_ext(ns_ext, arma::fill::zeros);
  arma::vec b_ext(nt_ext, arma::fill::zeros);
  a_ext.head(ns) = a;
  b_ext.head(nt) = b;
  a_ext.tail(static_cast<arma::uword>(nb_dummies)).fill(std::max(dummy_a, 0.0));
  b_ext.tail(static_cast<arma::uword>(nb_dummies)).fill(std::max(dummy_b, 0.0));

  arma::mat cost_ext(ns_ext, nt_ext, arma::fill::zeros);
  cost_ext.submat(0, 0, ns - 1, nt - 1) = cost;
  const double penalty = partial_penalty_cost(cost);
  cost_ext.submat(ns, nt, ns_ext - 1, nt_ext - 1).fill(penalty);

  const TransportSimplexResult res = transport_simplex_solve(
    cost_ext, a_ext, b_ext, lp_max_iter, lp_tol
  );
  return res.plan.submat(0, 0, ns - 1, nt - 1);
}

// [[Rcpp::export]]
Rcpp::List cpp_partial_fgw_exact_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double m,
    double alpha,
    bool symmetric,
    const arma::mat& init_plan,
    int max_iter,
    double tol,
    int nb_dummies,
    int lp_max_iter,
    double lp_tol) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  const bool skip_M = (alpha >= 1.0) && (M.n_elem == 0);
  if (!skip_M && (M.n_rows != C1.n_rows || M.n_cols != C2.n_rows)) {
    Rcpp::stop("`M` must have shape nrow(C1) x nrow(C2).");
  }
  if (p.n_elem != C1.n_rows || q.n_elem != C2.n_rows) {
    Rcpp::stop("`p` and `q` must match `C1` and `C2`.");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }

  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const double lin_w = 1.0 - alpha;
  const double quad_w = alpha;
  const bool use_lin = !skip_M && lin_w != 0.0;

  arma::mat G = init_plan;
  if (G.n_elem == 0) {
    const double denom = arma::accu(p) * arma::accu(q);
    G = (p * q.t()) * (m / denom);
  } else if (G.n_rows != ns || G.n_cols != nt) {
    Rcpp::stop("`init_plan` has incompatible shape.");
  }

  GwSquareTermsResult gw = gw_square_terms_exact(C1, C2, G, symmetric);
  double lin_loss = use_lin ? lin_w * arma::accu(M % G) : 0.0;
  double cost = lin_loss + quad_w * gw.loss;

  std::vector<double> loss_trace;
  loss_trace.reserve(static_cast<std::size_t>(max_iter) + 1u);
  loss_trace.push_back(cost);

  double rel_delta = std::numeric_limits<double>::infinity();
  double abs_delta = std::numeric_limits<double>::infinity();
  int it = 0;

  for (int k = 0; k < max_iter; ++k) {
    arma::mat Mi = quad_w * gw.grad;
    if (use_lin) {
      Mi += lin_w * M;
    }
    const arma::mat Gc = partial_transport_direction(
      Mi, p, q, m, nb_dummies, lp_max_iter, lp_tol
    );
    const arma::mat delta = Gc - G;
    const GwSquareTermsResult gw_c = gw_square_terms_exact(C1, C2, Gc, symmetric);
    const arma::mat grad_delta = gw_c.grad - gw.grad;

    const double a_ls = quad_w * 0.5 * arma::accu(grad_delta % delta);
    double b_ls = quad_w * arma::accu(gw.grad % delta);
    if (use_lin) {
      b_ls += lin_w * arma::accu(M % delta);
    }
    double step = solve_1d_linesearch_quad(a_ls, b_ls);
    step = std::min(1.0, std::max(0.0, step));
    const double new_cost = cost + a_ls * step * step + b_ls * step;

    G += step * delta;
    gw.grad += step * grad_delta;

    abs_delta = std::abs(new_cost - cost);
    rel_delta = abs_delta / (std::abs(new_cost) + 1e-15);
    cost = new_cost;
    loss_trace.push_back(cost);
    it = k + 1;
    if (rel_delta <= tol) {
      break;
    }
  }

  gw = gw_square_terms_exact(C1, C2, G, symmetric);
  lin_loss = use_lin ? lin_w * arma::accu(M % G) : 0.0;
  const double quad_loss = quad_w * gw.loss;

  return Rcpp::List::create(
    Rcpp::Named("plan") = G,
    Rcpp::Named("objective") = lin_loss + quad_loss,
    Rcpp::Named("lin_loss") = lin_loss,
    Rcpp::Named("quad_loss") = quad_loss,
    Rcpp::Named("gw_loss") = gw.loss,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = rel_delta,
    Rcpp::Named("abs_error") = abs_delta,
    Rcpp::Named("loss_trace") = loss_trace
  );
}

// [[Rcpp::export]]
Rcpp::List cpp_partial_fgw_entropic_square(
    const arma::mat& M,
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::vec& p,
    const arma::vec& q,
    double m,
    double reg,
    double alpha,
    bool symmetric,
    const arma::mat& init_plan,
    int max_iter,
    double tol,
    int inner_max_iter,
    double inner_tol,
    int check_every) {
  if (C1.n_rows != C1.n_cols || C2.n_rows != C2.n_cols) {
    Rcpp::stop("`C1` and `C2` must be square.");
  }
  const bool skip_M = (alpha >= 1.0) && (M.n_elem == 0);
  if (!skip_M && (M.n_rows != C1.n_rows || M.n_cols != C2.n_rows)) {
    Rcpp::stop("`M` must have shape nrow(C1) x nrow(C2).");
  }
  if (p.n_elem != C1.n_rows || q.n_elem != C2.n_rows) {
    Rcpp::stop("`p` and `q` must match `C1` and `C2`.");
  }
  if (!std::isfinite(reg) || reg <= 0.0) {
    Rcpp::stop("`reg` must be positive.");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }
  if (check_every <= 0) {
    check_every = 1;
  }

  const arma::uword ns = C1.n_rows;
  const arma::uword nt = C2.n_rows;
  const double lin_w = 1.0 - alpha;
  const double quad_w = alpha;
  const bool use_lin = !skip_M && lin_w != 0.0;

  arma::mat G = init_plan;
  if (G.n_elem == 0) {
    const double denom = arma::accu(p) * arma::accu(q);
    G = (p * q.t()) * (m / denom);
  } else if (G.n_rows != ns || G.n_cols != nt) {
    Rcpp::stop("`init_plan` has incompatible shape.");
  }

  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  std::vector<double> err_trace;

  for (int k = 0; k < max_iter; ++k) {
    arma::mat Gprev;
    const bool do_check = ((k + 1) % check_every == 0) || (k == 0);
    if (do_check) {
      Gprev = G;
    }
    const GwSquareTermsResult gw = gw_square_terms_exact(C1, C2, G, symmetric);
    arma::mat M_entr = quad_w * gw.grad;
    if (use_lin) {
      M_entr += lin_w * M;
    }
    const EntropicPartialOtResult inner = entropic_partial_wasserstein_core(
      p, q, M_entr, reg, m, inner_max_iter, inner_tol, false, false
    );
    G = inner.plan;
    if (do_check) {
      err = std::sqrt(arma::accu((G - Gprev) % (G - Gprev)));
      err_trace.push_back(err);
      if (err <= tol) {
        it = k + 1;
        break;
      }
    }
    it = k + 1;
  }

  const GwSquareTermsResult gw_end = gw_square_terms_exact(C1, C2, G, symmetric);
  const double lin_loss = use_lin ? lin_w * arma::accu(M % G) : 0.0;
  const double quad_loss = quad_w * gw_end.loss;

  return Rcpp::List::create(
    Rcpp::Named("plan") = G,
    Rcpp::Named("objective") = lin_loss + quad_loss,
    Rcpp::Named("lin_loss") = lin_loss,
    Rcpp::Named("quad_loss") = quad_loss,
    Rcpp::Named("gw_loss") = gw_end.loss,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = err,
    Rcpp::Named("err_trace") = err_trace
  );
}

inline std::pair<double, double> ucoot_cost_kl(
    const arma::mat& X,
    const arma::mat& Y,
    const arma::mat& X_sqr,
    const arma::mat& Y_sqr,
    const arma::mat& M_samp,
    const arma::mat& M_feat,
    const arma::vec& wx_samp,
    const arma::vec& wy_samp,
    const arma::vec& wx_feat,
    const arma::vec& wy_feat,
    const arma::mat& wxy_samp,
    const arma::mat& wxy_feat,
    const arma::mat& pi_samp,
    const arma::mat& pi_feat,
    double rho_x,
    double rho_y,
    double eps_samp,
    double eps_feat,
    bool joint) {
  const arma::vec pi1_samp = arma::sum(pi_samp, 1);
  const arma::vec pi2_samp = arma::sum(pi_samp, 0).t();
  const arma::vec pi1_feat = arma::sum(pi_feat, 1);
  const arma::vec pi2_feat = arma::sum(pi_feat, 0).t();

  const double A_sqr = arma::dot(X_sqr * pi1_feat, pi1_samp);
  const double B_sqr = arma::dot(Y_sqr * pi2_feat, pi2_samp);
  const arma::mat AB = (X * pi_feat * Y.t()) % pi_samp;
  const double linear_cost = A_sqr + B_sqr - 2.0 * arma::accu(AB);

  double ucoot_cost = linear_cost;
  if (M_samp.n_elem > 0) {
    ucoot_cost += arma::accu(pi_samp % M_samp);
  }
  if (M_feat.n_elem > 0) {
    ucoot_cost += arma::accu(pi_feat % M_feat);
  }
  if (std::isfinite(rho_x) && rho_x != 0.0) {
    ucoot_cost += rho_x * div_between_product_kl(pi1_samp, pi1_feat, wx_samp, wx_feat);
  }
  if (std::isfinite(rho_y) && rho_y != 0.0) {
    ucoot_cost += rho_y * div_between_product_kl(pi2_samp, pi2_feat, wy_samp, wy_feat);
  }
  if (joint) {
    if (eps_samp != 0.0) {
      ucoot_cost += eps_samp * div_between_product_kl(pi_samp, pi_feat, wxy_samp, wxy_feat);
    }
  } else {
    if (eps_samp != 0.0) {
      ucoot_cost += eps_samp * div_to_product_kl(pi_samp, wx_samp, wy_samp, pi1_samp, pi2_samp, true);
    }
    if (eps_feat != 0.0) {
      ucoot_cost += eps_feat * div_to_product_kl(pi_feat, wx_feat, wy_feat, pi1_feat, pi2_feat, true);
    }
  }
  return {linear_cost, ucoot_cost};
}

// [[Rcpp::export]]
Rcpp::List cpp_ucoot_kl(
    const arma::mat& X,
    const arma::mat& Y,
    const arma::vec& wx_samp,
    const arma::vec& wx_feat,
    const arma::vec& wy_samp,
    const arma::vec& wy_feat,
    const arma::vec& reg_marginals,
    const arma::vec& epsilon,
    const arma::mat& M_samp,
    const arma::mat& M_feat,
    const arma::mat& init_pi_samp,
    const arma::mat& init_pi_feat,
    bool joint,
    bool rescale_plan,
    int max_iter,
    double tol,
    int max_iter_ot,
    double tol_ot,
    bool use_warm_start) {
  if (X.n_rows != wx_samp.n_elem || X.n_cols != wx_feat.n_elem) {
    Rcpp::stop("`wx_samp` / `wx_feat` must match `X`.");
  }
  if (Y.n_rows != wy_samp.n_elem || Y.n_cols != wy_feat.n_elem) {
    Rcpp::stop("`wy_samp` / `wy_feat` must match `Y`.");
  }
  if (reg_marginals.n_elem < 2 || epsilon.n_elem < 2) {
    Rcpp::stop("`reg_marginals` and `epsilon` must have length 2.");
  }
  if (max_iter < 1) {
    Rcpp::stop("`max_iter` must be >= 1.");
  }
  if (max_iter_ot < 1) {
    Rcpp::stop("`max_iter_ot` must be >= 1.");
  }

  const arma::uword nx_samp = X.n_rows;
  const arma::uword nx_feat = X.n_cols;
  const arma::uword ny_samp = Y.n_rows;
  const arma::uword ny_feat = Y.n_cols;
  const double rho_x = reg_marginals(0);
  const double rho_y = reg_marginals(1);
  double eps_samp = epsilon(0);
  double eps_feat = epsilon(1);
  if (joint) {
    eps_feat = eps_samp;
  }
  if (!(eps_samp > 0.0) || !(eps_feat > 0.0)) {
    Rcpp::stop("KL Sinkhorn UCOOT requires positive `epsilon` values.");
  }

  arma::mat M_samp_use = M_samp;
  arma::mat M_feat_use = M_feat;
  if (M_samp_use.n_elem == 0) {
    M_samp_use.zeros(nx_samp, ny_samp);
  } else if (M_samp_use.n_rows != nx_samp || M_samp_use.n_cols != ny_samp) {
    Rcpp::stop("`M_samp` has incompatible shape.");
  }
  if (M_feat_use.n_elem == 0) {
    M_feat_use.zeros(nx_feat, ny_feat);
  } else if (M_feat_use.n_rows != nx_feat || M_feat_use.n_cols != ny_feat) {
    Rcpp::stop("`M_feat` has incompatible shape.");
  }

  const arma::mat X_sqr = X % X;
  const arma::mat Y_sqr = Y % Y;
  const arma::mat Xt = X.t();
  const arma::mat Yt = Y.t();
  const arma::mat X_sqr_t = X_sqr.t();
  const arma::mat Y_sqr_t = Y_sqr.t();
  const arma::mat wxy_samp = wx_samp * wy_samp.t();
  const arma::mat wxy_feat = wx_feat * wy_feat.t();
  const arma::vec log_wx_samp = arma::log(wx_samp + kTiny);
  const arma::vec log_wy_samp = arma::log(wy_samp + kTiny);
  const arma::vec log_wx_feat = arma::log(wx_feat + kTiny);
  const arma::vec log_wy_feat = arma::log(wy_feat + kTiny);
  const double eps_cost_feat = joint ? eps_samp : 0.0;
  const double eps_cost_samp = joint ? eps_feat : 0.0;

  arma::mat pi_samp = init_pi_samp;
  arma::mat pi_feat = init_pi_feat;
  if (pi_samp.n_elem == 0) {
    pi_samp = wxy_samp;
  } else if (pi_samp.n_rows != nx_samp || pi_samp.n_cols != ny_samp) {
    Rcpp::stop("`init_pi_samp` has incompatible shape.");
  }
  if (pi_feat.n_elem == 0) {
    pi_feat = wxy_feat;
  } else if (pi_feat.n_rows != nx_feat || pi_feat.n_cols != ny_feat) {
    Rcpp::stop("`init_pi_feat` has incompatible shape.");
  }

  SinkhornUnbalancedWorkspace ws_samp;
  SinkhornUnbalancedWorkspace ws_feat;
  arma::mat uot_cost_feat(nx_feat, ny_feat);
  arma::mat scratch_feat(nx_feat, ny_samp);
  arma::mat uot_cost_samp(nx_samp, ny_samp);
  arma::mat scratch_samp(nx_samp, ny_feat);
  arma::vec pi1_s(nx_samp);
  arma::vec pi2_s(ny_samp);
  arma::vec A_feat(nx_feat);
  arma::vec B_feat(ny_feat);
  arma::vec pi1_f(nx_feat);
  arma::vec pi2_f(ny_feat);
  arma::vec A_samp(nx_samp);
  arma::vec B_samp(ny_samp);

  std::vector<int> inner_iters_feat;
  std::vector<int> inner_iters_samp;
  std::vector<int> inner_warm_feat;
  std::vector<int> inner_warm_samp;
  std::vector<int> inner_fallback_feat;
  std::vector<int> inner_fallback_samp;
  inner_iters_feat.reserve(static_cast<std::size_t>(max_iter));
  inner_iters_samp.reserve(static_cast<std::size_t>(max_iter));
  inner_warm_feat.reserve(static_cast<std::size_t>(max_iter));
  inner_warm_samp.reserve(static_cast<std::size_t>(max_iter));
  inner_fallback_feat.reserve(static_cast<std::size_t>(max_iter));
  inner_fallback_samp.reserve(static_cast<std::size_t>(max_iter));

  Rcpp::NumericVector err_trace;
  double err = std::numeric_limits<double>::infinity();
  int it = 0;
  double feat_ms = 0.0;
  double samp_ms = 0.0;

  for (; it < max_iter; ++it) {
    const auto t_feat0 = std::chrono::steady_clock::now();
    const double mass_samp = arma::accu(pi_samp);
    uot_cost_matrix_kl_joint_inplace(
      X_sqr_t, Y_sqr_t, Xt, Y, M_feat_use, pi_samp,
      wx_samp, wy_samp, log_wx_samp, log_wy_samp,
      rho_x, rho_y, eps_cost_feat,
      uot_cost_feat, scratch_feat, pi1_s, pi2_s, A_feat, B_feat
    );
    const double eps_feat_ot = joint ? (mass_samp * eps_feat) : eps_feat;
    SinkhornUnbalancedResult su_feat = sinkhorn_unbalanced_kl(
      uot_cost_feat, wx_feat, wy_feat, wxy_feat,
      rho_x * mass_samp, rho_y * mass_samp, eps_feat_ot,
      max_iter_ot, tol_ot, pi_feat, ws_feat, use_warm_start
    );
    inner_iters_feat.push_back(su_feat.iters);
    inner_warm_feat.push_back(su_feat.warm_started ? 1 : 0);
    inner_fallback_feat.push_back(su_feat.warm_fallback ? 1 : 0);
    pi_feat = std::move(su_feat.plan);
    if (rescale_plan) {
      const double mass_feat = arma::accu(pi_feat);
      if (mass_feat > 0.0) {
        pi_feat *= std::sqrt(mass_samp / mass_feat);
      }
    }
    feat_ms += std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - t_feat0
    ).count();

    const auto t_samp0 = std::chrono::steady_clock::now();
    const double mass_feat = arma::accu(pi_feat);
    uot_cost_matrix_kl_joint_inplace(
      X_sqr, Y_sqr, X, Yt, M_samp_use, pi_feat,
      wx_feat, wy_feat, log_wx_feat, log_wy_feat,
      rho_x, rho_y, eps_cost_samp,
      uot_cost_samp, scratch_samp, pi1_f, pi2_f, A_samp, B_samp
    );
    const double eps_samp_ot = joint ? (mass_feat * eps_feat) : eps_feat;
    SinkhornUnbalancedResult su_samp = sinkhorn_unbalanced_kl(
      uot_cost_samp, wx_samp, wy_samp, wxy_samp,
      rho_x * mass_feat, rho_y * mass_feat, eps_samp_ot,
      max_iter_ot, tol_ot, pi_samp, ws_samp, use_warm_start
    );
    inner_iters_samp.push_back(su_samp.iters);
    inner_warm_samp.push_back(su_samp.warm_started ? 1 : 0);
    inner_fallback_samp.push_back(su_samp.warm_fallback ? 1 : 0);
    arma::mat pi_next = std::move(su_samp.plan);
    if (rescale_plan) {
      const double mass_next = arma::accu(pi_next);
      if (mass_next > 0.0) {
        pi_next *= std::sqrt(mass_feat / mass_next);
      }
    }
    err = arma::accu(arma::abs(pi_next - pi_samp));
    pi_samp = std::move(pi_next);
    samp_ms += std::chrono::duration<double, std::milli>(
      std::chrono::steady_clock::now() - t_samp0
    ).count();

    err_trace.push_back(err);
    if (err < tol) {
      ++it;
      break;
    }
  }

  if (!pi_samp.is_finite() || !pi_feat.is_finite()) {
    Rcpp::stop("Encountered non-finite values in coupling matrices. Adjust hyperparameters.");
  }

  const auto costs = ucoot_cost_kl(
    X, Y, X_sqr, Y_sqr, M_samp_use, M_feat_use,
    wx_samp, wy_samp, wx_feat, wy_feat, wxy_samp, wxy_feat,
    pi_samp, pi_feat, rho_x, rho_y, eps_samp, eps_feat, joint
  );
  const int inner_total = std::accumulate(inner_iters_feat.begin(), inner_iters_feat.end(), 0) +
    std::accumulate(inner_iters_samp.begin(), inner_iters_samp.end(), 0);

  return Rcpp::List::create(
    Rcpp::Named("pi_samp") = pi_samp,
    Rcpp::Named("pi_feat") = pi_feat,
    Rcpp::Named("linear_cost") = costs.first,
    Rcpp::Named("ucoot_cost") = costs.second,
    Rcpp::Named("iterations") = it,
    Rcpp::Named("error") = err,
    Rcpp::Named("err_trace") = err_trace,
    Rcpp::Named("inner_iters_feat") = Rcpp::IntegerVector(inner_iters_feat.begin(), inner_iters_feat.end()),
    Rcpp::Named("inner_iters_samp") = Rcpp::IntegerVector(inner_iters_samp.begin(), inner_iters_samp.end()),
    Rcpp::Named("inner_warm_feat") = Rcpp::LogicalVector(inner_warm_feat.begin(), inner_warm_feat.end()),
    Rcpp::Named("inner_warm_samp") = Rcpp::LogicalVector(inner_warm_samp.begin(), inner_warm_samp.end()),
    Rcpp::Named("inner_warm_fallback_feat") = Rcpp::LogicalVector(inner_fallback_feat.begin(), inner_fallback_feat.end()),
    Rcpp::Named("inner_warm_fallback_samp") = Rcpp::LogicalVector(inner_fallback_samp.begin(), inner_fallback_samp.end()),
    Rcpp::Named("inner_iters_total") = inner_total,
    Rcpp::Named("feat_ms") = feat_ms,
    Rcpp::Named("samp_ms") = samp_ms,
    Rcpp::Named("warm_start") = use_warm_start
  );
}
