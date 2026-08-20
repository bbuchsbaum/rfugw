#ifndef RFUGW_APPROXIMATION_CACHE_H
#define RFUGW_APPROXIMATION_CACHE_H

#include <armadillo>

namespace rfugw {

// Precision-generic cache contracts. Construction can remain backend-specific,
// but float and double paths share dimensions, validity, and C1 preparation.
template <typename T>
struct LowRankC2Cache {
  bool valid = false;
  arma::Mat<T> A2t_scaled;
  arma::Mat<T> B2t;
};

using LowRankC2CacheD = LowRankC2Cache<double>;
using LowRankC2CacheF = LowRankC2Cache<float>;

template <typename T>
inline bool lowrank_cache_compatible_t(
    const LowRankC2Cache<T>& cache,
    const arma::Mat<T>& C2) {
  return cache.valid &&
    cache.A2t_scaled.n_cols == C2.n_rows &&
    cache.B2t.n_cols == C2.n_cols &&
    cache.A2t_scaled.n_rows == cache.B2t.n_rows;
}

template <typename T>
inline bool lowrank_c1_cache_compatible_t(
    const arma::Mat<T>& A1_scaled,
    const arma::Mat<T>& B1t,
    const arma::Mat<T>& C1) {
  return A1_scaled.n_rows == C1.n_rows &&
    B1t.n_cols == C1.n_cols &&
    A1_scaled.n_cols == B1t.n_rows &&
    A1_scaled.n_cols > 0;
}

template <typename T>
struct SquareC2CacheBase {
  bool valid = false;
  arma::Mat<T> hC2;
  arma::Row<T> right_term;
};

struct SquareC2CacheD : SquareC2CacheBase<double> {};

struct SquareC2CacheF : SquareC2CacheBase<float> {
  arma::fmat C2;
  arma::fvec q;
};

template <typename InputT, typename CacheT>
inline bool square_cache_dimensions_compatible_t(
    const CacheT& cache,
    const arma::Mat<InputT>& C2,
    const arma::Col<InputT>& q) {
  return cache.valid &&
    C2.n_rows == C2.n_cols &&
    q.n_elem == C2.n_rows &&
    cache.hC2.n_rows == C2.n_rows &&
    cache.hC2.n_cols == C2.n_cols &&
    cache.right_term.n_elem == C2.n_rows;
}

template <typename T, typename CacheT>
inline void init_matrices_square_from_cache_t(
    const arma::Mat<T>& C1,
    const arma::Col<T>& p,
    const CacheT& cache,
    arma::Mat<T>& constC,
    arma::Mat<T>& hC1) {
  const arma::Mat<T> fC1 = C1 % C1;
  hC1 = C1;
  const arma::Col<T> left = fC1 * p;
  constC.set_size(C1.n_rows, cache.hC2.n_rows);
  for (arma::uword j = 0; j < cache.hC2.n_rows; ++j) {
    constC.col(j) = left;
    constC.col(j) += cache.right_term[j];
  }
}

}  // namespace rfugw

#endif
