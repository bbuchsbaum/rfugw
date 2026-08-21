// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include "gw_square.h"
#include <utility>

// This translation unit is the R adapter for canonical square-loss algebra.
// The actual forward/reverse tensor, loss, and gradient contract lives in the
// R-independent header and is shared by symmetric and general paths.

namespace rfugw {

arma::mat gw_square_forward_tensor(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    const arma::vec& pG,
    const arma::vec& qG) {
  const arma::vec source = (C1 % C1) * pG;
  const arma::vec target = (C2 % C2) * qG;
  arma::mat tensor = source * arma::ones<arma::rowvec>(C2.n_rows);
  tensor.each_row() += target.t();
  tensor -= 2.0 * C1 * G * C2.t();
  return tensor;
}

arma::mat gw_square_reverse_tensor(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    const arma::vec& pG,
    const arma::vec& qG) {
  const arma::mat C1t = C1.t();
  const arma::vec source = (C1t % C1t) * pG;
  const arma::vec target = (C2 % C2).t() * qG;
  arma::mat tensor = source * arma::ones<arma::rowvec>(C2.n_rows);
  tensor.each_row() += target.t();
  tensor -= 2.0 * C1t * G * C2;
  return tensor;
}

double gw_square_loss(
    const arma::mat& forward,
    const arma::mat& reverse,
    const arma::mat& G,
    bool symmetric) {
  const double forward_loss = arma::accu(forward % G);
  if (symmetric) {
    return forward_loss;
  }
  return 0.5 * (forward_loss + arma::accu(reverse % G));
}

arma::mat gw_square_gradient(
    const arma::mat& forward,
    const arma::mat& reverse,
    bool symmetric) {
  if (symmetric) {
    return 2.0 * forward;
  }
  return forward + reverse;
}

GwSquareTerms gw_square_terms(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    bool symmetric) {
  const arma::vec pG = arma::sum(G, 1);
  const arma::vec qG = arma::sum(G, 0).t();
  arma::mat forward = gw_square_forward_tensor(C1, C2, G, pG, qG);
  arma::mat reverse = symmetric
    ? forward
    : gw_square_reverse_tensor(C1, C2, G, pG, qG);

  GwSquareTerms out;
  out.loss = gw_square_loss(forward, reverse, G, symmetric);
  out.gradient = gw_square_gradient(forward, reverse, symmetric);
  out.forward = std::move(forward);
  out.reverse = std::move(reverse);
  return out;
}

}  // namespace rfugw

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

  rfugw::GwSquareTerms out = rfugw::gw_square_terms(C1, C2, G, symmetric);
  return Rcpp::List::create(
    Rcpp::Named("loss") = out.loss,
    Rcpp::Named("grad") = out.gradient,
    Rcpp::Named("forward_tensor") = out.forward,
    Rcpp::Named("reverse_tensor") = out.reverse
  );
}
