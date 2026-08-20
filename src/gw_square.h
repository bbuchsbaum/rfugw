#pragma once  // canonical square-loss GW algebra

#include <armadillo>
namespace rfugw {

struct GwSquareTerms {
  arma::mat forward;
  arma::mat reverse;
  double loss;
  arma::mat gradient;
};

arma::mat gw_square_forward_tensor(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    const arma::vec& pG,
    const arma::vec& qG);

arma::mat gw_square_reverse_tensor(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    const arma::vec& pG,
    const arma::vec& qG);

double gw_square_loss(
    const arma::mat& forward,
    const arma::mat& reverse,
    const arma::mat& G,
    bool symmetric);

arma::mat gw_square_gradient(
    const arma::mat& forward,
    const arma::mat& reverse,
    bool symmetric);

GwSquareTerms gw_square_terms(
    const arma::mat& C1,
    const arma::mat& C2,
    const arma::mat& G,
    bool symmetric);

}  // namespace rfugw
