
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]


// [[Rcpp::export]]
arma::mat arma_dist_al(const arma::mat& X) {
  int n = X.n_rows;

  // Compute squared norms of rows: ||x_i||^2
  arma::vec G = arma::sum(arma::square(X), 1); // Sum of squares along rows

  // Compute D^2 = G + G^T - 2 * X * X^T
  arma::mat D = -2.0 * X * X.t(); // -2 * dot products
  D.each_col() += G;               // Add G to each column
  D.each_row() += G.t();           // Add G^T to each row

  // Take square root and ensure non-negative values
  D = arma::sqrt(arma::abs(D));

  return D;
}
