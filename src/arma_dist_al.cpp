
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]


//' @name arma_dist_al
//'
//' @title Calculate a distance matrix.
//'
//' @description Computes a distance matrix for the supplied matrix of n-dimensional points
//'              using Euclidean distance.
//'
//' @author Phil Davies.
//'
//' @param X An n-dimensional matrix of points.
//'
//' @return The distance matrix for the point matrix X.
//'
//' @examples
//' # Distance matrix for the two 3D points.
//' spbalDAS::arma_dist_al(X = rbind(c(1, 2, 3), c(4, 3, 2)))
//' # Distance matrix for three 2D points.
//' spbalDAS::arma_dist_al(X = matrix(c(2, 3, 0, 9, 4, 5), nrow = 3, byrow = TRUE))
//'
//' @export
// [[Rcpp::export(rng = false)]]
arma::mat arma_dist_al(const arma::mat& X) {
  int n = X.n_rows;

  // Compute squared norms of rows: ||x_i||^2
  arma::vec G = arma::sum(arma::square(X), 1); // Sum of squares along rows

  // Compute D^2 = G + G^T - 2 * X * X^T
  arma::mat D = -2.0 * X * X.t();  // -2 * dot products
  D.each_col() += G;               // Add G to each column
  D.each_row() += G.t();           // Add G^T to each row

  // Take square root and ensure non-negative values
  D = arma::sqrt(arma::abs(D));

  return D;
}
