// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppHungarian)]]

#include <RcppArmadillo.h>
#include <RcppHungarian.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List hungarian_arma(const arma::mat cost) {
  size_t n_rows = cost.n_rows;
  size_t n_cols = cost.n_cols;

  std::vector<std::vector<double>> DistMatrix(n_rows, std::vector<double>(n_cols));
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      DistMatrix[i][j] = cost(i, j);
    }
  }

  std::vector<int> Assignment;
  HungarianAlgorithm hungarian;

  double total_cost = hungarian.Solve(DistMatrix, Assignment);

  IntegerMatrix pairs(n_rows, 2);
  for (size_t i = 0; i < n_rows; ++i) {
    pairs(i, 0) = i + 1;
    pairs(i, 1) = Assignment[i] + 1;
  }

  return List::create(Named("cost")  = total_cost,
                      Named("pairs") = pairs);
}
