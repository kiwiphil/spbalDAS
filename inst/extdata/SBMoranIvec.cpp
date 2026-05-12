
// =============================================================================
// SBMoranIvec.cpp
// Fast RcppArmadillo equivalent of the original MATLAB / R function
// =============================================================================

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Fast C++ implementation of modified Moran's I sample spread measure
 //'
 //' @param C Integer vector - current sample indices (1-based)
 //' @param A Integer vector - candidate indices to evaluate (1-based)
 //' @param W NumericMatrix - square spatial weight matrix
 //'
 //' @return Numeric vector of length(A) with the spread measure for each candidate
 // [[Rcpp::export]]
 Rcpp::NumericVector xxxSBMoranIvec(Rcpp::IntegerVector C,
                                 Rcpp::IntegerVector A,
                                 Rcpp::NumericMatrix W) {

   const arma::uword p = C.size();          // current sample size

   // 0-based indices (Armadillo uvec)
   arma::uvec C0 = Rcpp::as<arma::uvec>(C) - 1;
   arma::uvec A0 = Rcpp::as<arma::uvec>(A) - 1;

   // Zero-copy view of the weight matrix (no data duplication)
   const arma::mat W_arma(W.begin(), W.nrow(), W.ncol(), false);

   // Sum of all weights inside the current sample C
   double WCij = (p == 0) ? 0.0 : arma::accu(W_arma.submat(C0, C0));

   // Extract sub-matrices and compute row/column sums in one go
   arma::mat AC = W_arma.submat(A0, C0);           // |A| x |C|
   arma::mat CA = W_arma.submat(C0, A0);           // |C| x |A|

   arma::vec contrib = arma::sum(AC, 1) + arma::sum(CA, 0).t();

   return Rcpp::wrap(WCij + contrib);
 }
