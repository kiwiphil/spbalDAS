#include <RcppHungarian.h>
// [[Rcpp::depends(RcppHungarian)]] // Necessary if using sourceCpp or not building a package

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>




#include <random>
#include <vector>
#include <numeric>


using namespace Rcpp;
using namespace arma;

// Forward declarations
arma::vec moran_contribution_vec(const arma::uvec& C,
                                 const arma::uvec& A_candidates,
                                 const arma::mat& W);

arma::mat solve_linear_assignment(const arma::mat& cost_matrix);

arma::uvec arma_setdiff(arma::uvec x, arma::uvec y);

List hungarian_arma(arma::mat Cost);
//List HungarianSolver;
//HungarianAlgorithm solver;

std::vector<int> custom_sample(int n, int k, bool replace = false) {
  std::vector<int> pool(n);
  std::iota(pool.begin(), pool.end(), 0);  // 0 .. n-1

  std::mt19937 rng{std::random_device{}()};   // or your own seed

  if (replace) {
    std::vector<int> result(k);
    std::uniform_int_distribution<int> dist(0, n-1);
    for(int &v : result) v = dist(rng);
    return result;
  } else {
    // Fisher-Yates shuffle + take first k
    std::shuffle(pool.begin(), pool.end(), rng);
    return {pool.begin(), pool.begin() + k};
  }
}


//' Dynamic Assignment Sampling (DAS) - Rcpp/Armadillo version
//'
//' @param pop n × d numeric matrix of point coordinates
//' @param n desired final sample size
//' @param verbose print progress messages
//'
//' @return integer matrix (growing to n columns) containing indices (1-based)
// [[Rcpp::export]]
 IntegerMatrix DAS_cpp(const NumericMatrix& pop_,
                       int n,
                       bool verbose = false)
 {
   arma::mat pop = Rcpp::as<arma::mat>(pop_);
   int N = pop.n_rows;

   if (n >= N) {
     stop("Number of samples must be smaller than population size N");
   }

   // Initial number of candidates per "cluster"
   int J1 = std::max(2, static_cast<int>(N / 2));
   std::vector<int> Jvec = {J1, J1};

   // Random permutation of 1..N (Rcpp style)
   IntegerVector idx = seq(1, N);
   idx = sample(idx, N, false);                    // random permutation

   // Take first J1 elements (1-based indices)
   IntegerVector initial = idx[Range(0, J1 - 1)];

   // === FIXED: Proper way to create initial SampleMatrix ===
   // Option A - cleanest and most common
   arma::Mat<int> SampleMatrix = arma::conv_to<arma::Mat<int>>::from(
     Rcpp::as<arma::uvec>(initial)
   );

   // Option B - more explicit (also correct)
   // arma::Mat<int> SampleMatrix(J1, 1);
   // SampleMatrix.col(0) = arma::conv_to<arma::Col<int>>::from(initial);

   // Option C - zero-copy style (very efficient)
   // arma::Mat<int> SampleMatrix(initial.begin(), J1, 1, false);

   if (verbose) {
     Rcout << "Initial candidates (" << J1 << " points):\n"
           << SampleMatrix << "\n";
   }

   // ===============================================
   // Compute spatial weight matrix W = 1 / distance
   // ===============================================
   //arma::mat D = arma::sqrt(arma::sum(arma::square(pop.each_row() - pop.each_col()), 1));
   //arma::mat W = 1.0 / (D + arma::datum::eps);   // avoid division by zero
   //W.diag().zeros();                             // usually no self-distance
   //int n = pop.n_rows;

   // Compute squared norms of rows: ||x_i||^2
   arma::vec G = arma::sum(arma::square(pop), 1); // Sum of squares along rows

   // Compute D^2 = G + G^T - 2 * X * X^T
   arma::mat D = -2.0 * pop * pop.t(); // -2 * dot products
   D.each_col() += G;               // Add G to each column
   D.each_row() += G.t();           // Add G^T to each row

   // Take square root and ensure non-negative values
   D = arma::sqrt(arma::abs(D));

   arma::mat W = 1.0 / (D + arma::datum::eps);   // avoid division by zero
   W.diag().zeros();

   Rcout << "Have W.\n";

   // Main iterative loop
   for (int i = 2; i <= n; ++i) {
     int curr_J = Jvec[i-2];   // J[i-1] in 1-based R logic

     // Randomly select curr_J rows from current SampleMatrix
     IntegerVector scrs = seq(0, SampleMatrix.n_rows-1);
     IntegerVector row_idx = sample(scrs, curr_J, false);
     //Rcpp::RcppArmadillo::sample(seq(1, N), k)
     arma::Mat<int> Ci(curr_J, SampleMatrix.n_cols);

     for (int k = 0; k < curr_J; ++k) {
       Ci.row(k) = SampleMatrix.row(row_idx[k]);
     }

     // All points currently used
     arma::uvec used = arma::unique(arma::vectorise(arma::conv_to<arma::umat>::from(Ci)));
     Rcout << "in main loop.\n";
     // Remaining points Ω \ Ci
     arma::uvec all_points = arma::regspace<arma::uvec>(1, N);
     arma::uvec OmegaNotCi = arma_setdiff(all_points, used + 1); // +1 because 1-based

     int n_remaining = OmegaNotCi.n_elem;
     int n_to_take = std::min(curr_J, n_remaining);

     if (n_to_take < curr_J && verbose) {
       Rcout << "Warning: only " << n_to_take << " remaining points available at step " << i << " n_remaining " << n_remaining << "\n";
     }

     // Random sample from remaining
     arma::uvec perm = arma::randperm(n_remaining, n_to_take);
     arma::uvec Ai = OmegaNotCi.elem(perm);

     // Create CAi = [Ci | Ai]  (Ci has many columns, Ai one column)
     arma::Mat<int> CAi(curr_J, Ci.n_cols + 1);
     CAi.cols(0, Ci.n_cols-1) = Ci;
     CAi.col(Ci.n_cols) = arma::conv_to<arma::Col<int>>::from(Ai);

     // ===============================================
     // Compute cost matrix using Moran-like contribution
     // ===============================================
     arma::vec internal = moran_contribution_vec(
       arma::conv_to<arma::uvec>::from(arma::vectorise(Ci.row(0))), // dummy - just for type
       arma::conv_to<arma::uvec>::from(arma::vectorise(CAi.col(CAi.n_cols-1))),
       W
     );

     arma::mat Cost(curr_J, curr_J, arma::fill::zeros);

     for (int r = 0; r < curr_J; ++r) {
       arma::uvec C_r = arma::conv_to<arma::uvec>::from(arma::vectorise(Ci.row(r)));
       Cost.row(r) = moran_contribution_vec(C_r, Ai, W).t();
     }

     // Solve assignment problem (minimum cost assignment)
     //arma::mat assignment_result = solve_linear_assignment(Cost);
     List res = hungarian_arma(Cost);
     arma::mat assignment_result = res[1];

     // Reconstruct improved SampleMatrix
     arma::Mat<int> new_SampleMatrix = CAi;
     arma::uvec perm_indices = arma::conv_to<arma::uvec>::from(assignment_result.col(1) - 1); // 0-based

     new_SampleMatrix.col(new_SampleMatrix.n_cols-1) =
       arma::conv_to<arma::Col<int>>::from(Ai.elem(perm_indices));

     SampleMatrix = new_SampleMatrix;

     // Update number of candidates for next round
     int next_J = std::max(1, std::min(Jvec.back(), static_cast<int>(N / (i + 1))));
     Jvec.push_back(next_J);

     if (verbose && (i % 5 == 0 || i == n)) {
       Rcout << "Step " << i << " / " << n << "  -  candidates = " << next_J << "\n";
     }
   }

   return wrap(SampleMatrix);
 }

 // [[Rcpp::export]]
 arma::uvec arma_setdiff(arma::uvec x, arma::uvec y) {
   // Armadillo uvec (unsigned integer vectors) can be used directly with
   // std::set_difference, but must be sorted first.
   x = arma::sort(x);
   y = arma::sort(y);

   // The result will be stored in a std::vector first
   std::vector<unsigned int> result_vec;

   // Perform the set difference
   std::set_difference(
     x.begin(), x.end(),
     y.begin(), y.end(),
     std::back_inserter(result_vec)
   );

   // Convert the std::vector back to an arma::uvec and return
   arma::uvec result_arma(result_vec.size());
   for (size_t i = 0; i < result_vec.size(); ++i) {
     result_arma[i] = result_vec[i];
   }

   return result_arma;
 }


 // Helper: Moran-like interaction term (internal + incoming + outgoing)
 arma::vec moran_contribution_vec(const arma::uvec& C,
                                  const arma::uvec& A_candidates,
                                  const arma::mat& W)
 {
   if (C.is_empty()) return arma::vec(A_candidates.n_elem, arma::fill::zeros);

   double internal = arma::accu(W.submat(C, C));

   arma::vec out_sum  = arma::sum(W.submat(A_candidates, C), 1);
   arma::vec in_sum   = arma::sum(W.submat(C, A_candidates), 0).t();

   return internal + out_sum + in_sum;
   Rcout << "exit moran_contribution_vec.\n";
 }

 // Placeholder for linear assignment solver
 // You should replace this with real implementation:
 //   - RcppHungarian
 //   - lemon::MinCostFlow
 //   - LAPJV / Jonker-Volgenant
 //   - your own implementation
 arma::mat solve_linear_assignment(const arma::mat& cost_matrix)
 {
   Rcout << "entry solve_linear_assignment.\n";
   // Very naive greedy approach — REPLACE THIS!
   int m = cost_matrix.n_rows;
   arma::mat result(m, 2);

   arma::uvec perm = arma::sort_index(arma::min(cost_matrix, 1), "ascend");

   for (int i = 0; i < m; ++i) {
     result(i, 0) = i + 1;           // left side (1-based)
     result(i, 1) = perm(i) + 1;     // right side (1-based)
   }

   return result;
 }
