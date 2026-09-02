// [[Rcpp::depends(RcppArmadillo, RcppThread)]]
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <chrono>
#include <vector>
#include <unordered_set>

// ScopedTimer (unchanged)
struct ScopedTimer {
  std::string label;
  std::chrono::high_resolution_clock::time_point start;
  //auto start_time = std::chrono::high_resolution_clock::now();
  ScopedTimer(std::string l) : label(std::move(l)),
  start(std::chrono::high_resolution_clock::now()) {}
  ~ScopedTimer() {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now() - start).count();
    RcppThread::Rcout << label << ": " << us << " \u00b5s\n";
  }
};

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <limits>

// ------------------------------------------------------------------
// Types
// ------------------------------------------------------------------
using cost = double;
using row  = int;
using col  = int;

const cost BIG = 1e30;

// ------------------------------------------------------------------
// Jonker-Volgenant LAPJV (square only)
// ------------------------------------------------------------------
//' @name native_lapjv
//'
//' @title Solve the Linear Assignment Problem.
//'
//' @description Solves the linear assignment problem using the Jonker-Volgenant algorithm.
//'
//' @author Phil Davies.
//'
//' @param Cost_R The cost matrix.
//'
//' @return A list containing two variables, cost and assignment.
//'
//' @examples
//' # Cost matrix for three 3D points. $cost = 1, $assignment = 2, 1, 0
//' spbalDAS::native_lapjv(Cost_R = rbind(c(1, 2, 0), c(2, 0, 1), c(1, 4, 19)))
//' # Cost matrix for three 3D points. $cost = 10, $assignment = 1, 0, 2
//' spbalDAS::native_lapjv(Cost_R = rbind(c(4, 2, 8), c(2, 3, 7), c(3, 1, 6)))
//' # Assign 4 machines to 4 jobs to minimize total setup time. $cost = 15, $assignment = 1, 3, 2, 0
//' spbalDAS::native_lapjv(Cost_R = rbind(c(14, 5, 8, 7), c(2, 12, 6, 5), c(7, 8, 3, 9), c(2, 4, 6, 10)))
//' # $cost = 55, $assignment = 2, 3, 4, 0, 1
//' spbalDAS::native_lapjv(Cost_R = rbind(c(10, 4, 6, 10, 12), c(11, 7, 7, 9, 14), c(13, 8, 12, 14, 15), c(14, 16, 13, 17, 17), c(17, 11, 17, 20, 19)))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List native_lapjv(Rcpp::NumericMatrix Cost_R) {

  arma::mat assigncost = Rcpp::as<arma::mat>(Cost_R);
  int dim = assigncost.n_rows;

  if (dim == 0 || dim != assigncost.n_cols) {
    Rcpp::stop("spbal(native_lapjv) requires a square cost matrix (rows == cols).\n"
                 "Your matrix is %d x %d.\n"
                 "In your DAS code C is always square, so this is fine.",
                 dim, assigncost.n_cols);
  }

  // Output arrays
  std::vector<col> rowsol(dim, -1);
  std::vector<row> colsol(dim, -1);
  std::vector<cost> u(dim, 0.0);
  std::vector<cost> v(dim, 0.0);

  // Working arrays
  std::vector<row> freeunassigned(dim);
  std::vector<col> collist(dim);
  std::vector<col> matches(dim, 0);
  std::vector<cost> d(dim);
  std::vector<row> pred(dim);

  int numfree = 0;

  // ====================== COLUMN REDUCTION ======================
  for (col j = dim; j--;) {
    cost minv = assigncost(0, j);
    row imin = 0;
    for (row i = 1; i < dim; ++i) {
      if (assigncost(i, j) < minv) {
        minv = assigncost(i, j);
        imin = i;
      }
    }
    v[j] = minv;
    if (++matches[imin] == 1) {
      rowsol[imin] = j;
      colsol[j] = imin;
    } else if (v[j] < v[rowsol[imin]]) {
      col j1 = rowsol[imin];
      rowsol[imin] = j;
      colsol[j] = imin;
      colsol[j1] = -1;
    } else {
      colsol[j] = -1;
    }
  }

  // ====================== REDUCTION TRANSFER ======================
  for (row i = 0; i < dim; ++i) {
    if (matches[i] == 0) {
      freeunassigned[numfree++] = i;
    } else if (matches[i] == 1) {
      col j1 = rowsol[i];
      cost minv = BIG;
      for (col j = 0; j < dim; ++j) {
        if (j != j1) {
          cost h = assigncost(i, j) - v[j];
          if (h < minv) minv = h;
        }
      }
      v[j1] -= minv;
    }
  }

  // ====================== AUGMENTING ROW REDUCTION (twice) ======================
  for (int loop = 0; loop < 2; ++loop) {
    int k = 0;
    int prvnumfree = numfree;
    numfree = 0;
    while (k < prvnumfree) {
      row i = freeunassigned[k++];
      cost umin = assigncost(i, 0) - v[0];
      col j1 = 0;
      cost usubmin = BIG;
      col j2 = 0;
      for (col j = 1; j < dim; ++j) {
        cost h = assigncost(i, j) - v[j];
        if (h < usubmin) {
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          } else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
        }
      }
      row i0 = colsol[j1];
      if (umin < usubmin) {
        v[j1] -= (usubmin - umin);
      } else if (i0 >= 0) {
        j1 = j2;
        i0 = colsol[j2];
      }
      rowsol[i] = j1;
      colsol[j1] = i;
      if (i0 >= 0) {
        if (umin < usubmin)
          freeunassigned[--k] = i0;
        else
          freeunassigned[numfree++] = i0;
      }
    }
  }

  // ====================== AUGMENT SOLUTION ======================
  for (int f = 0; f < numfree; ++f) {
    row freerow = freeunassigned[f];

    for (col j = 0; j < dim; ++j) {
      d[j] = assigncost(freerow, j) - v[j];
      pred[j] = freerow;
      collist[j] = j;
    }

    int low = 0;
    int up = 0;
    bool unassignedfound = false;
    col endofpath = 0;
    int last = 0;

    do {
      if (up == low) {
        last = low - 1;
        cost minv = d[collist[up++]];
        for (int k = up; k < dim; ++k) {
          col j = collist[k];
          cost h = d[j];
          if (h <= minv) {
            if (h < minv) {
              up = low;
              minv = h;
            }
            collist[k] = collist[up];
            collist[up++] = j;
          }
        }
        for (int k = low; k < up; ++k) {
          if (colsol[collist[k]] < 0) {
            endofpath = collist[k];
            unassignedfound = true;
            break;
          }
        }
      }
      if (!unassignedfound) {
        col j1 = collist[low++];
        row i = colsol[j1];
        cost h = assigncost(i, j1) - v[j1] - d[j1];
        for (int k = up; k < dim; ++k) {
          col j = collist[k];
          cost v2 = assigncost(i, j) - v[j] - h;
          if (v2 < d[j]) {
            pred[j] = i;
            if (v2 == d[j1]) {
              if (colsol[j] < 0) {
                endofpath = j;
                unassignedfound = true;
                break;
              } else {
                collist[k] = collist[up];
                collist[up++] = j;
              }
            }
            d[j] = v2;
          }
        }
      }
    } while (!unassignedfound);

    for (int k = last + 1; k--;) {
      col j1 = collist[k];
      v[j1] += d[j1] - d[collist[last]];
    }

    row i;
    do {
      i = pred[endofpath];
      colsol[endofpath] = i;
      col j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    } while (i != freerow);
  }

  // ====================== TOTAL COST ======================
  cost lapcost = 0.0;
  for (row i = 0; i < dim; ++i) {
    col j = rowsol[i];
    u[i] = assigncost(i, j) - v[j];
    lapcost += assigncost(i, j);
  }

  // ====================== RETURN ======================
  Rcpp::IntegerVector assignment(dim);
  for (int i = 0; i < dim; ++i)
    assignment[i] = rowsol[i];

  return Rcpp::List::create(
    Rcpp::Named("cost")       = lapcost,
    Rcpp::Named("assignment") = assignment
  );
}


// ====================== STABLE LAPJV (TreeDist) ======================
//inline double xhungarian_solve(const arma::mat& cost, std::vector<int>& Assignment) {
//  ScopedTimer t("  LAPJV");

//  static Rcpp::Function LAPJV("LAPJV");   // cached

//  Rcpp::List res = LAPJV(Rcpp::wrap(cost));
//  double total_cost = Rcpp::as<double>(res["score"]);

//  Rcpp::IntegerVector matching = res["matching"];
//  Assignment.resize(matching.size());
//  for (size_t i = 0; i < Assignment.size(); ++i)
//    Assignment[i] = matching[i] - 1;

//  return total_cost;
//}


//' @name hungarian_solve
//'
//' @title Wrapper that calls native_lapjv.
//'
//' @description ...
//'
//' @author Phil Davies.
//'
//' @param cost The cost matrix.
//' @param Assignment The assignment vector.
//'
//' @return The total cost.
//'
// ====================== WRAPPER THAT CALLS native_lapjv ======================
inline double hungarian_solve(const arma::mat& cost, std::vector<int>& Assignment) {

  //ScopedTimer t(" Native LAPJV ");

  // call the native_lapjv function
  Rcpp::List res = native_lapjv(Rcpp::wrap(cost));

  // Extract total cost
  double total_cost = Rcpp::as<double>(res["cost"]);

  // Extract assignment (already 0-based from native_lapjv)
  Rcpp::IntegerVector matching = res["assignment"];
  Assignment.resize(matching.size());
  for (size_t i = 0; i < Assignment.size(); ++i) {
    Assignment[i] = matching[i];
  }

  return total_cost;
}


//' @name LinearAssignmentProblem_cpp
//'
//' @title Solve the Linear Assignment Problem.
//'
//' @description Solves the linear assignment problem using the hungarian algorithm.
//'
//' @author Phil Davies.
//'
//' @param W_R A weight matrix.
//' @param target_n The number of sample points required.
//' @param N The number of points in the population.
//' @param initial_J ...
//' @param n_threads The maximum number of threads to use. Default value is 1.
//' @param verbose When set true will display debug and informational messages.
//'
//' @return The sample matrix.
//'
//' @examples
//' # Distance matrix for the two 3D points.
//' spbalDAS::arma_dist_al(X = rbind(c(1, 2, 3), c(4, 3, 2)))
//' # Distance matrix for three 2D points.
//' spbalDAS::arma_dist_al(X = matrix(c(2, 3, 0, 9, 4, 5), nrow = 3, byrow = TRUE))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix LinearAssignmentProblem_cpp(Rcpp::NumericMatrix W_R,
                                                arma::uword target_n,
                                                arma::uword N,
                                                arma::uword initial_J,
                                                int n_threads = 1,
                                                bool verbose = false) {

  auto start_time = std::chrono::high_resolution_clock::now();   // ← for summary

  if(verbose){
    ScopedTimer total_timer("LinearAssignmentProblem_cpp TOTAL");
  }

  const arma::mat W = Rcpp::as<arma::mat>(W_R);

  std::vector<arma::uword> J(target_n);
  J[0] = initial_J;
  for (arma::uword k = 1; k < target_n; ++k)
    J[k] = std::max<arma::uword>(1u, std::min(J[k-1], N / (k + 1)));

  // Initial sample
  std::vector<std::vector<int>> SampleMatrix(initial_J, std::vector<int>(1));
  arma::uvec init_perm = arma::randperm(N, initial_J);
  for (arma::uword r = 0; r < initial_J; ++r)
    SampleMatrix[r][0] = static_cast<int>(init_perm[r] + 1);

  for (arma::uword i = 1; i < target_n; ++i) {
    if(verbose){
      ScopedTimer step_timer("Step i=" + std::to_string(i+1) + " (J=" + std::to_string(J[i]) + ")");
    }

    arma::uword curr_J = J[i];
    arma::uword prev_cols = i;

    // row selection + newSample + setdiff + Ai
    arma::uvec row_sel = arma::randperm(SampleMatrix.size()).subvec(0, curr_J - 1);
    std::vector<std::vector<int>> newSample(curr_J, std::vector<int>(i + 1));

    for (arma::uword c = 0; c < prev_cols; ++c)
      for (arma::uword r = 0; r < curr_J; ++r)
        newSample[r][c] = SampleMatrix[row_sel[r]][c];

    std::unordered_set<int> used;
    for (const auto& row : newSample)
      for (arma::uword c = 0; c < prev_cols; ++c)
        used.insert(row[c]);

    std::vector<int> OmegaNotCi;
    OmegaNotCi.reserve(N - used.size());
    for (int k = 1; k <= static_cast<int>(N); ++k)
      if (used.find(k) == used.end()) OmegaNotCi.push_back(k);

      arma::uvec r2 = arma::randperm(OmegaNotCi.size());
      for (arma::uword r = 0; r < curr_J; ++r)
        newSample[r][i] = OmegaNotCi[r2[r]];

      SampleMatrix = std::move(newSample);

      arma::ivec A_all(curr_J);
      // Cost matrix
      arma::mat C(curr_J, curr_J);
      {
        if(verbose){
          ScopedTimer t("  Build cost matrix C");
        }

        //arma::ivec A_all(curr_J);
        for (arma::uword r = 0; r < curr_J; ++r) A_all(r) = SampleMatrix[r][i];
        arma::uvec A0 = arma::conv_to<arma::uvec>::from(A_all - 1);

        if (curr_J >= 256 && n_threads > 1) {
          RcppThread::parallelFor(0, curr_J, [&](arma::uword r) {
            arma::uvec C0(prev_cols);
            for (arma::uword c = 0; c < prev_cols; ++c)
              C0(c) = SampleMatrix[r][c] - 1;

            double WCij = (prev_cols == 0) ? 0.0 : arma::accu(W.submat(C0, C0));
            arma::rowvec contrib = 2.0 * arma::sum(W(A0, C0), 1).t();
            C.row(r) = WCij + contrib;
          }, n_threads);
        } else {
          for (arma::uword r = 0; r < curr_J; ++r) {
            arma::uvec C0(prev_cols);
            for (arma::uword c = 0; c < prev_cols; ++c)
              C0(c) = SampleMatrix[r][c] - 1;

            double WCij = (prev_cols == 0) ? 0.0 : arma::accu(W.submat(C0, C0));
            arma::rowvec contrib = 2.0 * arma::sum(W(A0, C0), 1).t();
            C.row(r) = WCij + contrib;
          }
        }
      }

      // LAPJV
      {
        std::vector<int> Assignment;
        hungarian_solve(C, Assignment);

        for (arma::uword r = 0; r < curr_J; ++r)
          //SampleMatrix[r][i] = OmegaNotCi[Assignment[r]];
          SampleMatrix[r][i] = A_all[Assignment[r]];
      }
  }

  // ====================== SUMMARY LINE ======================
  auto end_time = std::chrono::high_resolution_clock::now();
  auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
  Rcpp::Rcout << "DAS completed in " << total_ms << " ms\n";

  // Return final matrix...
  arma::uword final_J = J.back();
  Rcpp::IntegerMatrix result(final_J, target_n);
  for (arma::uword r = 0; r < final_J; ++r)
    for (arma::uword c = 0; c < target_n; ++c)
      result(r, c) = SampleMatrix[r][c];

  return result;
}
