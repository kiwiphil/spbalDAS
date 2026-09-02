// [[Rcpp::depends(RcppArmadillo, RcppThread)]]
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <chrono>
#include <cmath>
#include <limits>
#include <thread>
#include <unordered_set>
#include <vector>

struct ScopedTimer {
  std::string label;
  std::chrono::high_resolution_clock::time_point start;
  ScopedTimer(std::string l) : label(std::move(l)),
  start(std::chrono::high_resolution_clock::now()) {}
  ~ScopedTimer() {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(
      std::chrono::high_resolution_clock::now() - start).count();
    RcppThread::Rcout << label << ": " << us << " \u00b5s\n";
  }
};

using cost = double;
using row  = int;
using col  = int;

const cost BIG = 1e30;

// Declared in arma_dist_al.cpp
arma::mat arma_dist_al(const arma::mat& X);

static inline int resolve_n_threads(int n_threads) {
  if (n_threads == 0) {
    unsigned hw = std::thread::hardware_concurrency();
    return static_cast<int>(hw == 0 ? 1 : hw);
  }
  return n_threads;
}

static inline double inv_dist(const arma::mat& pop, arma::uword u, arma::uword v) {
  double s = 0.0;
  const arma::uword q = pop.n_cols;
  for (arma::uword d = 0; d < q; ++d) {
    const double t = pop(u, d) - pop(v, d);
    s += t * t;
  }
  return (s <= 0.0) ? 0.0 : 1.0 / std::sqrt(s);
}

static arma::mat build_inverse_distance_weights(const arma::mat& pop) {
  arma::mat W = arma_dist_al(pop);
  W.transform([](double d) { return (d <= 0.0) ? 0.0 : 1.0 / d; });
  W.elem(arma::find_nonfinite(W)).zeros();
  W.diag().zeros();
  return W;
}

// Jonker-Volgenant on an Armadillo buffer. Assignment[i] = column assigned to row i (0-based).
static double lapjv_solve(const arma::mat& assigncost, std::vector<int>& Assignment) {
  const int dim = static_cast<int>(assigncost.n_rows);
  if (dim == 0 || dim != static_cast<int>(assigncost.n_cols)) {
    Rcpp::stop("spbalDAS(lapjv_solve) requires a square cost matrix (rows == cols). "
               "Your matrix is %d x %d.",
               dim, static_cast<int>(assigncost.n_cols));
  }

  std::vector<col> rowsol(dim, -1);
  std::vector<row> colsol(dim, -1);
  std::vector<cost> u(dim, 0.0);
  std::vector<cost> v(dim, 0.0);

  std::vector<row> freeunassigned(dim);
  std::vector<col> collist(dim);
  std::vector<col> matches(dim, 0);
  std::vector<cost> d(dim);
  std::vector<row> pred(dim);

  int numfree = 0;

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

  cost lapcost = 0.0;
  for (row i = 0; i < dim; ++i) {
    col j = rowsol[i];
    u[i] = assigncost(i, j) - v[j];
    lapcost += assigncost(i, j);
  }

  Assignment = std::vector<int>(rowsol.begin(), rowsol.end());
  return lapcost;
}

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
  std::vector<int> assignment;
  const double lapcost = lapjv_solve(assigncost, assignment);
  return Rcpp::List::create(
    Rcpp::Named("cost")       = lapcost,
    Rcpp::Named("assignment") = assignment
  );
}

static double hungarian_solve(const arma::mat& cost, std::vector<int>& Assignment) {
  return lapjv_solve(cost, Assignment);
}

static void fill_cost_from_W(arma::mat& C,
                             const arma::mat& W,
                             const std::vector<std::vector<int>>& SampleMatrix,
                             const arma::uvec& A0,
                             arma::uword curr_J,
                             arma::uword prev_cols,
                             int n_threads) {
  auto fill_row = [&](arma::uword r) {
    arma::uvec C0(prev_cols);
    for (arma::uword c = 0; c < prev_cols; ++c)
      C0(c) = static_cast<arma::uword>(SampleMatrix[r][c] - 1);
    // Row-constant WCij is omitted: it cannot change the JV assignment.
    C.row(r) = 2.0 * arma::sum(W(A0, C0), 1).t();
  };
  if (curr_J >= 256 && n_threads > 1) {
    RcppThread::parallelFor(0, curr_J, fill_row, n_threads);
  } else {
    for (arma::uword r = 0; r < curr_J; ++r) fill_row(r);
  }
}

static void fill_cost_onthefly(arma::mat& C,
                               const arma::mat& pop,
                               const std::vector<std::vector<int>>& SampleMatrix,
                               const arma::ivec& A_all,
                               arma::uword curr_J,
                               arma::uword prev_cols,
                               int n_threads) {
  auto fill_row = [&](arma::uword r) {
    for (arma::uword k = 0; k < curr_J; ++k) {
      const arma::uword a = static_cast<arma::uword>(A_all(k) - 1);
      double acc = 0.0;
      for (arma::uword c = 0; c < prev_cols; ++c) {
        const arma::uword u = static_cast<arma::uword>(SampleMatrix[r][c] - 1);
        acc += inv_dist(pop, u, a);
      }
      C(r, k) = 2.0 * acc;
    }
  };
  if (curr_J >= 256 && n_threads > 1) {
    RcppThread::parallelFor(0, curr_J, fill_row, n_threads);
  } else {
    for (arma::uword r = 0; r < curr_J; ++r) fill_row(r);
  }
}

//' @name LinearAssignmentProblem_cpp
//'
//' @title Grow a DAS master sample by successive linear assignments.
//'
//' @description Implements Algorithm 1 of Robertson, Price and Reale (2024).
//' Costs are inverse Euclidean distances. The unused row-constant
//' \eqn{s_j^\top W s_j} term is omitted. Jonker-Volgenant is called on the
//' Armadillo cost buffer (no wrap through R).
//'
//' @author Phil Davies.
//'
//' @param pop Population coordinates, N rows by q columns.
//' @param target_n The number of sample points required.
//' @param initial_J Number of candidate samples at the first iteration.
//' @param n_threads Maximum threads for the cost-matrix fill. 0 uses every
//'   logical CPU. Parallel fill is used only when the current J is at least 256
//'   and n_threads is greater than 1.
//' @param verbose When TRUE, print step timings.
//' @param cache_W When TRUE, build one N by N inverse-distance matrix in C++.
//'   When FALSE, evaluate inverse distances on the fly (less RAM for large N).
//'
//' @return Integer matrix of candidate samples (1-based population indices).
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix LinearAssignmentProblem_cpp(const arma::mat& pop,
                                                arma::uword target_n,
                                                arma::uword initial_J,
                                                int n_threads = 1,
                                                bool verbose = false,
                                                bool cache_W = true) {

  const arma::uword N = pop.n_rows;
  n_threads = resolve_n_threads(n_threads);

  if (target_n < 1 || target_n >= N) {
    Rcpp::stop("spbalDAS(LinearAssignmentProblem_cpp) target_n must satisfy 1 <= n < N.");
  }
  if (initial_J < 2 || initial_J > N / 2) {
    Rcpp::stop("spbalDAS(LinearAssignmentProblem_cpp) initial_J must satisfy 2 <= J1 <= floor(N/2).");
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  arma::mat W;
  if (cache_W) {
    if (verbose) {
      ScopedTimer t("Build inverse-distance weights");
      W = build_inverse_distance_weights(pop);
    } else {
      W = build_inverse_distance_weights(pop);
    }
  }

  std::vector<arma::uword> J(target_n);
  J[0] = initial_J;
  for (arma::uword k = 1; k < target_n; ++k)
    J[k] = std::max<arma::uword>(1u, std::min(J[k - 1], N / (k + 1)));

  std::vector<std::vector<int>> SampleMatrix(initial_J, std::vector<int>(1));
  arma::uvec init_perm = arma::randperm(N, initial_J);
  for (arma::uword r = 0; r < initial_J; ++r)
    SampleMatrix[r][0] = static_cast<int>(init_perm[r] + 1);

  for (arma::uword i = 1; i < target_n; ++i) {
    if (verbose) {
      ScopedTimer step_timer("Step i=" + std::to_string(i + 1) +
                             " (J=" + std::to_string(J[i]) + ")");
    }

    arma::uword curr_J = J[i];
    arma::uword prev_cols = i;

    arma::uvec row_sel = arma::randperm(SampleMatrix.size()).subvec(0, curr_J - 1);
    std::vector<std::vector<int>> newSample(curr_J, std::vector<int>(i + 1));

    for (arma::uword c = 0; c < prev_cols; ++c)
      for (arma::uword r = 0; r < curr_J; ++r)
        newSample[r][c] = SampleMatrix[row_sel[r]][c];

    std::unordered_set<int> used;
    for (const auto& grow : newSample)
      for (arma::uword c = 0; c < prev_cols; ++c)
        used.insert(grow[c]);

    std::vector<int> OmegaNotCi;
    OmegaNotCi.reserve(N - used.size());
    for (int k = 1; k <= static_cast<int>(N); ++k)
      if (used.find(k) == used.end()) OmegaNotCi.push_back(k);

    arma::uvec r2 = arma::randperm(OmegaNotCi.size());
    for (arma::uword r = 0; r < curr_J; ++r)
      newSample[r][i] = OmegaNotCi[r2[r]];

    SampleMatrix = std::move(newSample);

    arma::ivec A_all(curr_J);
    for (arma::uword r = 0; r < curr_J; ++r)
      A_all(r) = SampleMatrix[r][i];

    arma::mat C(curr_J, curr_J);
    auto build_C = [&]() {
      if (cache_W) {
        arma::uvec A0 = arma::conv_to<arma::uvec>::from(A_all - 1);
        fill_cost_from_W(C, W, SampleMatrix, A0, curr_J, prev_cols, n_threads);
      } else {
        fill_cost_onthefly(C, pop, SampleMatrix, A_all, curr_J, prev_cols, n_threads);
      }
    };
    if (verbose) {
      ScopedTimer t("  Build cost matrix C");
      build_C();
    } else {
      build_C();
    }

    std::vector<int> Assignment;
    hungarian_solve(C, Assignment);
    // Columns of C correspond to A_all (the shuffled draw), not OmegaNotCi in order.
    for (arma::uword r = 0; r < curr_J; ++r)
      SampleMatrix[r][i] = A_all(Assignment[r]);
  }

  if (verbose) {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time).count();
    Rcpp::Rcout << "DAS completed in " << total_ms << " ms\n";
  }

  arma::uword final_J = J.back();
  Rcpp::IntegerMatrix result(final_J, target_n);
  for (arma::uword r = 0; r < final_J; ++r)
    for (arma::uword c = 0; c < target_n; ++c)
      result(r, c) = SampleMatrix[r][c];

  return result;
}
