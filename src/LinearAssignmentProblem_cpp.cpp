// [[Rcpp::depends(RcppArmadillo, RcppThread)]]
#include <RcppArmadillo.h>
#include <RcppThread.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <thread>
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

// pop is N x q, column-major: coordinate d of unit u is mem[u + d * N].
static inline double inv_dist_mem(const double* pop, arma::uword N, arma::uword q,
                                  arma::uword u, arma::uword v) {
  double s = 0.0;
  for (arma::uword d = 0; d < q; ++d) {
    const double t = pop[u + d * N] - pop[v + d * N];
    s += t * t;
  }
  return (s <= 0.0) ? 0.0 : 1.0 / std::sqrt(s);
}

static inline double inv_dist2(const double* x, const double* y, arma::uword u, arma::uword v) {
  const double dx = x[u] - x[v];
  const double dy = y[u] - y[v];
  const double s = dx * dx + dy * dy;
  return (s <= 0.0) ? 0.0 : 1.0 / std::sqrt(s);
}

// Inverse distances in one pass (no separate D then 1/D copy).
static arma::mat build_inverse_distance_weights(const arma::mat& X) {
  arma::vec G = arma::sum(arma::square(X), 1);
  arma::mat W = -2.0 * (X * X.t());
  W.each_col() += G;
  W.each_row() += G.t();
  const arma::uword n = W.n_rows;
  double* w = W.memptr();
  const arma::uword nn = n * n;
  for (arma::uword i = 0; i < nn; ++i) {
    const double d2 = w[i];
    w[i] = (d2 <= 0.0) ? 0.0 : 1.0 / std::sqrt(d2);
  }
  W.diag().zeros();
  return W;
}

static inline bool use_parallel(arma::uword curr_J, arma::uword prev_cols, int n_threads) {
  if (n_threads <= 1 || curr_J < 8) return false;
  // ~J^2 * prev_cols fused multiply-adds in the spread fill.
  return (curr_J * curr_J * prev_cols) >= 2048u;
}

// Kuhn–Munkres (Hungarian) O(n^3). Assignment[i] = column assigned to row i (0-based).
// Used instead of Jonker–Volgenant: the JV port can cycle on floating-point DAS
// costs (overlapping candidates produce equal reduced costs). This method
// marks one column per inner iteration and therefore always terminates.
//
// 1-based internals follow https://cp-algorithms.com/graph/hungarian-algorithm.html
static double lapjv_solve(const arma::mat& assigncost, std::vector<int>& Assignment) {
  const int n = static_cast<int>(assigncost.n_rows);
  if (n == 0 || n != static_cast<int>(assigncost.n_cols)) {
    Rcpp::stop("spbalDAS(lapjv_solve) requires a square cost matrix (rows == cols). "
               "Your matrix is %d x %d.",
               n, static_cast<int>(assigncost.n_cols));
  }
  if (n == 1) {
    Assignment.assign(1, 0);
    return assigncost(0, 0);
  }

  const double INF = 1e100;
  std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  std::vector<double> minv(n + 1);
  std::vector<char> used(n + 1);

  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::fill(minv.begin(), minv.end(), INF);
    std::fill(used.begin(), used.end(), 0);
    int guard = 0;
    do {
      if (++guard > n + 2) {
        Rcpp::stop("spbalDAS: Hungarian assignment stalled (n=%d).", n);
      }
      used[j0] = 1;
      const int i0 = p[j0];
      int j1 = 0;
      double delta = INF;
      for (int j = 1; j <= n; ++j) {
        if (used[j]) continue;
        const double cur = assigncost(i0 - 1, j - 1) - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);
    do {
      const int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }

  Assignment.resize(static_cast<std::size_t>(n));
  for (int j = 1; j <= n; ++j)
    Assignment[static_cast<std::size_t>(p[j] - 1)] = j - 1;

  double lapcost = 0.0;
  for (int i = 0; i < n; ++i)
    lapcost += assigncost(i, Assignment[static_cast<std::size_t>(i)]);
  return lapcost;
}


//' @name native_lapjv
//'
//' @title Solve the Linear Assignment Problem.
//'
//' @description Solves the linear assignment problem using the Kuhn–Munkres
//'   (Hungarian) algorithm. The previous Jonker–Volgenant port could cycle on
//'   floating-point DAS costs.
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

template <typename Fn>
static void run_fill(arma::uword curr_J, arma::uword prev_cols, int n_threads,
                     Fn fill_row) {
  if (!use_parallel(curr_J, prev_cols, n_threads)) {
    for (arma::uword r = 0; r < curr_J; ++r) fill_row(r);
    return;
  }
  const int nt = std::max(1, std::min(n_threads, static_cast<int>(curr_J)));
  std::atomic<arma::uword> next{0};
  std::vector<std::thread> pool;
  pool.reserve(static_cast<std::size_t>(nt));
  for (int t = 0; t < nt; ++t) {
    pool.emplace_back([&]() {
      for (;;) {
        const arma::uword r = next.fetch_add(1, std::memory_order_relaxed);
        if (r >= curr_J) break;
        fill_row(r);
      }
    });
  }
  for (auto& th : pool) th.join();
}

// C_rm is row-major curr_J x curr_J. sample[r * stride + c] is 0-based.
static void fill_spread_from_W(double* C_rm,
                               const arma::mat& W,
                               const int* sample,
                               arma::uword stride,
                               const int* A0,
                               arma::uword curr_J,
                               arma::uword prev_cols,
                               int n_threads) {
  auto fill_row = [&](arma::uword r) {
    double* row = C_rm + static_cast<std::size_t>(r) * curr_J;
    std::fill(row, row + curr_J, 0.0);
    const int* s = sample + static_cast<std::size_t>(r) * stride;
    for (arma::uword c = 0; c < prev_cols; ++c) {
      const double* wcol = W.colptr(static_cast<arma::uword>(s[c]));
      for (arma::uword k = 0; k < curr_J; ++k)
        row[k] += wcol[static_cast<arma::uword>(A0[k])];
    }
    for (arma::uword k = 0; k < curr_J; ++k) row[k] *= 2.0;
  };
  run_fill(curr_J, prev_cols, n_threads, fill_row);
}

static void fill_spread_onthefly(double* C_rm,
                                 const double* pop_mem,
                                 arma::uword N,
                                 arma::uword q,
                                 const double* xs,
                                 const double* ys,
                                 bool q2,
                                 const int* sample,
                                 arma::uword stride,
                                 const int* A0,
                                 arma::uword curr_J,
                                 arma::uword prev_cols,
                                 int n_threads) {
  auto fill_row = [&](arma::uword r) {
    double* row = C_rm + static_cast<std::size_t>(r) * curr_J;
    const int* s = sample + static_cast<std::size_t>(r) * stride;
    for (arma::uword k = 0; k < curr_J; ++k) {
      const arma::uword a = static_cast<arma::uword>(A0[k]);
      double acc = 0.0;
      if (q2) {
        for (arma::uword c = 0; c < prev_cols; ++c)
          acc += inv_dist2(xs, ys, static_cast<arma::uword>(s[c]), a);
      } else {
        for (arma::uword c = 0; c < prev_cols; ++c)
          acc += inv_dist_mem(pop_mem, N, q, static_cast<arma::uword>(s[c]), a);
      }
      row[k] = 2.0 * acc;
    }
  };
  run_fill(curr_J, prev_cols, n_threads, fill_row);
}

// MATLAB LinearAssignment.m:
//   b = sum(X(Ci(i,:), ii));
//   B(i,:) = B(i,:) + (b + X(Ai,ii)').^2;
static void fill_balance_cost(double* B_rm,
                              const arma::mat& X,
                              const int* sample,
                              arma::uword stride,
                              const int* A0,
                              arma::uword curr_J,
                              arma::uword prev_cols,
                              int n_threads) {
  const arma::uword q = X.n_cols;
  const arma::uword N = X.n_rows;
  const double* xmem = X.memptr();
  auto fill_row = [&](arma::uword r) {
    std::vector<double> b(q, 0.0);
    const int* s = sample + static_cast<std::size_t>(r) * stride;
    for (arma::uword c = 0; c < prev_cols; ++c) {
      const arma::uword u = static_cast<arma::uword>(s[c]);
      for (arma::uword d = 0; d < q; ++d)
        b[d] += xmem[u + d * N];
    }
    double* row = B_rm + static_cast<std::size_t>(r) * curr_J;
    for (arma::uword k = 0; k < curr_J; ++k) {
      const arma::uword a = static_cast<arma::uword>(A0[k]);
      double acc = 0.0;
      for (arma::uword d = 0; d < q; ++d) {
        const double t = b[d] + xmem[a + d * N];
        acc += t * t;
      }
      row[k] = acc;
    }
  };
  run_fill(curr_J, prev_cols, n_threads, fill_row);
}

static Rcpp::IntegerMatrix run_assignment_loop(const arma::mat& pop,
                                               arma::uword target_n,
                                               arma::uword initial_J,
                                               int n_threads,
                                               bool verbose,
                                               bool cache_W,
                                               const arma::mat* aux,
                                               double alpha) {
  const arma::uword N = pop.n_rows;
  const arma::uword q = pop.n_cols;
  n_threads = resolve_n_threads(n_threads);

  if (target_n < 1 || target_n >= N) {
    Rcpp::stop("spbalDAS: n must satisfy 1 <= n < N.");
  }
  if (initial_J < 2 || initial_J > N / 2) {
    Rcpp::stop("spbalDAS: J1 must satisfy 2 <= J1 <= floor(N/2).");
  }
  if (alpha < 0.0 || alpha > 1.0) {
    Rcpp::stop("spbalDAS: alpha must be in [0, 1].");
  }
  if (aux != nullptr) {
    if (aux->n_rows != N) {
      Rcpp::stop("spbalDAS: aux must have one row per population unit.");
    }
    if (aux->n_cols < 1) {
      Rcpp::stop("spbalDAS: aux must have at least one column.");
    }
  }

  const bool use_spread = (alpha > 0.0);
  const bool use_balance = (aux != nullptr && alpha < 1.0);

  auto start_time = std::chrono::high_resolution_clock::now();

  arma::mat W;
  if (use_spread && cache_W) {
    if (verbose) {
      ScopedTimer t("Build inverse-distance weights");
      W = build_inverse_distance_weights(pop);
    } else {
      W = build_inverse_distance_weights(pop);
    }
  }

  const double* pop_mem = pop.memptr();
  const bool q2 = (q == 2);
  const double* xs = q2 ? pop.colptr(0) : nullptr;
  const double* ys = q2 ? pop.colptr(1) : nullptr;

  std::vector<arma::uword> J(target_n);
  J[0] = initial_J;
  for (arma::uword k = 1; k < target_n; ++k)
    J[k] = std::max<arma::uword>(1u, std::min(J[k - 1], N / (k + 1)));

  const arma::uword stride = target_n;
  std::vector<int> sample(initial_J * stride, 0);
  arma::uvec init_perm = arma::randperm(N, initial_J);
  for (arma::uword r = 0; r < initial_J; ++r)
    sample[r * stride] = static_cast<int>(init_perm[r]);

  arma::uword n_rows = initial_J;
  std::vector<int> sample_next;
  std::vector<char> used(N);
  std::vector<int> omega;
  std::vector<int> A0;
  std::vector<int> Assignment;
  std::vector<double> B_rm;
  arma::mat C;

  for (arma::uword i = 1; i < target_n; ++i) {
    const arma::uword curr_J = J[i];
    const arma::uword prev_cols = i;
    const int dim = static_cast<int>(curr_J);

    arma::uvec row_sel = arma::randperm(n_rows).subvec(0, curr_J - 1);
    sample_next.assign(curr_J * stride, 0);
    for (arma::uword r = 0; r < curr_J; ++r) {
      std::memcpy(&sample_next[r * stride],
                  &sample[static_cast<std::size_t>(row_sel[r]) * stride],
                  prev_cols * sizeof(int));
    }

    std::fill(used.begin(), used.end(), 0);
    for (arma::uword r = 0; r < curr_J; ++r) {
      const int* s = &sample_next[r * stride];
      for (arma::uword c = 0; c < prev_cols; ++c)
        used[static_cast<arma::uword>(s[c])] = 1;
    }
    omega.clear();
    omega.reserve(N);
    for (arma::uword u = 0; u < N; ++u)
      if (!used[u]) omega.push_back(static_cast<int>(u));

    arma::uvec r2 = arma::randperm(omega.size());
    A0.resize(curr_J);
    for (arma::uword r = 0; r < curr_J; ++r) {
      A0[r] = omega[r2[r]];
      sample_next[r * stride + i] = A0[r];
    }
    sample.swap(sample_next);
    n_rows = curr_J;

    auto t_fill0 = std::chrono::high_resolution_clock::now();
    C.set_size(curr_J, curr_J);
    if (curr_J == 1) {
      C(0, 0) = 0.0;
      Assignment.assign(1, 0);
    } else {
      std::vector<double> C_rm(static_cast<std::size_t>(dim) * dim, 0.0);
      if (use_spread) {
        if (cache_W) {
          fill_spread_from_W(C_rm.data(), W, sample.data(), stride, A0.data(),
                             curr_J, prev_cols, n_threads);
        } else {
          fill_spread_onthefly(C_rm.data(), pop_mem, N, q, xs, ys, q2,
                               sample.data(), stride, A0.data(),
                               curr_J, prev_cols, n_threads);
        }
      }
      if (use_balance) {
        B_rm.resize(static_cast<std::size_t>(dim) * dim);
        fill_balance_cost(B_rm.data(), *aux, sample.data(), stride, A0.data(),
                          curr_J, prev_cols, n_threads);
        if (!use_spread) {
          C_rm.swap(B_rm);
        } else {
          const std::size_t nn = static_cast<std::size_t>(dim) * dim;
          const double om = 1.0 - alpha;
          for (std::size_t t = 0; t < nn; ++t)
            C_rm[t] = alpha * C_rm[t] + om * B_rm[t];
        }
      }
      for (int r = 0; r < dim; ++r)
        for (int k = 0; k < dim; ++k)
          C(r, k) = C_rm[static_cast<std::size_t>(r) * dim + k];
    }
    auto t_fill1 = std::chrono::high_resolution_clock::now();
    if (curr_J > 1) {
      if (!C.is_finite()) {
        C.elem(arma::find_nonfinite(C)).fill(BIG);
      }
      lapjv_solve(C, Assignment);
    }
    auto t_jv1 = std::chrono::high_resolution_clock::now();

    for (arma::uword r = 0; r < curr_J; ++r)
      sample[r * stride + i] = A0[static_cast<arma::uword>(Assignment[r])];

    if (verbose) {
      const bool noisy = (target_n <= 80) || curr_J >= 32 || i < 8 ||
                         ((i + 1) % 25 == 0) || (i + 1 == target_n);
      if (noisy) {
        const auto fill_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
          t_fill1 - t_fill0).count();
        const auto jv_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
          t_jv1 - t_fill1).count();
        Rcpp::Rcout << "Step i=" << (i + 1)
                    << " J=" << curr_J
                    << " fill=" << fill_ms << "ms"
                    << " lapjv=" << jv_ms << "ms" << std::endl;
      }
    }
  }

  if (verbose) {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time).count();
    Rcpp::Rcout << (aux == nullptr ? "DAS" : "doubleDAS")
                << " completed in " << total_ms << " ms"
                << " (n_threads=" << n_threads << ")\n";
  }

  const arma::uword final_J = J.back();
  Rcpp::IntegerMatrix result(final_J, target_n);
  for (arma::uword r = 0; r < final_J; ++r)
    for (arma::uword c = 0; c < target_n; ++c)
      result(r, c) = sample[r * stride + c] + 1;

  return result;
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
//'   logical CPU. Parallel fill runs when J^2 * (current sample size) is at
//'   least 2048 and n_threads is greater than 1.
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
  return run_assignment_loop(pop, target_n, initial_J, n_threads, verbose,
                             cache_W, nullptr, 1.0);
}

//' @name DoubleAssignmentProblem_cpp
//'
//' @title Grow a doubly balanced DAS master sample.
//'
//' @description Port of MATLAB \code{doubleDAS}: successive assignment with
//' cost \eqn{\alpha SB + (1-\alpha) B}, where \eqn{SB} is the inverse-distance
//' spread used by DAS and \eqn{B} is the squared auxiliary-sum balance term.
//' Uses the same candidate growth, cached/on-the-fly weights, and
//' Jonker-Volgenant solver as \code{LinearAssignmentProblem_cpp}.
//'
//' @author Blair Robertson (Matlab), Phil Davies (R/C++).
//'
//' @param pop Population coordinates, N rows by q columns.
//' @param aux Auxiliary matrix, N rows.
//' @param alpha Mix of spread and balance, in \eqn{[0,1]}. \code{1} is
//'   spatially balanced, \code{0} is approximately balanced.
//' @param target_n Sample size.
//' @param initial_J Number of candidate samples at the first iteration.
//' @param n_threads Maximum threads for the cost-matrix fill.
//' @param verbose When TRUE, print step timings.
//' @param cache_W When TRUE, build one N by N inverse-distance matrix.
//'
//' @return Integer matrix of candidate samples (1-based population indices).
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix DoubleAssignmentProblem_cpp(const arma::mat& pop,
                                                const arma::mat& aux,
                                                double alpha,
                                                arma::uword target_n,
                                                arma::uword initial_J,
                                                int n_threads = 1,
                                                bool verbose = false,
                                                bool cache_W = true) {
  return run_assignment_loop(pop, target_n, initial_J, n_threads, verbose,
                             cache_W, &aux, alpha);
}
