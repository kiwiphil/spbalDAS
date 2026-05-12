// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <algorithm>

// ====================== ROBUST NATIVE LAPJV (Jonker-Volgenant) ======================
struct LapJV {
  std::vector<double> u, v, d, cost;
  std::vector<int> pred, match;
  std::vector<bool> visited;
  size_t n;
  const double EPS = 1e-10;

  LapJV(size_t max_n = 1024) {
    u.resize(max_n);
    v.resize(max_n);
    d.resize(max_n);
    cost.resize(max_n * max_n);
    pred.resize(max_n);
    match.resize(max_n);
    visited.resize(max_n);
  }

  double solve(const arma::mat& CostMatrix, std::vector<int>& Assignment) {
    n = CostMatrix.n_rows;
    if (n == 0) {
      Assignment.clear();
      return 0.0;
    }

    // Copy cost (row-major)
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
        cost[i * n + j] = CostMatrix(i, j);

    u.assign(n, 0.0);
    v.assign(n, 0.0);
    match.assign(n, -1);
    Assignment.assign(n, -1);

    // Column reduction
    for (size_t j = 0; j < n; ++j) {
      double minv = cost[j];
      for (size_t i = 1; i < n; ++i)
        if (cost[i * n + j] < minv) minv = cost[i * n + j];
        v[j] = minv;
    }

    // Row reduction + initial matching
    for (size_t i = 0; i < n; ++i) {
      double minv = std::numeric_limits<double>::infinity();
      size_t minj = 0;
      for (size_t j = 0; j < n; ++j) {
        double val = cost[i * n + j] - v[j];
        if (val < minv) {
          minv = val;
          minj = j;
        }
      }
      u[i] = minv;
      if (match[minj] == -1) {
        match[minj] = i;
        Assignment[i] = minj;
      }
    }

    // Main augmentation loop
    for (size_t i = 0; i < n; ++i) {
      if (Assignment[i] != -1) continue;

      std::fill(visited.begin(), visited.end(), false);
      std::fill(d.begin(), d.end(), std::numeric_limits<double>::infinity());
      std::fill(pred.begin(), pred.end(), -1);

      size_t cur = i;
      bool found = false;

      while (!found) {
        visited[cur] = true;
        double delta = std::numeric_limits<double>::infinity();
        size_t minj = n;

        for (size_t j = 0; j < n; ++j) {
          if (visited[j]) continue;
          double val = cost[cur * n + j] - u[cur] - v[j];
          if (val < d[j] - EPS) {
            d[j] = val;
            pred[j] = cur;
          }
          if (d[j] < delta) {
            delta = d[j];
            minj = j;
          }
        }

        for (size_t j = 0; j < n; ++j) {
          if (visited[j]) {
            u[match[j]] += delta;
            v[j] -= delta;
          } else {
            d[j] -= delta;
          }
        }

        if (match[minj] == -1) {
          size_t j = minj;
          while (j != n) {
            size_t p = pred[j];
            match[j] = p;
            Assignment[p] = j;
            j = match[p];
          }
          found = true;
        } else {
          cur = match[minj];
        }
      }
    }

    double total = 0.0;
    for (size_t i = 0; i < n; ++i)
      total += cost[i * n + Assignment[i]];

    return total;
  }
};

// ====================== EXPORTED FUNCTION FOR sourceCpp ======================
// [[Rcpp::export]]
Rcpp::List native_lapjv(Rcpp::NumericMatrix Cost_R) {
  arma::mat Cost = Rcpp::as<arma::mat>(Cost_R);

  static LapJV solver(1024);          // pre-allocated
  std::vector<int> Assignment;

  double total_cost = solver.solve(Cost, Assignment);

  return Rcpp::List::create(
    Rcpp::Named("cost")       = total_cost,
    Rcpp::Named("assignment") = Assignment   // 0-based
  );
}
