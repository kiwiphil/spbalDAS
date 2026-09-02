// [[Rcpp::depends(RcppArmadillo, RcppHungarian)]]
#include <RcppArmadillo.h>
#include <RcppHungarian.h>

// ====================== ULTRA-FAST LAPJV (Jonker-Volgenant from TreeDist) ======================
inline double qqhungarian_solve(const arma::mat& cost, std::vector<int>& Assignment) {
  //ScopedTimer t("  LAPJV solve + assignment");   // ← will now show the real fast time

  // Call TreeDist::LAPJV directly (returns list with $cost and $assignment)
  Rcpp::Function LAPJV("LAPJV");                 // looks up once per call – negligible
  Rcpp::List res = LAPJV(Rcpp::wrap(cost));

  double total_cost = Rcpp::as<double>(res["cost"]);

  // assignment is 1-based → convert to 0-based for your existing code
  Rcpp::IntegerVector assign_r = res["assignment"];
  Assignment.resize(assign_r.size());
  for (size_t i = 0; i < Assignment.size(); ++i) {
    Assignment[i] = assign_r[i] - 1;
  }

  return total_cost;
}

// ====================== FAST DIRECT HUNGARIAN (arma::mat → DistMatrix) ======================
inline double zzhungarian_solve(const arma::mat& cost, std::vector<int>& Assignment) {
  const size_t nr = cost.n_rows;
  const size_t nc = cost.n_cols;

  Assignment.assign(nr, -1);
  if (nr == 0 || nc == 0) return 0.0;

  // Pre-reserve everything — minimal allocations
  std::vector<std::vector<double>> DistMatrix;
  DistMatrix.reserve(nr);

  const double* src = cost.memptr();   // Armadillo column-major
  const size_t stride = nr;

  for (size_t i = 0; i < nr; ++i) {
    DistMatrix.emplace_back(nc);               // exact size, no reallocation
    double* dst = DistMatrix.back().data();
    const double* row_start = src + i;
    for (size_t j = 0; j < nc; ++j) {
      dst[j] = row_start[j * stride];        // strided read (cache-friendly)
    }
  }

  HungarianAlgorithm solver;
  return solver.Solve(DistMatrix, Assignment);
}

// ===============================================
// Fast direct Hungarian solver for arma::mat
// (recommended for large matrices > ~300×300)
// ===============================================
inline double xhungarian_solve(const arma::mat& cost,
                              std::vector<int>& Assignment)
{
  const arma::uword nr = cost.n_rows;
  const arma::uword nc = cost.n_cols;

  // Pre-size assignment (HungarianAlgorithm expects it sized to #rows)
  Assignment.assign(nr, -1);

  if (nr == 0 || nc == 0) {
    return 0.0;
  }

  // === Fastest practical conversion (row-wise, cache-friendly) ===
  std::vector<std::vector<double>> DistMatrix(nr, std::vector<double>(nc));

  //for (arma::uword i = 0; i < nr; ++i) {
  //  // Armadillo row(i) returns a strided view → std::copy is heavily optimised
  //  std::copy(cost.row(i).cbegin(), cost.row(i).cend(), DistMatrix[i].begin());
  //}

  // inside hungarian_solve, replace the loop with:
  const double* src = cost.memptr();   // column-major
  const size_t stride = nr;

  for (size_t i = 0; i < nr; ++i) {
    double* dst = DistMatrix[i].data();
    const double* row_start = src + i;
    for (size_t j = 0; j < nc; ++j) {
      dst[j] = row_start[j * stride];
    }
  }

  HungarianAlgorithm solver;
  return solver.Solve(DistMatrix, Assignment);
}



