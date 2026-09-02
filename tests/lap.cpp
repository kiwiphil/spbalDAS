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
// Jonker-Volgenant LAPJV (square only - stable)
// ------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List native_lapjv(Rcpp::NumericMatrix Cost_R) {
  arma::mat assigncost = Rcpp::as<arma::mat>(Cost_R);
  int dim = assigncost.n_rows;

  if (dim == 0 || dim != assigncost.n_cols) {
    Rcpp::stop("native_lapjv requires a square cost matrix (rows == cols).\n"
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
