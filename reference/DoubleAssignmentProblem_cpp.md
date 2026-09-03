# Grow a doubly balanced DAS master sample.

Port of MATLAB `doubleDAS`: successive assignment with cost \\\alpha
SB + (1-\alpha) B\\, where \\SB\\ is the inverse-distance spread used by
DAS and \\B\\ is the squared auxiliary-sum balance term. Uses the same
candidate growth, cached/on-the-fly weights, and Jonker-Volgenant solver
as `LinearAssignmentProblem_cpp`.

## Usage

``` r
DoubleAssignmentProblem_cpp(
  pop,
  aux,
  alpha,
  target_n,
  initial_J,
  n_threads = 1L,
  verbose = FALSE,
  cache_W = TRUE
)
```

## Arguments

- pop:

  Population coordinates, N rows by q columns.

- aux:

  Auxiliary matrix, N rows.

- alpha:

  Mix of spread and balance, in \\\[0,1\]\\. `1` is spatially balanced,
  `0` is approximately balanced.

- target_n:

  Sample size.

- initial_J:

  Number of candidate samples at the first iteration.

- n_threads:

  Maximum threads for the cost-matrix fill.

- verbose:

  When TRUE, print step timings.

- cache_W:

  When TRUE, build one N by N inverse-distance matrix.

## Value

Integer matrix of candidate samples (1-based population indices).

## Author

Blair Robertson (Matlab), Phil Davies (R/C++).
