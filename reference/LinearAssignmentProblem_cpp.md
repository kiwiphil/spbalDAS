# Grow a DAS master sample by successive linear assignments.

Implements Algorithm 1 of Robertson, Price and Reale (2024). Costs are
inverse Euclidean distances. The unused row-constant \\s_j^\top W s_j\\
term is omitted. Jonker-Volgenant is called on the Armadillo cost buffer
(no wrap through R).

## Usage

``` r
LinearAssignmentProblem_cpp(
  pop,
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

- target_n:

  The number of sample points required.

- initial_J:

  Number of candidate samples at the first iteration.

- n_threads:

  Maximum threads for the cost-matrix fill. 0 uses every logical CPU.
  Parallel fill runs when J^2 \* (current sample size) is at least 2048
  and n_threads is greater than 1.

- verbose:

  When TRUE, print step timings.

- cache_W:

  When TRUE, build one N by N inverse-distance matrix in C++. When
  FALSE, evaluate inverse distances on the fly (less RAM for large N).

## Value

Integer matrix of candidate samples (1-based population indices).

## Author

Phil Davies.
