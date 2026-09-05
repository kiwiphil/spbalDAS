# Doubly balanced Dynamic Assignment Sampling (doubleDAS).

Port of the MATLAB `doubleDAS` sampler. Candidate samples grow by
successive linear assignment with cost \$\$C = \alpha\\SB +
(1-\alpha)\\B\$\$ where \\SB\\ is the inverse-distance spread used by
[`DAS`](https://kiwiphil.github.io/spbalDAS/reference/DAS.md) (MATLAB
`SBMoranIvec`) and \\B\\ is the squared sum of auxiliary variables after
appending a candidate unit (MATLAB `LinearAssignment`).

`alpha = 1` is spatially balanced (same objective as `DAS`). `alpha = 0`
is approximately balanced on `aux`. Values in between are doubly
balanced.

Reuses `arma_dist_al` weights, the DAS candidate-growth schedule, and
`native_lapjv` (Jonker-Volgenant) in place of MATLAB `matchpairs`.

## Usage

``` r
doubleDAS(
  pop,
  aux,
  n,
  alpha = 0.5,
  J1 = NULL,
  n_threads = 1,
  verbose = FALSE,
  cache_W = NULL
)
```

## Arguments

- pop:

  Numeric matrix of population coordinates, one unit per row.

- aux:

  Numeric matrix of auxiliary variables, one row per unit.

- n:

  Sample size.

- alpha:

  Mix in \\\[0, 1\]\\. Default 0.5.

- J1:

  Maximum number of candidate samples, `2 <= J1 <= floor(N/2)`. Default
  `min(500, floor(N/2))`, matching the MATLAB.

- n_threads:

  Threads for the cost-matrix fill. See
  [`DAS`](https://kiwiphil.github.io/spbalDAS/reference/DAS.md).

- verbose:

  Print step timings. Default FALSE.

- cache_W:

  Cache the \\N \times N\\ inverse-distance matrix when `alpha > 0`.
  Default `N <= 2500`. At large `N` this is usually slower than
  on-the-fly distances.

## Value

Integer matrix of candidate samples (1-based population indices). Each
row is a sample of size `n`; the first `k` columns of any row are a
size-`k` doubleDAS sample.

## Author

Blair Robertson (Matlab); R/C++ port by Phil Davies.

## Examples

``` r
set.seed(511)
pop <- matrix(runif(200), ncol = 2)
aux <- matrix(runif(300), ncol = 3)
samp <- doubleDAS(pop, aux, n = 12, alpha = 0.5)
#> [2026-09-05 04:33:20.408] Starting doubleDAS N=100 n=12 J1=50 alpha=0.500 cache_W=TRUE
#> [2026-09-05 04:33:20.409] Finished doubleDAS.
samp[1, ]
#>  [1]  42  43  50  98  61 100  91  54  56   8  77  30
```
