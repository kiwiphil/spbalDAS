# Dynamic Assignment Sampling (DAS).

DAS draws well-spread master samples over an auxiliary space (Robertson,
Price and Reale 2024). Each row of the returned matrix is a candidate
sample of size `n`; the first `k` columns of any row are a size-`k` DAS
sample.

## Usage

``` r
DAS(pop, n, J1 = NULL, n_threads = 1, verbose = FALSE, cache_W = NULL)
```

## Arguments

- pop:

  Numeric matrix of population coordinates, one unit per row.

- n:

  The number of sample points to draw from the population.

- J1:

  Maximum number of candidate samples, `2 <= J1 <= floor(N/2)`. Default
  is `min(200, floor(N/2))`. Values above 500 are allowed but trigger a
  warning (large assignment problems dominate run time).

- n_threads:

  Threads used when filling the cost matrix. `0` uses every logical CPU.
  Parallel fill runs when \\J^2 i\\ is at least 2048 (so `J1 = 250` does
  use extra cores). Default is 1.

- verbose:

  Print step timings. Default is FALSE.

- cache_W:

  If `TRUE`, build one \\N \times N\\ inverse-distance matrix in C++. If
  `FALSE`, compute inverse distances on the fly. Default is `TRUE` when
  `N <= 2500`, otherwise `FALSE`. At large `N` the cached matrix is
  RAM-bound and usually slower.

## Value

An integer matrix with one candidate sample per row (1-based population
indices).

## Author

This function was first written by Blair Robertson in Matlab and later
re-written in R/C++ by Phil Davies.

## Examples

``` r
set.seed(511)
pop <- matrix(runif(200), ncol = 2)
sampMx <- DAS(pop = pop, n = 15)
#> [2026-09-05 04:33:19.822] Starting DAS N=100 n=15 J1=50 cache_W=TRUE n_threads=1
#> [2026-09-05 04:33:19.823] Finished DAS.
sampMx[1:3, ]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13] [,14]
#> [1,]   14   18   34   86   62   78   38   10   68    47    27    61    35    31
#> [2,]   82   93   16   67   71   49   76    3   80    48    52    17    64    13
#> [3,]   37   56   28   40   97   98   79   43   95    81     6    63    19    50
#>      [,15]
#> [1,]    96
#> [2,]    99
#> [3,]    54
```
