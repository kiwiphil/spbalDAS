# spbalDAS

R/C++ implementation of **Dynamic Assignment Sampling** (Robertson,
Price and Reale 2024) for well-spread master samples over an auxiliary
space.

## Install

``` r

# from a local checkout
install.packages(".", repos = NULL, type = "source")
```

## Draw a sample

``` r

library(spbalDAS)
set.seed(511)
pop <- matrix(runif(400), ncol = 2)
samp <- DAS(pop, n = 25)          # J1 defaults to min(200, N/2)
samp[1, ]                         # one candidate sample of size 25
```

[`DAS()`](https://kiwiphil.github.io/spbalDAS/reference/DAS.md) returns
a matrix: each row is a candidate sample of size `n` (1-based population
indices). The first `k` columns of any row are a size-`k` DAS sample.

Doubly balanced samples (MATLAB `doubleDAS`) mix spatial spread with
balance on auxiliary variables:

``` r

aux <- cbind(pop, runif(nrow(pop)))
dsamp <- doubleDAS(pop, aux, n = 25, alpha = 0.5)
# alpha = 1  -> spatially balanced (same objective as DAS)
# alpha = 0  -> approximately balanced on aux
```

### Options that affect speed

| Argument | Default | Role |
|----|----|----|
| `J1` | `min(200, floor(N/2))` | Candidate samples. Each step solves a (J J) assignment ((O(J^3))) and you do this about `n` times. |
| `cache_W` | `N <= 2500` | Cached (N N) weights. For large `N` (thousands+) use `FALSE`; random access to an 800 MiB matrix is slower than computing distances. |
| `n_threads` | `1` | Parallel cost fill when (J^2 i ). `0` uses every logical CPU. |
| `verbose` | `FALSE` | Step timings (`fill=` vs `lapjv=`). |

## What 0.1.4 changed

Assignment uses Kuhn–Munkres (Hungarian). The Jonker–Volgenant port
could still hang around step 8–9 on real DAS runs even after a
tie-break.

## What 0.1.3 changed

Cost fill is cheaper, and `n_threads = 0` parallelises `J1` around 250.
Prefer `cache_W = FALSE` when `N` is larger than a few thousand.

## What 0.1.2 changed

[`doubleDAS()`](https://kiwiphil.github.io/spbalDAS/reference/doubleDAS.md)
is a MATLAB-faithful port: (SB + (1-)B), using the DAS C++ loop and
`native_lapjv` instead of `matchpairs`.

## What 0.1.1 changed

1.  Inverse distances are not built three times in R.
2.  The unused row-constant (s_j^W s_j) is dropped from the assignment
    cost.
3.  `J1` is optional and capped at (N/2).
4.  Jonker–Volgenant runs on the Armadillo buffer; `n_threads` /
    `verbose` are honoured.
5.  Assignment write-back uses the shuffled draw (correct ()).

## References

Robertson, B. L., Price, C. J. and Reale, M. (2024). Well-spread samples
with dynamic sample sizes. *Biometrics* 80(2), ujae026.
