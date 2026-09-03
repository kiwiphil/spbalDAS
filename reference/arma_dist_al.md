# Calculate a distance matrix.

Computes a distance matrix for the supplied matrix of n-dimensional
points using Euclidean distance.

## Usage

``` r
arma_dist_al(X)
```

## Arguments

- X:

  An n-dimensional matrix of points.

## Value

The distance matrix for the point matrix X.

## Author

Phil Davies.

## Examples

``` r
# Distance matrix for the two 3D points.
spbalDAS::arma_dist_al(X = rbind(c(1, 2, 3), c(4, 3, 2)))
#>          [,1]     [,2]
#> [1,] 0.000000 3.316625
#> [2,] 3.316625 0.000000
# Distance matrix for three 2D points.
spbalDAS::arma_dist_al(X = matrix(c(2, 3, 0, 9, 4, 5), nrow = 3, byrow = TRUE))
#>          [,1]     [,2]     [,3]
#> [1,] 0.000000 6.324555 2.828427
#> [2,] 6.324555 0.000000 5.656854
#> [3,] 2.828427 5.656854 0.000000
```
