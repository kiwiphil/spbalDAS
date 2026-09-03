# Solve the Linear Assignment Problem.

Solves the linear assignment problem using the Kuhn–Munkres (Hungarian)
algorithm. The previous Jonker–Volgenant port could cycle on
floating-point DAS costs.

## Usage

``` r
native_lapjv(Cost_R)
```

## Arguments

- Cost_R:

  The cost matrix.

## Value

A list containing two variables, cost and assignment.

## Author

Phil Davies.

## Examples

``` r
# Cost matrix for three 3D points. $cost = 1, $assignment = 2, 1, 0
spbalDAS::native_lapjv(Cost_R = rbind(c(1, 2, 0), c(2, 0, 1), c(1, 4, 19)))
#> $cost
#> [1] 1
#> 
#> $assignment
#> [1] 2 1 0
#> 
# Cost matrix for three 3D points. $cost = 10, $assignment = 1, 0, 2
spbalDAS::native_lapjv(Cost_R = rbind(c(4, 2, 8), c(2, 3, 7), c(3, 1, 6)))
#> $cost
#> [1] 10
#> 
#> $assignment
#> [1] 1 0 2
#> 
```
