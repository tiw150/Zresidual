# Shapiro-Francia normality test for Z-residuals

Applies the Shapiro-Francia test column-wise to a matrix of Z-residuals.

## Usage

``` r
sf_test.zresid(Zresidual)
```

## Arguments

- Zresidual:

  A numeric vector or matrix of Z-residuals.

## Value

A numeric vector of p-values, one per column.

## Examples

``` r
if (requireNamespace("nortest", quietly = TRUE)) {
  set.seed(1)
  z <- matrix(rnorm(60), ncol = 2)
  sf_test.zresid(z)
}
#> [1] 0.1577635 0.9264574
```
