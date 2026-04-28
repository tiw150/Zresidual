# Shapiro-Wilk normality test for Z-residuals

Applies the Shapiro-Wilk test column-wise to a matrix of Z-residuals.

## Usage

``` r
sw.test.zresid(Zresidual, ...)
```

## Arguments

- Zresidual:

  A numeric vector or matrix of Z-residuals.

- ...:

  Optional named arguments. Supported arguments include `max_n` and
  `seed`.

## Value

A numeric vector of p-values, one per column.

## Examples

``` r
set.seed(1)
z <- matrix(rnorm(60), ncol = 2)
sw.test.zresid(z)
#> [1] 0.1702537 0.9482203
```
