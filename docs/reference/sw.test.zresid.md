# Shapiro-Wilk Normality Test for Z-Residuals

Performs the Shapiro-Wilk test for normality on each column of a matrix
of Z-residuals.

## Usage

``` r
# S3 method for class 'test.zresid'
sw(Zresidual, ...)
```

## Arguments

- Zresidual:

  A numeric matrix of Z-residuals, where each column represents a
  separate set of residuals (e.g., from different posterior predictive
  draws or variables).

- ...:

  Additional arguments.

## Value

A numeric vector of Shapiro-Wilk p-values, one for each column of
`Zresidual`.

## Details

Infinite or non-finite values are handled by replacing negative infinity
with -1e10 and positive infinity with 1e10. Any NaN or remaining
infinite values are reported via messages.

## See also

[`shapiro.test`](https://rdrr.io/r/stats/shapiro.test.html),
[`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)

## Examples

``` r
if (FALSE) { # \dontrun{
Zres <- matrix(rnorm(100), ncol=5)
sw.pvals <- sw.test.zresid(Zres)
} # }
```
