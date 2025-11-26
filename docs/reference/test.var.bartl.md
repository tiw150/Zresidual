# Bartlett Test for Homogeneity of Variances of Z-Residual

Performs Bartlett's test to assess whether the variance of Z-residuals
differs across levels of a covariate. Numeric covariates with many
distinct values are binned, and small or empty bins are removed before
testing.

## Usage

``` r
test.var.bartl(Zresidual, fitted.value, k.bl = 10)
```

## Arguments

- Zresidual:

  Numeric vector of Z-residuals.

- fitted.value:

  Numeric or factor covariate to test against.

- k.bl:

  Integer; the number of bins to discretize a numeric covariate (default
  10).

## Value

Numeric p-value from Bartlett's test for homogeneity of variances.

## Details

The function handles covariates as follows:

- If `fitted.value` is a factor or has fewer than `k.bl` unique values,
  it is treated as categorical.

- Otherwise, numeric covariates are binned into `k.bl` bins.

- Bins with fewer than 3 observations are removed.

- If insufficient bins remain, the covariate is log-transformed and
  binned again.

Bartlett's test is then applied to the Z-residuals grouped by the factor
or binned covariate.

## See also

[`bartlett.test`](https://rdrr.io/r/stats/bartlett.test.html), `split`

## Examples

``` r
if (FALSE) { # \dontrun{
Zres <- rnorm(100)
x <- runif(100)
test.var.bartl(Zres, x, k.bl = 5)
} # }
```
