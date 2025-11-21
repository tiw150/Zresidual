# ANOVA Test for Z-Residuals

Performs an ANOVA test to assess whether Z-residuals differ across
levels of a covariate. Numeric covariates with many distinct values are
discretized into bins. Small or empty bins are removed before testing.

## Usage

``` r
test.nl.aov(Zresidual, fitted.value, k.anova = 10)
```

## Arguments

- Zresidual:

  Numeric vector of Z-residuals.

- fitted.value:

  Numeric or factor covariate to test against.

- k.anova:

  Integer; the maximum number of bins to discretize a numeric covariate
  (default 10).

## Value

Numeric p-value from the ANOVA F-test for the effect of the covariate on
Z-residuals.

## Details

The function handles covariates as follows:

- If `fitted.value` is a factor or has fewer than `k.anova` unique
  values, it is treated as categorical.

- Otherwise, numeric covariates are binned into `k.anova` bins using
  [`cut`](https://rdrr.io/r/base/cut.html).

- Bins with fewer than 3 observations are removed.

- If insufficient bins remain, the covariate is log-transformed and
  binned again.

ANOVA is then performed with `lm(Zresidual ~ binned_covariate)` and the
p-value for the first term is returned.

## See also

[`anova`](https://rdrr.io/r/stats/anova.html),
[`lm`](https://rdrr.io/r/stats/lm.html)

## Examples

``` r
if (FALSE) { # \dontrun{
Zres <- rnorm(100)
x <- runif(100)
test.nl.aov(Zres, x, k.anova = 5)
} # }
```
