# Goodness-of-fit test for censored Z-residuals

Perform a normality goodness-of-fit test for censored Z-residuals using
[`gofTestCensored`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html).
This is typically used to assess whether censored Z-residuals are
approximately standard normal under the fitted survival model.

## Usage

``` r
gof.censored.zresidual(censored.Zresidual)
```

## Arguments

- censored.Zresidual:

  Numeric vector (or one-column matrix) of censored Z-residuals. It must
  carry an attribute `"censored.status"` giving the censoring indicator
  (1 = event, 0 = censored) for each observation.

## Value

A single numeric value: the p-value from
[`EnvStats::gofTestCensored`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html)
using the Shapiro–Francia test (`test = "sf"`) for a right-censored
normal distribution (`distribution = "norm"`). Larger p-values indicate
no strong evidence against normality of the censored Z-residuals.

## Details

Infinite residual values (both positive and negative) are truncated to
large finite values (`+/- 1e10`) before the test is applied. The
censoring indicator is taken from the `"censored.status"` attribute of
`censored.Zresidual`.

## See also

[`gofTestCensored`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html)
