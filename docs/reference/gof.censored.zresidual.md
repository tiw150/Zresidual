# Goodness-of-fit test for censored Z-residuals

Applies a goodness-of-fit test to censored Z-residuals using
[`EnvStats::gofTestCensored`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html).
This is typically used to assess whether censored Z-residuals are
approximately standard normal under a fitted survival model.

## Usage

``` r
gof.censored.zresidual(censored.Zresidual, test = "sf", ...)
```

## Arguments

- censored.Zresidual:

  Numeric vector or one-column matrix of censored Z-residuals. It must
  carry an attribute `"censored.status"`.

- test:

  Character string passed to
  [`EnvStats::gofTestCensored()`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html).
  The default is `"sf"`.

- ...:

  Additional arguments passed to
  [`EnvStats::gofTestCensored()`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html).

## Value

A single numeric p-value.

## See also

[`gofTestCensored`](https://alexkowa.github.io/EnvStats/reference/gofTestCensored.html)

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE) &&
    requireNamespace("EnvStats", quietly = TRUE)) {
  set.seed(1)
  n <- 30
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.3 * x))
  t_cens  <- rexp(n, rate = 0.5)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)

  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
  rz <- surv_residuals(fit, data = dat, residual.type = "censored Z-residual")
  gof.censored.zresidual(rz)
}
#> [1] 0.999449
```
