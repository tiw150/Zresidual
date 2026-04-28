# ANOVA test diagnostic test for Z-residuals

Computes an ANOVA-style diagnostic p-value for each column of a
Z-residual matrix against an index, linear predictor, covariate, or
user-supplied grouping variable.

## Usage

``` r
aov.test.zresid(
  Zresidual,
  X = c("lp", "covariate"),
  k.anova = 10,
  zcov = NULL,
  ...
)
```

## Arguments

- Zresidual:

  A numeric vector or matrix of Z-residuals.

- X:

  X-axis specification or grouping variable.

- k.anova:

  Maximum number of bins for numeric `X`.

- zcov:

  Optional metadata returned by
  [`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- ...:

  Reserved for forward compatibility.

## Value

A numeric vector of p-values.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n <- 30
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.2 * x))
  t_cens  <- rexp(n, rate = 0.5)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)
  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
  z <- Zresidual(fit,  data = dat, nrep = 2, seed = 1)
  info <- Zcov(fit, data = dat)
  aov.test.zresid(z, X = "lp", zcov = info)
}
#> Warning: NaNs produced
#> [1] 0.6839322 0.5311146
```
