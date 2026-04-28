# Compute residual diagnostics for supported survival models

Unified wrapper for residual diagnostics from supported survival models.
Depending on the fitted model class, the function dispatches internally
to the appropriate residual computation for Cox models, shared-frailty
Cox models, or parametric survival regression models.

## Usage

``` r
surv_residuals(
  fit.object,
  data,
  residual.type = c("censored Z-residual", "Cox-Snell", "martingale", "deviance")
)
```

## Arguments

- fit.object:

  A fitted survival model object.

- data:

  A data frame containing the variables needed to evaluate residuals.

- residual.type:

  Character string specifying the residual type.

## Value

A residual object whose structure depends on `residual.type`. The result
is usually numeric with attributes used by downstream plotting
functions.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n <- 30
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.3 * x))
  t_cens  <- rexp(n, rate = 0.5)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)

  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
  r <- surv_residuals(fit, data = dat, residual.type = "martingale")
  head(as.vector(r))
}
#> [1]  0.9017249 -0.9004452  0.7137561  0.4771257  0.3631152  0.8328093
```
