# Residuals for Accelerated failure time model models

Compute several types of residuals (censored Z-residuals, Cox–Snell,
martingale, and deviance) for accelerated failure time models fitted
with [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

## Usage

``` r
residual.survreg(
  survreg_fit,
  newdata,
  residual.type = c("censored Z-residual", "Cox-Snell", "martingale", "deviance")
)
```

## Arguments

- survreg_fit:

  A fitted [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  model object. The model should be a right-censored survival regression
  model (i.e., with a `Surv(time, status)` response), using one of the
  supported distributions.

- newdata:

  A `data.frame` containing the variables required by
  `survreg_fit$terms`, including the
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response and
  all covariates. Residuals are evaluated for the observations in
  `newdata`.

- residual.type:

  Character string specifying the type of residual to compute. Must be
  one of `"censored Z-residual"`, `"Cox-Snell"`, `"martingale"`, or
  `"deviance"`. The default is the full vector
  `c("censored Z-residual", "Cox-Snell", "martingale", "deviance")`, but
  in practice a single value should be supplied.

## Value

A numeric matrix of dimension \\n \times 1\\, where \\n\\ is the number
of observations in `newdata`. The single column is named according to
`residual.type`. Several attributes are attached:

- `Survival.Prob`: vector of survival probabilities \\S_i(t_i)\\.

- `linear.pred`: vector of linear predictors \\\eta_i\\.

- `covariates`: model matrix of covariates used in the linear predictor.

- `censored.status`: event indicator (1 = event, 0 = censored).

- `object.model.frame`: the `model.frame` constructed from
  `survreg_fit$terms` and `newdata`.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  data(lung)

  ## Weibull survival regression
  fit_weib <- survreg(Surv(time, status) ~ age + sex,
                      data = lung, dist = "weibull")

  ## Censored Z-residuals
  r_z <- residual.survreg(fit_weib, newdata = lung,
                          residual.type = "censored Z-residual")

  ## Cox–Snell residuals
  r_cs <- residual.survreg(fit_weib, newdata = lung,
                           residual.type = "Cox-Snell")

  ## Martingale residuals
  r_m <- residual.survreg(fit_weib, newdata = lung,
                          residual.type = "martingale")

  ## Deviance residuals
  r_d <- residual.survreg(fit_weib, newdata = lung,
                          residual.type = "deviance")
} # }
```
