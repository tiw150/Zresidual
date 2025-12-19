# Residuals for Cox proportional hazards models

Compute several types of residuals (censored Z-residuals, Cox–Snell,
martingale, and deviance) for Cox proportional hazards models fitted
with [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html).

## Usage

``` r
residual.coxph(
  fit_coxph,
  newdata,
  residual.type = c("censored Z-residual", "Cox-Snell", "martingale", "deviance")
)
```

## Arguments

- fit_coxph:

  A fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model
  object without a shared frailty term. The response should be a
  right-censored survival object, typically `Surv(time, status)` or
  `Surv(tstart, tstop, status)`.

- newdata:

  A `data.frame` containing the variables required by
  `fit_coxph$formula`, including the
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) response and
  all covariates. Residuals are evaluated for the observations in
  `newdata`.

- residual.type:

  Character string specifying the type of residual to compute. Must be
  one of `"censored Z-residual"`, `"Cox-Snell"`, `"martingale"`, or
  `"deviance"`. The default is the full vector
  `c("censored Z-residual", "Cox-Snell", "martingale", "deviance")`, but
  in typical use a single value should be supplied.

## Value

A numeric matrix of dimension \\n \times 1\\, where \\n\\ is the number
of observations in `newdata`. The single column is named according to
`residual.type`. Several attributes are attached:

- `Survival.Prob`: vector of survival probabilities \\S_i(t_i)\\.

- `linear.pred`: vector of linear predictors \\\eta_i\\.

- `covariates`: model matrix of covariates (columns used in the linear
  predictor).

- `censored.status`: event indicator (1 = event, 0 = censored).

- `object.model.frame`: the `model.frame` constructed from
  `fit_coxph$formula` and `newdata`.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  data(lung)

  ## Cox PH model
  fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)

  ## Censored Z-residuals
  r_z <- residual.coxph(fit_cox, newdata = lung,
                        residual.type = "censored Z-residual")

  ## Cox–Snell residuals
  r_cs <- residual.coxph(fit_cox, newdata = lung,
                         residual.type = "Cox-Snell")

  ## Martingale residuals
  r_m <- residual.coxph(fit_cox, newdata = lung,
                        residual.type = "martingale")

  ## Deviance residuals
  r_d <- residual.coxph(fit_cox, newdata = lung,
                        residual.type = "deviance")
} # }
```
