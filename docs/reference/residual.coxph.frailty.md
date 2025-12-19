# Residuals for shared frailty Cox proportional hazards models

Compute several types of residuals (censored Z-residuals, Cox–Snell,
martingale, and deviance) for shared frailty Cox proportional hazards
models fitted with
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) using a
multiplicative frailty term (e.g., `frailty(group)`).

## Usage

``` r
residual.coxph.frailty(
  fit_coxph,
  traindata,
  newdata,
  residual.type = c("censored Z-residual", "Cox-Snell", "martingale", "deviance")
)
```

## Arguments

- fit_coxph:

  A fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) object
  for a shared frailty Cox model, typically specified with a term such
  as `frailty(group)` in the model formula. The object must contain
  cluster-level frailty estimates in `fit_coxph$frail`.

- traindata:

  A `data.frame` used to reconstruct the baseline hazard and frailty
  effects. It must contain the survival response, all fixed-effect
  covariates, and the shared frailty grouping factor appearing in
  `fit_coxph$formula`.

- newdata:

  A `data.frame` for which residuals are to be computed. It must contain
  the survival response, covariates, and the shared frailty grouping
  factor used in the original model. The factor levels of the grouping
  variable must be compatible with those in `traindata`.

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

- `Survival.Prob`: vector of survival probabilities \\S\_{ij}(t_i)\\ for
  each observation in `newdata`.

- `linear.pred`: vector of fixed-effect linear predictors \\\eta\_{ij}\\
  (excluding the frailty term).

- `covariates`: model matrix of fixed-effect covariates used in the
  linear predictor.

- `censored.status`: event indicator (1 = event, 0 = censored).

- `object.model.frame`: the `model.frame` constructed from
  `fit_coxph$formula` and `newdata`.

## Details

The function is designed for an out-of-sample / cross-validation
setting:

- `traindata` is used to reconstruct the baseline cumulative hazard and
  the relationship between covariates and shared frailty.

- `newdata` is the dataset for which residuals are computed.

Both datasets must contain the survival response, all fixed-effect
covariates, and the same frailty grouping factor (with compatible factor
levels). The grouping variable is extracted from the
[`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html) term in
`fit_coxph$formula` and must be a factor in both `traindata` and
`newdata`. The internal implementation treats models with many groups
(`gpnumber > 5`) and few groups (`gpnumber <= 5`) slightly differently,
reflecting how frailty coefficients are stored in the fitted `coxph`
object.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  data(lung)
  lung$inst <- factor(lung$inst)

  ## Shared frailty Cox model
  set.seed(1)
  idx <- sample(seq_len(nrow(lung)), size = floor(0.7 * nrow(lung)))
  train_dat <- lung[idx, ]
  test_dat  <- lung[-idx, ]

  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
                     data = train_dat)

  ## Censored Z-residuals on the test set
  r_z <- residual.coxph.frailty(fit_frail,
                                traindata = train_dat,
                                newdata   = test_dat,
                                residual.type = "censored Z-residual")

  ## Cox–Snell residuals
  r_cs <- residual.coxph.frailty(fit_frail,
                                 traindata = train_dat,
                                 newdata   = test_dat,
                                 residual.type = "Cox-Snell")
} # }
```
