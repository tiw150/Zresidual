# Z-residuals for shared frailty Cox proportional hazards models

This function calculates randomized Z-residuals for shared frailty
survival models fitted with
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) using a
multiplicative frailty term (e.g., `frailty(group)`).

## Usage

``` r
Zresidual.coxph.frailty(fit_coxph, traindata, newdata, n.rep = 1)
```

## Arguments

- fit_coxph:

  A fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) object
  for a shared frailty model, typically specified with a term such as
  `frailty(group)` in the model formula. The object must contain the
  cluster-level frailty estimates in `fit_coxph$frail`.

- traindata:

  Optional `data.frame` used to reconstruct the baseline hazard and
  frailty effects when computing out-of-sample residuals. It must
  contain the survival response, all fixed-effect covariates, and the
  shared frailty grouping factor appearing in `fit_coxph$formula`. Set
  to `NULL` (default) when computing in-sample residuals from the fitted
  object only.

- newdata:

  Optional `data.frame` for which Z-residuals are to be computed in the
  out-of-sample / cross-validation setting. It must contain the survival
  response, covariates, and the shared frailty grouping factor used in
  the original model. When `newdata` is not `NULL`, `traindata` must
  also be supplied. In the in-sample mode, `newdata` must be `NULL`.

- n.rep:

  Integer; number of independent randomized Z-residual replicates to
  generate. Defaults to `1`. Each replicate corresponds to a different
  randomization of the censored observations.

## Value

A numeric matrix of dimension \\n \times\\ `n.rep`, where \\n\\ is the
number of observations in the relevant model frame (original data in the
in-sample mode, or `newdata` in the out-of-sample mode). Each column
corresponds to one set of randomized Z-residuals. The returned matrix
has the following attributes attached:

- `Survival.Prob`: vector of survival probabilities \\S\_{ij}(t_i)\\ for
  each observation.

- `linear.pred`: vector of fixed-effect linear predictors \\\eta\_{ij}\\
  (excluding the frailty term).

- `covariates`: data frame of fixed-effect covariates (model frame
  without the survival response and grouping factor).

- `censored.status`: event indicator (1 = event, 0 = censored).

- `object.model.frame`: the `model.frame` used to compute the residuals
  (training frame in the in-sample mode, prediction frame in the
  out-of-sample mode).

- `type`: character string `"survival"` (if set in the current
  implementation).

## Details

There are two main usage modes:

- **In-sample mode**: If both `traindata` and `newdata` are `NULL`, the
  function extracts the model frame, covariates, and shared frailty
  structure directly from `fit_coxph` and computes Z-residuals for the
  original data used to fit the shared frailty model. The grouping
  variable for the frailty term is assumed to be the last column of the
  model frame and must be a factor.

- **Out-of-sample / cross-validation mode**: If both `traindata` and
  `newdata` are provided, the baseline cumulative hazard and frailty
  effects are reconstructed using `traindata`, while Z-residuals are
  computed for `newdata`. In this case, both datasets must contain the
  survival response, all fixed effects covariates, and the same shared
  frailty grouping factor (with compatible factor levels). Mixed cases
  where only one of `traindata` or `newdata` is supplied are not
  supported.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  ## Shared frailty Cox model (in-sample residuals)
  data(lung)
  lung$inst <- factor(lung$inst)
  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
                     data = lung)

  z_in <- Zresidual.coxph.frailty(fit_frail)

  ## Simple train/test split for out-of-sample residuals
  set.seed(1)
  idx <- sample(seq_len(nrow(lung)), size = floor(0.7 * nrow(lung)))
  train_dat <- lung[idx, ]
  test_dat  <- lung[-idx, ]

  fit_frail_cv <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
                        data = train_dat)

  z_out <- Zresidual.coxph.frailty(fit_frail_cv,
                                   traindata = train_dat,
                                   newdata   = test_dat,
                                   n.rep     = 10)
} # }
```
