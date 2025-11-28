# Cross-Validated Z-Residuals for Shared Frailty Cox Models

Computes cross-validated (CV) Z-residuals for a shared frailty Cox
proportional hazards model fitted with
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) using a term such
as `frailty(group)`. The function performs K-fold cross-validation,
refits the frailty Cox model on each training subset, and computes
(possibly replicated) Z-residuals on the held-out fold via
[`Zresidual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/Zresidual.coxph.frailty.md).

## Usage

``` r
CV.Zresidual.coxph.frailty(fit.coxph, data, nfolds, foldlist, n.rep)
```

## Arguments

- fit.coxph:

  A fitted shared frailty Cox proportional hazards model from survival,
  created using [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  with a frailty term (e.g.
  `Surv(time, status) ~ x1 + x2 + frailty(group)`).

- data:

  Optional `data.frame` used for cross-validation. It must contain the
  survival response, all fixed-effect covariates, and the frailty
  grouping factor used in `fit.coxph$formula`. If `NULL`, the model
  frame of `fit.coxph` (via
  [`model.frame.coxph()`](https://tiw150.github.io/Zresidual/reference/model.frame.coxph.md))
  is used internally.

- nfolds:

  Integer. Number of cross-validation folds. If `NULL`, the number of
  folds is chosen heuristically based on the number of available cores,
  approximately as `10 %/% ncores * ncores`.

- foldlist:

  Optional list specifying fold indices. Each element should be an
  integer vector giving the row indices of the held-out (test) set for a
  given fold. If `NULL`, folds are created using `make_fold()`,
  stratifying primarily by the frailty grouping factor and ensuring that
  factor-by-censoring combinations are represented in each training set.

- n.rep:

  Integer. Number of repeated Z-residual samples per observation in each
  fold (i.e. number of Monte Carlo replications for censored
  observations in `Zresidual.coxph.frailty`).

## Value

A numeric matrix of dimension \\n \times\\ `n.rep`, where \\n\\ is the
number of rows in `data` (if supplied) or in the internal model frame of
`fit.coxph`. Columns are named `"CV.Z-residual 1"`, `"CV.Z-residual 2"`,
..., up to `n.rep`. The matrix has the following attributes attached:

- `"type"`: character string `"survival"`.

- `"Survival.Prob"`: vector of predicted survival probabilities at the
  observed time for each observation, assembled from all folds.

- `"linear.pred"`: vector of fixed-effect linear predictors.

- `"censored.status"`: event indicator (1 = event, 0 = censored).

- `"covariates"`: matrix or data frame of covariates used in the
  Z-residual computation (structure depends on whether `data` is
  supplied).

- `"object.model.frame"`: data frame representing the model frame
  underlying the CV residuals (assembled from all folds).

This object is intended for downstream diagnostic and visualization
tools (e.g. QQ-plots, residual-vs-covariate plots) tailored to shared
frailty Cox models.

## Details

The function operates in two modes:

- If `data` is not `NULL`, the folds are constructed (or interpreted) in
  terms of the rows of `data`, and each fold fit uses `data[-test, ]` as
  the training set and `data[test, ]` as the test set.

- If `data` is `NULL`, the internal model frame of `fit.coxph` is used.
  In this case the function reconstructs a data frame with explicit time
  and status columns plus covariates and the frailty grouping factor,
  consistent with `fit.coxph$formula`, before refitting the model in
  each fold.

In both cases, the workflow is:

1.  Extract the model frame and identify the frailty grouping factor
    (the grouping variable in
    [`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html)).

2.  Create or use supplied folds via `foldlist`. By default,
    `make_fold()` is called with the covariates, the grouping factor,
    and the event indicator to enforce reasonable balance.

3.  For each fold:

    - Fit the shared frailty Cox model to the training subset.

    - Compute Z-residuals on the held-out subset using
      [`Zresidual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/Zresidual.coxph.frailty.md).

    - If the model fit fails, fill the corresponding fold residuals with
      `NA`.

4.  Stack the per-fold residual matrices into a single \\n \times\\
    `n.rep` matrix, where \\n\\ is the total number of observations.

5.  Collect and re-assemble attributes (survival probabilities, linear
    predictors, censoring status, covariates, and model frame) in the
    original observation order.

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
`Zresidual.coxph.frailty`, `CV.Zresidual.coxph`, `make_fold`

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  ## Example data with cluster/grouping variable
  data(lung)
  set.seed(1)
  lung$cluster_id <- factor(sample(1:10, size = nrow(lung), replace = TRUE))

  ## Shared frailty Cox model
  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(cluster_id),
                     data = lung)

  ## 5-fold CV Z-residuals with 10 Monte Carlo replications
  cvz_frail <- CV.Zresidual.coxph.frailty(
    fit.coxph = fit_frail,
    data      = lung,
    nfolds    = 5,
    foldlist  = NULL,
    n.rep     = 10
  )
} # }
```
