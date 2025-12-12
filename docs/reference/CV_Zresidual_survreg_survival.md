# Cross-validated Z-residuals for parametric survival regression models (internal worker)

\`CV_Zresidual_survreg_survival()\` is the internal workhorse used by
\[CV.Zresidual.survreg()\] to compute cross-validated (CV) Z-residuals
for parametric survival regression models fitted by
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

It performs K-fold cross-validation, refits the `survreg` model on each
training subset, and computes Z-residuals on the held-out fold via
`Zresidual_survreg_survival()`.

## Usage

``` r
CV_Zresidual_survreg_survival(fit.survreg, data, nfolds, foldlist, n.rep, ...)
```

## Arguments

- fit.survreg:

  A fitted parametric survival regression model from survival, created
  using [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

- data:

  Optional `data.frame` used for cross-validation. It must contain the
  survival response and all covariates appearing in `fit.survreg$terms`.
  If `NULL`, the model frame of `fit.survreg` (via
  [`model.frame.survreg()`](https://tiw150.github.io/Zresidual/reference/model.frame.survreg.md))
  is used internally.

- nfolds:

  Integer. Number of cross-validation folds. If `NULL`, the number of
  folds is chosen heuristically by the caller (e.g. the S3 method
  [`CV.Zresidual.survreg()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.survreg.md)).

- foldlist:

  Optional list specifying fold indices. Each element should be an
  integer vector giving the row indices of the held-out (test) set for a
  given fold. If `NULL`, folds are created using `make_fold()`,
  stratifying by the survival response and censoring indicator.

- n.rep:

  Integer. Number of repeated Z-residual samples per observation in each
  fold (i.e. number of Monte Carlo replications for censored
  observations in `Zresidual_survreg_survival()`).

- ...:

  Further arguments passed to the internal worker function.

## Value

A numeric matrix of dimension \\n \times\\ `n.rep`, where \\n\\ is the
number of rows in `data` (if supplied) or in the internal model frame of
`fit.survreg`. Columns are typically named `"CV.Z-residual 1"`,
`"CV.Z-residual 2"`, ..., up to `n.rep`. The matrix usually carries
attributes such as:

- `"type"`: character string `"survival"`.

- `"Survival.Prob"`: vector of predicted survival probabilities at the
  observed time for each observation.

- `"linear.pred"`: vector of linear predictors on the `survreg` scale.

- `"censored.status"`: event indicator (1 = event, 0 = censored).

- `"covariates"`: data frame of covariates used for the residual
  computation.

- `"object.model.frame"`: data frame representing the model frame
  underlying the CV residuals.

## Details

The function works in two modes:

- If `data` is not `NULL`, folds are defined on the rows of `data`. For
  each fold, the model is refitted on `data[-test, ]` and Z-residuals
  are computed on `data[test, ]`.

- If `data` is `NULL`, the internal model frame of `fit.survreg` is
  used. The function reconstructs explicit time and status columns from
  the `Surv` response before refitting the `survreg` model within each
  fold.

For each fold, \`CV_Zresidual_survreg_survival()\` attempts to refit the
model using the same formula and distribution as in `fit.survreg`. If
the model fit fails (due to convergence or other errors/warnings), the
corresponding fold residuals are filled with `NA`.

The per-fold Z-residual matrices are then stacked into a single \\n
\times\\ `n.rep` matrix, where \\n\\ is the total number of
observations. Attributes such as survival probabilities, linear
predictors, censoring status, covariates, and the model frame are
reassembled in the original observation order.

## See also

[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html),
`Zresidual_survreg_survival()`,
[`CV.Zresidual.survreg()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.survreg.md),
`make_fold()`

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  data(lung)
  fit_wb <- survreg(Surv(time, status) ~ age + sex,
                    data = lung, dist = "weibull")

  cvz_wb <- CV_Zresidual_survreg_survival(
    fit.survreg = fit_wb,
    data        = lung,
    nfolds      = 5,
    foldlist    = NULL,
    n.rep       = 10
  )
} # }
```
