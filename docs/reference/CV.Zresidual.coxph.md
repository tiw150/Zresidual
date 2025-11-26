# Cross-Validated Z-Residuals for Cox Proportional Hazards Models

Computes cross-validated (CV) Z-residuals for a Cox proportional hazards
model (\`coxph\`), optionally including user-supplied data or performing
CV using the model's internal model frame. The function performs K-fold
cross-validation, fits the Cox model on each training subset, and
computes Z-residuals on the held-out fold.

## Usage

``` r
CV.Zresidual.coxph(fit.coxph, data, nfolds, foldlist, n.rep)
```

## Arguments

- fit.coxph:

  A fitted Cox proportional hazards model from survival, created using
  \`coxph()\`.

- data:

  Optional data frame used for cross-validation. If \`NULL\`, the model
  frame of \`fit.coxph\` is used internally.

- nfolds:

  Integer. Number of cross-validation folds. If \`NULL\`, a default
  number of folds is chosen heuristically based on the number of
  available cores (typically around 10 folds).

- foldlist:

  Optional list specifying fold indices. If \`NULL\`, folds are created
  using \`make_fold()\`, stratifying by survival time and censoring
  status.

- n.rep:

  Integer. Number of repeated Z-residual samples per observation.

## Value

A numeric matrix of dimension \\n \times\\ `n.rep`, where \\n\\ is the
number of rows in `data` (if supplied) or in the internal model frame of
`fit.coxph`. Columns are named \`"CV.Z-residual 1"\`, \`"CV.Z-residual
2"\`, …, up to `n.rep`. The matrix has the following attributes:

\* \`type\`: character string \`"survival"\`. \* \`Survival.Prob\`:
numeric vector of predicted survival probabilities at the observed time
for each observation. \* \`linear.pred\`: numeric vector of Cox linear
predictors. \* \`censored.status\`: event indicator (1 = event, 0 =
censored). \* \`covariates\`: data frame of covariates used in the
residual computation. \* \`object.model.frame\`: data frame representing
the model frame underlying the CV residuals.

## Details

The function:

1\. Extracts covariates and event information from either \`data\` or
the model frame of \`fit.coxph\`. 2. Creates or uses supplied folds. 3.
For each fold: - Fits the Cox model to the training subset. - Computes
Z-residuals on the held-out subset using \`Zresidual.coxph()\`. -
Handles failed model fits by filling the fold with \`NA\`. 4. Combines
results into a matrix of dimension \*n × n.rep\*, where \*n\* is the
number of observations.

## See also

\`coxph()\`, \`Zresidual.coxph()\`, \`make_fold()\`, \`CV.Zresidual()\`

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  fit <- coxph(Surv(time, status) ~ age + sex, data = lung)

  # 5-fold CV Z-residuals
  cvz <- CV.Zresidual.coxph(fit, data = lung, nfolds = 5, n.rep = 10)
} # }
```
