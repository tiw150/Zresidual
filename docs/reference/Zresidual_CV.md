# Cross-validated Z-residuals

Compute cross-validated Z-residuals for supported frequentist survival
models. Each fold is evaluated by refitting the model on its training
data and calling
[`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)
on the corresponding held-out data.

## Usage

``` r
Zresidual_CV(
  object,
  data,
  nfolds = NULL,
  foldlist = NULL,
  nrep = 1L,
  randomized = TRUE,
  type = NULL,
  log_pointpred = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A fitted survival model. Supported classes are `"coxph"` and
  `"survreg"`.

- data:

  A data frame containing all variables used to fit `object`.

- nfolds:

  Number of cross-validation folds. When both `nfolds` and `foldlist`
  are `NULL`, the backend uses `min(10, nrow(complete_case_data))`.

- foldlist:

  Optional list of held-out row indices. Indices refer to the
  complete-case data used internally. Supply either `nfolds` or
  `foldlist`, not both.

- nrep:

  Number of randomized Z-residual replicates for each observation.

- randomized:

  Logical. If `TRUE`, randomized residuals are generated for censored
  observations. If `FALSE`, one non-randomized residual is returned
  regardless of `nrep`.

- type:

  Optional model subtype passed to the model-specific backend.

- log_pointpred:

  Optional predictive backend function or function name. When `NULL`,
  the backend is resolved automatically.

- seed:

  Optional single integer seed used for fold generation and residual
  randomization.

- ...:

  Additional arguments passed to the model-specific CV backend and
  ultimately to
  [`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

## Value

A matrix-like object with class
`c("cvzresid", "zresid", "matrix", "array")`. Rows correspond to
complete observations and columns correspond to residual replicates.
Attributes include fold indices, original row indices, probability-scale
residuals, covariates, linear predictors, and a `zcov` metadata list.

## Details

Currently supported fitted-model classes are `"coxph"`, including
`"coxph.penal"` frailty models, and `"survreg"`.

`Zresidual_CV()` is the public dispatcher. Cox calculations are handled
by
[`Zresidual_CV_survival_coxph()`](https://tiw150.github.io/Zresidual/reference/Zresidual_CV_survival_coxph.md),
and parametric survival calculations are handled by
[`Zresidual_CV_survival_survreg()`](https://tiw150.github.io/Zresidual/reference/Zresidual_CV_survival_survreg.md).

For frailty Cox models, the current CV target is within-group
prediction: every frailty group appearing in a test fold must also
appear in that fold's training data.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)

lung2 <- stats::na.omit(
  survival::lung[, c("time", "status", "age", "sex", "ph.ecog")]
)

fit <- survival::coxph(
  survival::Surv(time, status == 2) ~ age + sex + ph.ecog,
  data = lung2
)

zcv <- Zresidual_CV(
  object = fit,
  data = lung2,
  nfolds = 5,
  nrep = 10,
  seed = 123
)

print(zcv)
plot(zcv)
qqnorm(zcv)
} # }
```
