# Cross-validated Z-residuals for Cox proportional hazards models

S3 method for
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md)
applied to Cox proportional hazards models fitted with
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). The
method automatically detects whether a frailty term is present in the
model formula and dispatches to internal implementations for standard
Cox models and shared frailty Cox models.

## Usage

``` r
# S3 method for class 'coxph'
CV.Zresidual(object, nfolds, foldlist = NULL, data = NULL, nrep = 1, ...)
```

## Arguments

- object:

  A fitted
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model
  object.

- nfolds:

  Integer. Number of folds for cross-validation.

- foldlist:

  Optional list specifying custom fold assignments. If `NULL`, folds are
  generated internally.

- data:

  Optional data frame used to refit the model during cross-validation.
  Required when `foldlist` is supplied or when the original model call
  does not contain the data explicitly.

- nrep:

  Integer. Number of repeated cross-validations to perform. Default is
  1.

- ...:

  Further arguments passed to the internal worker functions.

## Value

An object of class `"cvzresid"` containing cross-validated Z-residual
diagnostics for the Cox model.

## Details

Depending on the presence of a frailty term such as `frailty(group)` in
`object$terms`, the method calls:

- `CV_Zresidual_coxph_survival()` for standard Cox models;

- `CV_Zresidual_coxph_frailty_survival()` for shared frailty Cox models.

The returned object is tagged with class `"cvzresid"` in addition to any
classes returned by the internal worker.

## See also

`CV_Zresidual_coxph_survival()`,
`CV_Zresidual_coxph_frailty_survival()`, and the generic
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
  cv_out <- CV.Zresidual(fit, nfolds = 5, data = lung)
} # }
```
