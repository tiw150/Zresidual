# Cross-validated Z-residuals for parametric survival regression models

S3 method for
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md)
applied to parametric survival regression models fitted with
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html).

## Usage

``` r
# S3 method for class 'survreg'
CV.Zresidual(object, nfolds, foldlist = NULL, data = NULL, nrep = 1, ...)
```

## Arguments

- object:

  A fitted
  [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  model object.

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

  Further arguments passed to the internal worker function.

## Value

An object of class `"cvzresid"` containing cross-validated Z-residual
diagnostics for the parametric survival model.

## Details

This method delegates the actual cross-validation work to
[`CV_Zresidual_survreg_survival()`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md),
which performs fold construction, model refitting, and computation of
cross-validated Z-residuals.

The returned object is tagged with class `"cvzresid"` in addition to any
classes returned by the internal worker.

## See also

[`CV_Zresidual_survreg_survival()`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md)
and the generic
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  fit_weibull <- survreg(Surv(time, status) ~ age + sex,
                         data = lung, dist = "weibull")
  cv_out <- CV.Zresidual(fit_weibull, nfolds = 5, data = lung)
} # }
```
