# Cross-validated Z-residual diagnostics (generic)

Generic function for cross-validated Z-residual diagnostics. Method
dispatch is based on the class of `object`. Methods are currently
provided for survival models such as `coxph` and `survreg`.

## Usage

``` r
CV.Zresidual(object, nfolds, foldlist = NULL, data = NULL, nrep = 1, ...)
```

## Arguments

- object:

  A fitted model object.

- nfolds:

  Integer. Number of folds for cross-validation.

- foldlist:

  Optional list specifying custom fold assignments. If `NULL`, folds are
  generated internally by the method.

- data:

  Optional data frame used to refit the model during cross-validation,
  when required by the method.

- nrep:

  Integer. Number of repeated cross-validations to perform. Default is
  1.

- ...:

  Further arguments passed on to specific methods.

## Value

An object whose structure depends on the underlying method, typically
tagged with class `"cvzresid"` in addition to method-specific classes.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
  out <- CV.Zresidual(fit, nfolds = 5, data = lung)
} # }
```
