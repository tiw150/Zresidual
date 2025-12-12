# Cross-validated Z-residuals for Cox proportional hazards models

S3 method for the generic function
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md)
applied to Cox proportional hazards models fitted with
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). This
method computes the cross-validated Z-residuals, automatically detecting
whether a frailty term is present in the model formula. It then
dispatches to the appropriate internal implementation for standard Cox
models or shared frailty Cox models.

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

  Integer. The number of folds for cross-validation (K in K-fold CV).

- foldlist:

  Optional list specifying custom fold assignments. If `NULL`, folds are
  generated internally, often ensuring a balanced distribution of
  events.

- data:

  Optional data frame used to refit the model during cross-validation.
  It is required if the original model call did not contain the data
  explicitly, or when `foldlist` is supplied.

- nrep:

  Integer. Number of repeated cross-validations to perform. Default
  is 1. Repeating CV provides a more stable estimate of the residuals.

- ...:

  Further arguments passed to the internal worker functions.

## Value

A matrix of class `"cvzresid"` (and others inherited from the internal
worker) with dimension \\N \times nrep\\ (where \$N\$ is the number of
observations). The matrix columns contain the cross-validated
Z-residuals.

The object also includes the following diagnostic attributes:

- `Survival.Prob`: Cross-validated predicted survival probabilities.

- `linear.pred`: Cross-validated linear predictors
  (\\\mathbf{x}\hat{\mathbf{\beta}}\\).

- `censored.status`: The event indicator from the survival object.

- `covariates`: The covariates used in the model.

- `object.model.frame`: The model frame of the original data.

## Details

The method determines the correct worker function by examining
`attr(object$terms, "specials")$frailty`:

- \*\*Standard Cox Models\*\*: If no frailty term is found, it calls
  [`CV_Zresidual_coxph_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_coxph_survival.md).

- \*\*Shared Frailty Cox Models\*\*: If a frailty term (e.g.,
  `frailty(group)`) is present, it calls
  [`CV_Zresidual_coxph_frailty_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_coxph_frailty_survival.md).

The returned object is a matrix of Z-residuals with crucial diagnostic
information stored as attributes.

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`CV_Zresidual_coxph_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_coxph_survival.md),
[`CV_Zresidual_coxph_frailty_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_coxph_frailty_survival.md),
and the generic
[`CV.Zresidual`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  # Example 1: Standard Cox Model
  fit_std <- coxph(Surv(time, status) ~ age + sex, data = lung)
  cv_out_std <- CV.Zresidual(fit_std, nfolds = 5, data = lung)

  # Example 2: Shared Frailty Cox Model (assuming the lung data has a group column 'inst')
  # lung$inst_factor <- as.factor(lung$inst)
  # fit_frty <- coxph(Surv(time, status) ~ age + sex + frailty(inst_factor), data = lung)
  # cv_out_frty <- CV.Zresidual(fit_frty, nfolds = 5, data = lung)
} # }
```
