# Cross-validated Z-residuals for parametric survival regression models

S3 method for
[`CV.Zresidual()`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md)
applied to parametric survival regression models fitted with
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html). This method
performs K-fold cross-validation to obtain external Z-residuals for
model diagnostics.

## Usage

``` r
# S3 method for class 'survreg'
CV.Zresidual(object, nfolds, foldlist = NULL, data = NULL, nrep = 1, ...)
```

## Arguments

- object:

  A fitted [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  model object.

- nfolds:

  Integer. Number of folds for cross-validation.

- foldlist:

  Optional list specifying custom fold assignments. If `NULL`, folds are
  generated internally, typically stratified by the survival response
  and censoring indicator.

- data:

  Optional data frame used to refit the model during cross-validation.
  It is \*\*highly recommended\*\* to supply the original data here to
  ensure correct model refitting in each fold, especially when the
  original call was complex.

- nrep:

  Integer. Number of repeated Z-residual samples per observation to
  generate. Defaults to `1`. Each replicate involves re-randomizing the
  imputed survival probability for censored observations.

- ...:

  Further arguments passed to the internal worker function
  [`CV_Zresidual_survreg_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md).

## Value

An object of class `"cvzresid"` containing cross-validated Z-residual
diagnostics for the parametric survival model. It is a numeric matrix
with \$N\$ rows and `nrep` columns, accompanied by diagnostic attributes
(see
[`CV_Zresidual_survreg_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md)
for details).

## Details

This method delegates the actual cross-validation work to
[`CV_Zresidual_survreg_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md),
which handles the iterative refitting of the `survreg` model on \$K-1\$
folds and computes the randomized Z-residuals on the held-out fold.

The randomized Z-residual, \\Z\_{ij}\\, for the \$j\$-th observation in
the \$i\$-th fold is defined as: \$\$Z\_{ij} = ...\$\$ is computed based
on the predicted out-of-sample survival probability
\\\hat{S}\_{\text{train}\_i}(t_j)\\.

The returned object is tagged with class `"cvzresid"` in addition to any
classes returned by the internal worker.

## See also

[`CV_Zresidual_survreg_survival`](https://tiw150.github.io/Zresidual/reference/CV_Zresidual_survreg_survival.md),
the generic
[`CV.Zresidual`](https://tiw150.github.io/Zresidual/reference/CV.Zresidual.md),
and the \`survival\` fitting function
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  # Fit a Weibull model
  fit_weibull <- survreg(Surv(time, status) ~ age + sex,
                         data = lung, dist = "weibull")
  # Compute 5-fold cross-validated Z-residuals
  cv_out <- CV.Zresidual(fit_weibull, nfolds = 5, data = lung, nrep = 10)

  # Check the first few cross-validated residuals
  head(cv_out)
} # }
```
