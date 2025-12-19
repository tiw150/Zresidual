# Cross-validated Z-residuals for Shared Frailty Cox Models

Internal function to compute cross-validated Z-residuals for **shared
frailty** Cox proportional hazards models.

## Usage

``` r
CV_Zresidual_coxph_frailty_survival(
  fit.coxph,
  data,
  nfolds,
  foldlist,
  n.rep,
  ...
)
```

## Arguments

- fit.coxph:

  A fitted
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model
  object containing a frailty term.

- data:

  Optional data frame containing the data. Required if the original
  model was fit without specifying the `data` argument or if `foldlist`
  is supplied.

- nfolds:

  Integer. Number of folds for cross-validation (K in K-fold CV).

- foldlist:

  Optional list specifying custom fold assignments.

- n.rep:

  Integer. Number of repeated cross-validations.

- ...:

  Additional arguments.

## Value

A matrix containing the cross-validated Z-residuals (\\N \times nrep\\)
with diagnostic attributes: `Survival.Prob`, `linear.pred`,
`censored.status`, `covariates`, and `object.model.frame`.

## Details

This function implements the K-fold cross-validation procedure
specifically for shared frailty Cox models (i.e., those containing a
[`frailty()`](https://rdrr.io/pkg/survival/man/frailty.html) term). In
each fold, the model is refitted on the training data, and the
Z-residuals are calculated on the held-out test data.
