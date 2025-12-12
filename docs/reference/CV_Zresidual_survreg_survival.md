# Cross-validated Z-residuals for Survreg Models

Internal function to compute cross-validated Z-residuals for
\*\*parametric\*\* survival regression models fitted with
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

## Usage

``` r
CV_Zresidual_survreg_survival(fit.survreg, data, nfolds, foldlist, n.rep, ...)
```

## Arguments

- fit.survreg:

  A fitted [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html)
  model object.

- data:

  Optional `data.frame` used for cross-validation. Highly recommended if
  the original model was fit without specifying the `data` argument or
  if `foldlist` is supplied.

- nfolds:

  Integer. Number of folds for cross-validation (\$K\$ in K-fold CV).

- foldlist:

  Optional list specifying custom fold assignments. If `NULL`, folds are
  generated internally, typically stratified by the survival response.

- n.rep:

  Integer. Number of repeated Z-residual samples to generate per
  observation (Monte Carlo replications for censored observations).

- ...:

  Additional arguments passed to the residual calculation function
  `Zresidual_survreg_survival()`.

## Value

A numeric matrix containing the cross-validated Z-residuals (\\N \times
nrep\\), where \$N\$ is the total number of observations. The matrix
carries the following diagnostic attributes:

- `Survival.Prob`: Out-of-sample predicted survival probabilities.

- `linear.pred`: Out-of-sample linear predictors (on the `survreg`
  scale).

- `censored.status`: Event indicator (1 = event, 0 = censored).

- `covariates`: Data frame of covariates used.

- `object.model.frame`: The full model frame used for CV residual
  calculation.

- `type`: Character string, typically `"survival"`.

## Details

This function implements the K-fold cross-validation procedure. In each
fold, the parametric survival model (preserving the original
distribution, e.g., Weibull) is refitted on the training data. The
out-of-sample randomized Z-residuals are then calculated on the held-out
test data.

\*\*Data Handling Note:\*\* If `data` is `NULL`, the internal model
frame must be manually reconstructed into a standard data frame
(un-packing the `Surv` object) before the `survreg` model can be
successfully refitted on the training subset of each fold. Failed model
fits during cross-validation result in `NA` residuals for the
corresponding test fold.
