# Holdout predictive checks

Computes holdout predictive checks on a new dataset using the same
thin-waist predictive-check infrastructure as `ppc`.

## Usage

``` r
hpc(
  fit,
  newdata,
  log_pointpred = NULL,
  point_Zcov = NULL,
  postpred_simulate = NULL,
  test = NULL,
  x = NULL,
  ndraws = NULL,
  seed = NULL,
  type = NULL,
  randomized = TRUE,
  k_anova = 10,
  k_bl = 10,
  ...
)
```

## Arguments

- fit:

  A fitted model object.

- newdata:

  A data frame used for holdout predictive checking.

- log_pointpred:

  Optional pointwise predictive bridge.

- point_Zcov:

  Optional point-level covariate/moment bridge.

- postpred_simulate:

  Optional posterior predictive simulation bridge.

- test:

  Optional list of Z-residual test specifications.

- x:

  Optional covariate for the default ANOVA-style diagnostic.

- ndraws:

  Optional number of posterior draws to use.

- seed:

  Optional random seed.

- type:

  Optional component selector.

- randomized:

  Logical; if `TRUE`, use randomized residuals for discrete outcomes.

- k_anova:

  Maximum number of bins used when an ANOVA covariate is numeric.

- k_bl:

  Maximum number of bins used when a Bartlett covariate is numeric.

- ...:

  Additional arguments passed to bridge functions.

## Value

An object of class `"predcheck"`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- hpc(fit = fit_nb, newdata = dat_hold, x = "depth", ndraws = 500, seed = 1)
print(res)
} # }
```
