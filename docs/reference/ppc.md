# Posterior predictive checks

Computes posterior predictive checks using the package thin-waist
bridges: `log_pointpred`, `point_Zcov`, and `postpred_simulate`.

## Usage

``` r
ppc(
  fit,
  data,
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

- data:

  A data frame used for posterior predictive checking.

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
res <- ppc(fit = fit_nb, data = dat, x = "depth", ndraws = 500, seed = 1)
print(res)
} # }
```
