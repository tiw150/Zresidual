# Predictive checks

Performs predictive checks on evaluation data. The same computational
algorithm is used for in-sample and holdout evaluation.

## Usage

``` r
predcheck(
  fit,
  data,
  data_role = c("training", "holdout"),
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

  Data used to evaluate the fitted model.

- data_role:

  Role of `data`: `"training"` if the data were used to fit the model,
  or `"holdout"` otherwise.

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

  Logical; use randomized residuals for discrete outcomes.

- k_anova:

  Maximum number of bins for numeric ANOVA covariates.

- k_bl:

  Maximum number of bins for numeric Bartlett covariates.

- ...:

  Additional arguments passed to bridge functions.

## Value

An object of class `"predcheck"`.

## Examples

``` r
if (FALSE) { # \dontrun{
training_result <- predcheck(
  fit_full, full_data, data_role = "training"
)
holdout_result <- predcheck(
  fit_training, test_data, data_role = "holdout"
)
} # }
```
