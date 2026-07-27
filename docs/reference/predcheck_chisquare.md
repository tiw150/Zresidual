# Pearson chi-square predictive check

Computes a moment-based predictive check by comparing the observed
Pearson discrepancy with the same discrepancy computed from posterior
predictive replicated data.

## Usage

``` r
predcheck_chisquare(
  fit,
  data,
  data_role = c("training", "holdout"),
  point_Zcov = NULL,
  postpred_simulate = NULL,
  type = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  seed = NULL,
  variance_floor = 1e-06,
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

- point_Zcov:

  Optional point-level covariate/moment bridge.

- postpred_simulate:

  Optional posterior predictive simulation bridge.

- type:

  Optional component selector.

- ndraws:

  Optional number of posterior draws.

- draw_ids:

  Optional posterior draw indices.

- seed:

  Optional random seed.

- variance_floor:

  Minimum positive variance used for numerical stability.

- ...:

  Additional arguments passed to bridges.

## Value

An object of class `"moment_predcheck"`.
