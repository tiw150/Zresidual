# Z-residual predictive checks

Computes probability-space predictive diagnostics from pointwise
predictive log-survival and log-density/probability matrices returned by
`log_pointpred`. The final PPC/HPC p-value is computed by comparing the
test statistic/p-value from posterior predictive replicated data against
the same quantity from the observed data, using the same posterior draw.

## Usage

``` r
predcheck_zresid(
  fit,
  data,
  data_role = c("training", "holdout"),
  test = NULL,
  log_pointpred = NULL,
  point_Zcov = NULL,
  postpred_simulate = NULL,
  type = NULL,
  x = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  seed = NULL,
  randomized = TRUE,
  k_anova = 10,
  k_bl = 10,
  eps = 1e-12,
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

- test:

  A list of test specifications. Supported aliases include `"sw"`,
  `"aov"`, and `"bartlett"`.

- log_pointpred:

  Optional pointwise predictive bridge.

- point_Zcov:

  Optional point-level covariate bridge.

- postpred_simulate:

  Optional posterior predictive simulation bridge.

- type:

  Optional component selector.

- x:

  Optional covariate for the default ANOVA-style test.

- ndraws:

  Optional number of posterior draws.

- draw_ids:

  Optional posterior draw indices.

- seed:

  Optional random seed.

- randomized:

  Logical; if `TRUE`, use randomized residuals for discrete outcomes. If
  `FALSE`, use midpoint residuals with `U=0.5`.

- k_anova:

  Maximum number of bins for numeric ANOVA covariates.

- k_bl:

  Maximum number of bins for numeric Bartlett covariates.

- eps:

  Probability clipping constant.

- ...:

  Additional arguments passed to bridges.

## Value

An object of class `"z_predcheck"`.
