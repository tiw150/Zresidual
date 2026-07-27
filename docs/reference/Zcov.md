# Extract aligned covariates and optional point-level predictive quantities

`Zcov()` extracts model-aligned response information, covariates,
response type metadata, and linear predictors for Z-residual
diagnostics. When `point_details` is supplied, it also returns
point-level quantities such as mean, variance, and lp matrices for
predictive checks.

## Usage

``` r
Zcov(
  fit,
  data,
  type = NULL,
  point_details = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)

# S3 method for class 'brmsfit'
Zcov(
  fit,
  data,
  type = NULL,
  point_details = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)

# S3 method for class 'glm'
Zcov(
  fit,
  data,
  type = NULL,
  point_details = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  ...
)

# S3 method for class 'coxph'
Zcov(fit, data, type = NULL, point_details = NULL, ndraws = NULL, ...)

# S3 method for class 'survreg'
Zcov(fit, data, type = NULL, point_details = NULL, ndraws = NULL, ...)
```

## Arguments

- fit:

  A fitted model object.

- data:

  Data used to extract aligned metadata. Must be provided for most model
  classes.

- type:

  Optional component selector. For hurdle models, common values are
  `"hurdle"`, `"count"`, and `"zero"`.

- point_details:

  Optional character vector specifying point-level quantities to return.
  Supported values are `"mean"`, `"var"`, `"variance"`, `"lp"`, and
  `"covariate"`. The point-level matrices follow the same orientation as
  `log_pointpred`: rows are draws/replications and columns are
  observations. For frequentist models, a single-row matrix is returned.

- ndraws:

  Optional number of posterior draws for Bayesian models.

- draw_ids:

  Optional posterior draw indices for Bayesian models.

- ...:

  Additional arguments passed to class-specific methods.

## Value

A named list containing the original Z-residual metadata and,
optionally, point-level predictive quantities.

## Methods (by class)

- `Zcov(brmsfit)`: Method for `brmsfit` objects.

- `Zcov(glm)`: Method for `glm` objects.

- `Zcov(coxph)`: Method for `coxph` objects from the survival package.

- `Zcov(survreg)`: Method for `survreg` objects from the survival
  package.
