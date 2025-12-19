# Compute Z-residuals for a Bernoulli/Logistic brms model

Computes Z-residuals for a fitted Bayesian Bernoulli/Logistic (binary)
model fitted with
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) and
`family = bernoulli()`. Z-residuals are calculated using posterior
predictive methods and can be used for model diagnostics.

This is an internal workhorse for
[`Zresidual.bernoulli.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.bernoulli.brms.md)
and is not intended to be called directly by end users.

## Usage

``` r
Zresidual_bernoulli_brms(
  fit,
  method = "iscv",
  n.rep = 1,
  data = NULL,
  type = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted `brmsfit` model object with Bernoulli family.

- method:

  Character string specifying the residual calculation method: `"iscv"`
  for importance-sampled cross-validated randomized predictive p-values,
  `"rpost"` for randomized posterior predictive p-values, or `"mpost"`
  for middle-value posterior predictive p-values. Default is `"iscv"`.

- n.rep:

  Integer; the number of replicated Z-residual sets to generate. Default
  is `1`.

- data:

  Optional data frame used to override the data stored inside `fit` for
  prediction and diagnostic calculation. If `NULL`, the data embedded in
  `fit` are used.

- type:

  Optional character string controlling the residual type; the meaning
  is determined by the underlying implementation (if used).

- ...:

  Further arguments passed to lower-level helpers.

## Value

A numeric matrix of Z-residuals with attributes:

- `type`: Type of outcome (Bernoulli).

- `zero_id`: Indices of zero outcomes.

- `log_pmf`: Log-probability mass function values.

- `log_cdf`: Log-cumulative distribution function values.

- `covariates`: Model covariates.

- `linear.pred`: Linear predictor values from the fitted model.

## Details

The function typically performs the following steps:

1.  Extracts the observed response vector from the model data.

2.  Computes the log-PMF and log-CDF for the Bernoulli model using
    [`log_pred_dist_bern`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_bern.md).

3.  Generates posterior predictive p-values according to the specified
    `method`.

4.  Converts the p-values to Z-residuals via the negative quantile of
    the standard normal distribution.

The output is a matrix of Z-residuals with one column per replication.

## See also

[`log_pred_dist_bern`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_bern.md),
[`post_logrpp`](https://tiw150.github.io/Zresidual/reference/post_logrpp.md),
[`iscv_logrpp`](https://tiw150.github.io/Zresidual/reference/iscv_logrpp.md),
and the S3 wrapper
[`Zresidual.bernoulli.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.bernoulli.brms.md).
