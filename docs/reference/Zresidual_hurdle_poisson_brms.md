# Compute Z-residuals for hurdle or count Poisson brms models

Computes Z-residuals for fitted Bayesian hurdle or count models with a
Poisson distribution using a brms model with
`family = hurdle_poisson()`. Z-residuals can be calculated for zeros,
counts, or the overall hurdle distribution, and can be used for model
diagnostics.

This is an internal workhorse for
[`Zresidual.hurdle_poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.hurdle_poisson.brms.md)
and is not intended to be called directly by end users.

## Usage

``` r
Zresidual_hurdle_poisson_brms(fit, type, method = "iscv", n.rep = 1, ...)
```

## Arguments

- fit:

  A fitted brms model object for a hurdle or count Poisson outcome.

- type:

  Character string specifying which part of the model to calculate
  Z-residuals for: `"zero"` for the hurdle/zero portion, `"count"` for
  the truncated Poisson counts, `"hurdle"` for the full hurdle-Poisson
  model.

- method:

  Character string specifying the residual calculation method: `"iscv"`
  for importance-sampled cross-validated randomized predictive p-values,
  `"rpost"` for randomized posterior predictive p-values, or `"mpost"`
  for middle-value posterior predictive p-values. Default is `"iscv"`.

- n.rep:

  Integer; the number of replicated Z-residual sets to generate. Default
  is `1`.

## Value

A numeric matrix of Z-residuals with attributes such as:

- `type`: The requested model component.

- `zero_id`: Indices of zero outcomes.

- `log_pmf`: Log-probability mass function values.

- `log_cdf`: Log-cumulative distribution function values.

- `covariates`: Model covariates.

- `linear.pred`: Linear predictor values from the fitted model.

The S3 wrapper
[`Zresidual.hurdle_poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.hurdle_poisson.brms.md)
will additionally attach the class `"zresid"` to the returned object.

## Details

A typical implementation:

1.  Extracts the observed response vector from the model data.

2.  Computes the log-PMF and log-CDF for the specified part of the model
    using the corresponding `log_pred_dist_*` function, such as
    [`log_pred_dist_HP`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_HP.md)
    or
    [`log_pred_dist_TP`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_TP.md).

3.  Generates posterior predictive p-values according to the specified
    `method`.

4.  Converts the p-values to Z-residuals via the negative quantile of
    the standard normal distribution.

The output is a matrix of Z-residuals with one column per replication.

## See also

[`log_pred_dist_HP`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_HP.md),
[`log_pred_dist_TP`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_TP.md),
[`post_logrpp`](https://tiw150.github.io/Zresidual/reference/post_logrpp.md),
[`iscv_logrpp`](https://tiw150.github.io/Zresidual/reference/iscv_logrpp.md),
and the S3 wrapper
[`Zresidual.hurdle_poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.hurdle_poisson.brms.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Compute Z-residuals for counts
  zres_counts <- Zresidual_hurdle_poisson_brms(
    fit    = fit_hp,
    type   = "count",
    method = "iscv"
  )
} # }
```
