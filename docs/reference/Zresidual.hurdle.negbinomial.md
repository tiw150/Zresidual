# Compute Z-Residuals for Hurdle or Count Negative Binomial Models

Computes Z-residuals for fitted Bayesian hurdle or count models with a
negative binomial distribution. Z-residuals can be calculated for zeros,
counts, or the overall hurdle distribution, and can be used for model
diagnostics.

## Usage

``` r
Zresidual.hurdle.negbinomial(fit, type, method = "iscv", n.rep = 1)
```

## Arguments

- fit:

  A fitted brms model object for a hurdle or count negative binomial
  outcome.

- type:

  Character string specifying which part of the model to calculate
  Z-residuals for: `"zero"` for the hurdle/zero portion, `"count"` for
  the truncated negative binomial counts, `"hurdle"` for the full
  hurdle-negative binomial model.

- method:

  Character string specifying the residual calculation method: `"iscv"`
  for importance-sampled cross-validated randomized predictive p-values,
  or `"rpost"` for randomized posterior predictive p-values or `"mpost"`
  for middle-value posterior predictive p-values. Default is `"iscv"`.

- n.rep:

  Integer; the number of replicated Z-residual sets to generate. Default
  is 1.

## Value

A numeric matrix of Z-residuals with attributes:

- `type`: The specified model component

- `zero_id`: Indices of zero outcomes

- `log_pmf`: Log-probability mass function values

- `log_cdf`: Log-cumulative distribution function values

- `covariates`: Model covariates

- `linear.pred`: Linear predictor values from the fitted model

The returned object has class `c("zresid", "matrix")`.

## Details

The function performs the following steps:

1.  Extracts the observed response vector from the model data.

2.  Computes the log-PMF and log-CDF for the specified part of the model
    using the corresponding `log.pred.dist.*` function.

3.  Generates posterior predictive p-values according to the specified
    `method`.

4.  Converts the p-values to Z-residuals via the negative quantile of
    the standard normal distribution.

The output is a matrix of Z-residuals with one column per replication.

## See also

[`log.pred.dist.HNB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.HNB.md),
[`log.pred.dist.TNB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.TNB.md),
[`post_logrpp`](https://tiw150.github.io/Zresidual/reference/post_logrpp.md),
[`iscv_logrpp`](https://tiw150.github.io/Zresidual/reference/iscv_logrpp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute Z-residuals for counts
zres_counts <- Zresidual.hurdle.negbinomial(fit, type = "count", method = "iscv")

# Compute Z-residuals for the full hurdle model with 2 replicates
zres_hurdle <- Zresidual.hurdle.negbinomial(fit, type = "hurdle", n.rep = 2)
} # }
```
