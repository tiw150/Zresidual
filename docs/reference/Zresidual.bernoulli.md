# Compute Z-Residuals for a Bernoulli/Logistic Model

Computes Z-residuals for a fitted Bayesian Bernoulli/Logistic (binary)
model. Z-residuals are calculated using posterior predictive methods and
can be used for model diagnostics.

## Usage

``` r
Zresidual.bernoulli(fit, method = "iscv", n.rep = 1)
```

## Arguments

- fit:

  A fitted brms model object.

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

- `type`: Type of outcome (Bernoulli)

- `zero_id`: Indices of zero outcomes

- `log_pmf`: Log-probability mass function values

- `log_cdf`: Log-cumulative distribution function values

- `covariates`: Model covariates

- `linear.pred`: Linear predictor values from the fitted model

The returned object has class `c("zresid", "matrix")`.

## Details

The function performs the following steps:

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
[`iscv_logrpp`](https://tiw150.github.io/Zresidual/reference/iscv_logrpp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute Z-residuals using ISCV method
zres <- Zresidual.bernoulli(fit, method = "iscv", n.rep = 2)

# Compute Z-residuals using posterior predictive method
zres_post <- Zresidual.bernoulli(fit, method = "post")
} # }
```
