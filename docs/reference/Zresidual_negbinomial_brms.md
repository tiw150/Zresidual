# Compute Z-residuals for negative binomial brms models

Computes Z-residuals for fitted Bayesian negative binomial models.
Z-residuals are useful for model diagnostics, including checking fit and
overdispersion, and can be calculated using posterior or cross-validated
predictive p-values.

This is an internal workhorse for \[Zresidual.negbinomial.brms()\] and
is not intended to be called directly by end users.

## Usage

``` r
Zresidual_negbinomial_brms(fit, method = "iscv", n.rep = 1)
```

## Arguments

- fit:

  A fitted brms model object for a negative binomial outcome.

- method:

  Character string specifying the residual calculation method: `"iscv"`
  for importance-sampled cross-validated randomized predictive p-values,
  `"rpost"` for posterior predictive p-values, or `"mpost"` for marginal
  posterior predictive p-values. Default is `"iscv"`.

- n.rep:

  Integer; the number of replicated Z-residual sets to generate. Default
  is 1.

- ...:

  Further arguments passed to lower-level helper functions (if any).

## Value

A numeric matrix of Z-residuals with attributes such as:

- `zero_id`: Indices of zero outcomes.

- `log_pmf`: Log-probability mass function values.

- `log_cdf`: Log-cumulative distribution function values.

- `covariates`: Model covariates.

- `linear.pred`: Linear predictor values from the fitted model.

The S3 wrapper \[Zresidual.negbinomial.brms()\] will additionally attach
the class `"zresid"` to the returned object.

## Details

The function typically performs the following steps:

1.  Extracts the observed response vector from the model data.

2.  Computes the log-PMF and log-CDF for the negative binomial model
    using
    [`log_pred_dist_NB`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_NB.md).

3.  Generates randomized or posterior predictive p-values according to
    the specified `method`.

4.  Converts the p-values to Z-residuals via the negative quantile of
    the standard normal distribution.

The output is a matrix of Z-residuals with one column per replication.

## See also

[`log_pred_dist_NB`](https://tiw150.github.io/Zresidual/reference/log_pred_dist_NB.md),
[`post_logrpp`](https://tiw150.github.io/Zresidual/reference/post_logrpp.md),
[`iscv_logrpp`](https://tiw150.github.io/Zresidual/reference/iscv_logrpp.md),
and the S3 wrapper \[Zresidual.negbinomial.brms()\].

## Examples

``` r
if (FALSE) { # \dontrun{
  # Compute Z-residuals for a negative binomial model
  zres_nb <- Zresidual_negbinomial_brms(
    fit    = fit_nb,
    method = "iscv"
  )

  # Compute Z-residuals with 2 replicates using posterior predictive p-values
  zres_nb_post <- Zresidual_negbinomial_brms(
    fit    = fit_nb,
    method = "rpost",
    n.rep  = 2
  )
} # }
```
