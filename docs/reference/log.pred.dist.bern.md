# Compute Log Predictive Distributions for Logistic Regression of a 'brms' Fit

Calculates the log probability mass function (log-PMF) and log
cumulative distribution function (log-CDF) for a logistic model based on
a fitted 'brms' model.

## Usage

``` r
# S3 method for class 'pred.dist.bern'
log(fit)
```

## Arguments

- fit:

  A fitted model object from `brms`.

## Value

A list containing:

- `lpmf_hat`: numeric matrix of log-PMF values (posterior draws ×
  observations).

- `lcdf_hat`: numeric matrix of log-CDF values (posterior draws ×
  observations).

## Details

The function extracts the posterior predictions of the Bernoulli
component (structural zeros) using
[`posterior.pred()`](https://tiw150.github.io/Zresidual/reference/posterior.pred.md)
for the mean parameter ("mu").

For each observation:

- Computes the log-PMF of the observed binary outcome.

- Computes the log-CDF (upper-tail probability) of the observed outcome.

This produces matrices of size \\M \times N\\, where \\M\\ is the number
of posterior draws and \\N\\ is the number of observations.

## See also

[`posterior.pred`](https://tiw150.github.io/Zresidual/reference/posterior.pred.md),
`dbern`, `pbern`

## Examples

``` r
# Assuming 'fit' is a fitted brms logistic model
# pred_dist <- log.pred.dist.bern(fit)
# lpmf_hat <- pred_dist$lpmf_hat
# lcdf_hat <- pred_dist$lcdf_hat
```
