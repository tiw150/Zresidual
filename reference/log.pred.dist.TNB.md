# Compute Log Predictive Distributions for a Truncated Negative Binomial Model of a 'brms' Fit

This function calculates the log predictive mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted truncated negative binomial model (fitted using brms). The
function extracts posterior samples for the model’s distributional
parameters and evaluates the predictive distributions across all
posterior draws.

## Usage

``` r
# S3 method for class 'pred.dist.TNB'
log(fit)
```

## Arguments

- fit:

  A fitted brms truncated negative binomial model object. The model must
  include the distributional parameters `mu` (mean parameter) and
  `shape` (dispersion parameter).

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CDF values (posterior samples × observations).

## Details

For each posterior draw and observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`pdf.tnb()`](https://tiw150.github.io/Zresidual/reference/pdf.tnb.md).

- `lcdf_hat`: Log cumulative distribution function values using
  [`cdf.tnb()`](https://tiw150.github.io/Zresidual/reference/cdf.tnb.md)
  with `lower.tail = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y | trunc(lb = 1) ~ x1 + x2), family = negbinomial(), data = mydata)
pred_dist <- log.pred.dist.TNB(fit)
str(pred_dist)
} # }
```
