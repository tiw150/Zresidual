# Compute Log Predictive Distributions for a Poisson Model of a 'brms' Fit

This function calculates the log predictive mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted Poisson model (fitted using brms). The function extracts
posterior samples for the model’s mean parameter and evaluates the
predictive distributions across all posterior draws.

## Usage

``` r
# S3 method for class 'pred.dist.pois'
log(fit)
```

## Arguments

- fit:

  A fitted brms Poisson model object. The model must include the
  distributional parameter `mu` (mean parameter).

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CDF values (posterior samples × observations).

- `zero_id`:

  Indices of observations with zero counts.

## Details

For each posterior draw and observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`dpois()`](https://rdrr.io/r/stats/Poisson.html).

- `lcdf_hat`: Log cumulative distribution function values using
  [`ppois()`](https://rdrr.io/r/stats/Poisson.html) with
  `lower.tail = FALSE`.

The function also identifies indices of zero-valued observations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y ~ x1 + x2), family = poisson(), data = mydata)
pred_dist <- log.pred.dist.pois(fit)
str(pred_dist)
} # }
```
