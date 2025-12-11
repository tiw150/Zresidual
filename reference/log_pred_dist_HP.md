# Compute Log Predictive Distributions for a Hurdle Poisson Model of a 'brms' Fit

This function calculates the log predictive mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted hurdle Poisson model (fitted using brms). The function
extracts posterior samples for the model parameters and evaluates the
predictive distributions across all posterior draws.

## Usage

``` r
log_pred_dist_HP(fit)
```

## Arguments

- fit:

  A fitted brms hurdle Poisson model object. The model must include the
  distributional parameters `mu` (mean parameter) and the hurdle
  probability `zero`.

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CDF values (posterior samples × observations).

- `zero_id`:

  Indices of observations with zero counts.

- `count_id`:

  Indices of observations with positive counts.

## Details

For each posterior draw and observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`dhurdle.pois()`](https://tiw150.github.io/Zresidual/reference/dhurdle.pois.md).

- `lcdf_hat`: Log cumulative distribution function values using
  [`phurdle.pois()`](https://tiw150.github.io/Zresidual/reference/phurdle.pois.md)
  with `lower.tail = FALSE`.

The function also identifies indices of zero and positive count
responses.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y ~ x1 + x2, hu ~ x1), family = hurdle_poisson(), data = mydata)
pred_dist <- log_pred_dist_HP(fit)
str(pred_dist)
} # }
```
