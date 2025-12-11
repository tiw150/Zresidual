# Compute Log Predictive Distributions for a Truncated Poisson Model of a 'brms' Fit

This function calculates the log predictive mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted truncated Poisson model (fitted using brms). The function
extracts posterior samples for the model’s mean parameter and evaluates
the predictive distributions across all posterior draws.

## Usage

``` r
log_pred_dist_TP(fit)
```

## Arguments

- fit:

  A fitted brms truncated Poisson model object. The model must include
  the distributional parameter `mu` (mean parameter).

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CDF values (posterior samples × observations).

## Details

For each posterior draw and observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`pdf.tp()`](https://tiw150.github.io/Zresidual/reference/pdf.tp.md).

- `lcdf_hat`: Log cumulative distribution function values using
  `cdf.tp.li()` with `lower.tail = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y | trunc(lb = 1) ~ x1 + x2), family = poisson(), data = mydata)
pred_dist <- log_pred_dist_TP(fit)
str(pred_dist)
} # }
```
