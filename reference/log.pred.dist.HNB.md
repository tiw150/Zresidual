# Compute Log Predictive Distributions for a Hurdle Negative Binomial Model of a 'brms' Fit

This function calculates the log probability mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted hurdle negative binomial model (fitted using brms). The
function extracts posterior samples for the model parameters and
evaluates the predictive distributions across all posterior draws.

## Usage

``` r
# S3 method for class 'pred.dist.HNB'
log(fit)
```

## Arguments

- fit:

  A fitted brms hurdle negative binomial model object. The model must
  include the distributional parameters `mu`, `shape`, and the hurdle
  probability `zero`.

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CCDF values (posterior samples × observations).

- `zero_id`:

  Indices of observations with zero counts.

- `count_id`:

  Indices of observations with positive counts.

## Details

For each posterior draw and each observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`dhurdle.nb()`](https://tiw150.github.io/Zresidual/reference/dhurdle.nb.md).

- `lcdf_hat`: Log cumulative distribution function values using
  [`phurdle.nb()`](https://tiw150.github.io/Zresidual/reference/phurdle.nb.md)
  with `lower.tail = FALSE`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y ~ x1 + x2, hu ~ x1), family = hurdle_negbinomial(), data = mydata)
pred_dist <- log.pred.dist.HNB(fit)
str(pred_dist)
} # }
```
