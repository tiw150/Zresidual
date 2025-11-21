# Compute Log Predictive Distributions for a Negative Binomial Model of a 'brms' Fit

This function calculates the log predictive mass function (log-PMF) and
the log cumulative distribution function (log-CDF) for each observation
from a fitted negative binomial model (fitted using brms). The function
extracts posterior samples for the distributional parameters and
evaluates the predictive distributions across all posterior draws.

## Usage

``` r
# S3 method for class 'pred.dist.NB'
log(fit)
```

## Arguments

- fit:

  A fitted brms negative binomial model object. The model must include
  the distributional parameters `mu` (mean parameter) and `shape`
  (dispersion parameter).

## Value

A list with the following components:

- `lpmf_hat`:

  A matrix of log-PMF values (posterior samples × observations).

- `lcdf_hat`:

  A matrix of log-CDF values (posterior samples × observations).

- `zero_id`:

  Indices of observations with zero counts.

## Details

For each posterior draw and each observation, the function computes:

- `lpmf_hat`: Log predictive mass function values using
  [`dnbinom()`](https://rdrr.io/r/stats/NegBinomial.html).

- `lcdf_hat`: Log cumulative distribution function values using
  [`pnbinom()`](https://rdrr.io/r/stats/NegBinomial.html) with
  `lower.tail = FALSE`.

The function also identifies indices of zero-valued observations, which
can be useful for model diagnostic procedures such as Z-residuals or
posterior predictive checks.

## See also

[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html), and
[`posterior_predict`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
for related computations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
fit <- brm(bf(y ~ x1 + x2), family = negbinomial(), data = mydata)
pred_dist <- log.pred.dist.NB(fit)
str(pred_dist)
} # }
```
