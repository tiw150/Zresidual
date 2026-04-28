# Compute observation-wise predictive quantities for a fitted Gaussian `stats::glm` model on a user-supplied dataset.

Compute observation-wise predictive quantities for a fitted Gaussian
[`stats::glm`](https://rdrr.io/r/stats/glm.html) model on a
user-supplied dataset.

## Usage

``` r
log_pointpred_glm_gaussian(fit, data, ...)
```

## Arguments

- fit:

  A fitted [`stats::glm`](https://rdrr.io/r/stats/glm.html) object with
  Gaussian family.

- data:

  A data frame on which predictive quantities are evaluated. This
  argument must be explicitly provided.

- ...:

  Additional arguments reserved for future extensions.

## Value

A list with components:

- log_like:

  A 1 x n matrix of observation-wise log predictive densities.

- log_surv:

  A 1 x n matrix of log survival probabilities, \\\log P(Y \> y_i)\\.

- is_discrete:

  An integer vector of length n, equal to 0 for all observations because
  the Gaussian outcome is continuous.
