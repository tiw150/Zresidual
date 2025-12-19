# Compute Posterior Log Middle-value Predictive p-values (Log-MPP)

Calculates the posterior log middle-value (0.5) predictive p-values
(log-MPP) for a set of observations given the log-PMF and log-CDF
matrices obtained from posterior predictive samples.

## Usage

``` r
post_logmpp(log_cdf, log_pmf)
```

## Arguments

- log_cdf:

  Numeric matrix of log cumulative probabilities for each observation
  (posterior samples × observations).

- log_pmf:

  Numeric matrix of log probability mass function values for each
  observation (posterior samples × observations).

## Value

A numeric vector of log posterior mid-value predictive probabilities for
each observation.

## Details

The function computes a stabilized version of the posterior predictive
p-values:

1.  It sums the log-CDF and log-PMF across posterior samples using
    log-sum-exp for numerical stability.

2.  Probabilities exactly equal to 0 or 1 are replaced with small bounds
    (\\10^{-5}\\ and \\9 \cdot 10^{-5}\\) to avoid numerical issues in
    further computations.

## Examples

``` r
NULL
#> NULL
```
