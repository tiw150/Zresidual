# Compute Posterior Log Randomized Predictive p-values (Log-RPP)

Calculates the posterior log randomized predictive p-values (log-RPP)
for a set of observations given the log-PMF and log-CDF matrices from
posterior predictive samples.

## Usage

``` r
post_logrpp(log_cdf, log_pmf)
```

## Arguments

- log_cdf:

  Numeric matrix of log cumulative probabilities for each observation
  (posterior samples × observations).

- log_pmf:

  Numeric matrix of log probability mass function values for each
  observation (posterior samples × observations).

## Value

A numeric vector of log randomized posterior predictive p-values for
each observation.

## Details

The function computes a randomized posterior predictive p-values as
follows:

1.  A uniform random number is generated for each observation to
    randomize the probability for zero counts.

2.  Probabilities exactly equal to 0 or 1 are replaced with small bounds
    (\\10^{-5}\\ and \\9 \cdot 10^{-5}\\) to avoid numerical issues.

## See also

[`post_logmpp`](https://tiw150.github.io/Zresidual/reference/post_logmpp.md),
[`log.pred.dist.HNB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.HNB.md),
[`log.pred.dist.NB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.NB.md),
[`log.pred.dist.TNB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.TNB.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume log_cdf and log_pmf are matrices from posterior predictive draws
log_rpp <- post_logrpp(log_cdf, log_pmf)
head(exp(log_rpp))  # Randomized posterior predictive probabilities
} # }
```
