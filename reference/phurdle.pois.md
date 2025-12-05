# Cumulative Distribution Function of the Hurdle Poisson Distribution

Computes the cumulative distribution function (CDF) or its logarithm for
the hurdle Poisson (HP) distribution. The hurdle model combines a point
mass at zero with a truncated poisson distribution for positive counts.

## Usage

``` r
phurdle.pois(y, lambda, pi, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of observed count values.

- lambda:

  Numeric vector of mean parameters of the poisson distribution.

- pi:

  Numeric vector of hurdle probabilities (probability of structural
  zeros).

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le y)\\;
  otherwise, they are \\P(Y \> y)\\.

- log.p:

  Logical; if `TRUE`, probabilities are returned on the log scale.

## Value

A numeric vector of cumulative probabilities (or log-probabilities if
`log.p = TRUE`).

## Details

The hurdle poisson model assumes: \$\$ P(Y = 0) = \pi, \quad P(Y = y
\mid Y \> 0) = (1 - \pi) \frac{F\_{Pois}(y)-F\_{Pois}(0)}{1 -
F\_{Pois}(0)}, \quad y \> 0 \$\$ where \\F\_{Pois}(y)\\ is CDF of the
standard possion distribution.

The function computes the upper or lower tail probabilities for both
zeros and positive counts using the logarithmic form for numerical
stability. Internal helper functions (`log_diff_exp`, `log_sum_exp`) are
used to handle differences and sums of log-scale probabilities safely.

## See also

[`ppois`](https://rdrr.io/r/stats/Poisson.html),
[`dpois`](https://rdrr.io/r/stats/Poisson.html),
[`pdf.tp`](https://tiw150.github.io/Zresidual/reference/pdf.tp.md) for
the zero-truncated Poisson PMF.

## Examples

``` r
# Example: Hurdle Poisson CDF
y <- 0:5
lambda <- 2
pi <- 0.3
phurdle.pois(y, lambda, pi)
#> [1] 0.70000000 0.48087530 0.26175060 0.11566747 0.04262590 0.01340927

# Upper tail probabilities on log scale
phurdle.pois(y, lambda, pi, lower.tail = FALSE, log.p = TRUE)
#> [1] -0.3566749 -0.7321473 -1.3403631 -2.1570359 -3.1552932 -4.3118087
```
