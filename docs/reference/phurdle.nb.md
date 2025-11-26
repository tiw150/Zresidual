# Cumulative Distribution Function of the Hurdle Negative Binomial Distribution

Computes the cumulative distribution function (CDF) or its logarithm for
the hurdle negative binomial (HNB) distribution. The hurdle model
combines a point mass at zero with a truncated negative binomial
distribution for positive counts.

## Usage

``` r
phurdle.nb(y, mu, size, pi, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of observed count values.

- mu:

  Numeric vector of mean parameters of the negative binomial
  distribution.

- size:

  Numeric vector of shape (dispersion) parameters of the negative
  binomial distribution.

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

The hurdle negative binomial model assumes: \$\$ P(Y = 0) = \pi, \quad
P(Y = y \mid Y \> 0) = (1 - \pi) \frac{F\_{NB}(y)-F\_{NB}(0)}{1 -
F\_{NB}(0)}, \quad y \> 0 \$\$ where \\F\_{NB}(y)\\ is CDF of the
standard negative binomial distribution.

The function computes the upper or lower tail probabilities for both
zeros and positive counts using the logarithmic form for numerical
stability. Internal helper functions (`log_diff_exp`, `log_sum_exp`) are
used to handle differences and sums of log-scale probabilities safely.

## See also

[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`pdf.tnb`](https://tiw150.github.io/Zresidual/reference/pdf.tnb.md) for
the zero-truncated negative binomial PMF.

## Examples

``` r
# Example: Hurdle Negative Binomial CDF
y <- 0:5
mu <- 2
size <- 1.5
pi <- 0.3
phurdle.nb(y, mu, size, pi)
#> Error in phurdle.nb(y, mu, size, pi): could not find function "phurdle.nb"

# Upper tail probabilities on log scale
phurdle.nb(y, mu, size, pi, lower.tail = FALSE, log.p = TRUE)
#> Error in phurdle.nb(y, mu, size, pi, lower.tail = FALSE, log.p = TRUE): could not find function "phurdle.nb"
```
