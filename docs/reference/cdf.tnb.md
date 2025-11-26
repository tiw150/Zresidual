# Cumulative Distribution Function (CDF) of Zero-Truncated Negative Binomial Distribution

Computes the cumulative distribution function (CDF) for a zero-truncated
negative binomial (TNB) distribution. This function adjusts the standard
negative binomial CDF to account for truncation at zero (i.e., only
supports \\y \> 0\\).

## Usage

``` r
cdf.tnb(y, mu, size, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of quantiles (count values) for which to compute the
  CDF.

- mu:

  Mean parameter (\\\mu\\) of the negative binomial distribution.

- size:

  Dispersion (shape) parameter (\\r\\) of the negative binomial
  distribution.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(Y \le y)\\. If
  `FALSE`, probabilities are \\P(Y \> y)\\.

- log.p:

  Logical; if `TRUE`, probabilities \\p\\ are given as \\\log(p)\\.

## Value

A numeric vector of the same length as the input, containing:

- CDF values (\\P(Y \le y \mid Y \> 0)\\) if `lower.tail = TRUE`.

- Upper-tail probabilities (\\P(Y \> y \mid Y \> 0)\\) if
  `lower.tail = FALSE`.

- Log-probabilities if `log.p = TRUE`.

## Details

The function computes probabilities for the zero-truncated version of
the negative binomial distribution: \$\$P(Y \le y \mid Y \> 0) =
\frac{P(Y \le y) - P(Y = 0)}{1 - P(Y = 0)}.\$\$

Internally, this is implemented using the log-scale for numerical
stability. When `lower.tail = FALSE`, it computes the upper-tail
probabilities \\P(Y \> y \mid Y \> 0)\\ instead.

## See also

[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html) for the standard
negative binomial CDF.

## Examples

``` r
# Example: Compute the CDF for y = 1:5
mu <- 2
size <- 1
cdf.tnb(1:5, mu, size)
#> Error in cdf.tnb(1:5, mu, size): could not find function "cdf.tnb"

# Compute the upper-tail probabilities
cdf.tnb(1:5, mu, size, lower.tail = FALSE)
#> Error in cdf.tnb(1:5, mu, size, lower.tail = FALSE): could not find function "cdf.tnb"

# Log probabilities
cdf.tnb(1:5, mu, size, log.p = TRUE)
#> Error in cdf.tnb(1:5, mu, size, log.p = TRUE): could not find function "cdf.tnb"
```
