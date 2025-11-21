# Probability Mass Function of the Zero-Truncated Negative Binomial Distribution

Computes the probability mass function (PMF) or log-PMF for the
zero-truncated negative binomial (TNB) distribution. This version
excludes zeros and rescales the probabilities accordingly so that they
sum to one over positive counts only.

## Usage

``` r
pdf.tnb(y, mu, size, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of observed count values (`y > 0`).

- mu:

  Numeric vector of mean parameters for the negative binomial
  distribution.

- size:

  Numeric vector of shape (dispersion) parameters.

- log.p:

  Logical; if `TRUE`, returns log probabilities instead of
  probabilities.

## Value

A numeric vector of probabilities (or log-probabilities if
`log.p = TRUE`).

## Details

The zero-truncated negative binomial probability for an observation \\y
\> 0\\ is: \$\$ P(Y = y \mid Y \> 0) = \frac{P(Y = y)}{1 - P(Y = 0)}
\$\$ where \\P(Y = y)\\ and \\P(Y = 0)\\ are evaluated using the
standard negative binomial PMF and CDF, respectively. The implementation
uses [`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) and
[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html) for computation.

The function automatically vectorizes inputs, ensuring that the output
corresponds elementwise to each set of parameters.

## See also

[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html), and
[`cdf.tnb`](https://tiw150.github.io/Zresidual/reference/cdf.tnb.md) for
the corresponding cumulative function.

## Examples

``` r
# Example: Zero-truncated negative binomial probabilities
y <- 1:5
mu <- 2
size <- 1.5
pdf.tnb(y, mu, size)
#> Error in pdf.tnb(y, mu, size): could not find function "pdf.tnb"

# Log probabilities
pdf.tnb(y, mu, size, log.p = TRUE)
#> Error in pdf.tnb(y, mu, size, log.p = TRUE): could not find function "pdf.tnb"
```
