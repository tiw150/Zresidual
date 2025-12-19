# Probability Mass Function of the Zero-Truncated Poisson Distribution

Computes the probability mass function (PMF) or log-PMF for the
zero-truncated Poisson (TP) distribution. This version excludes zeros
and rescales the probabilities so that they sum to one over positive
counts only.

## Usage

``` r
pdf.tp(y, lambda, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of observed count values (`y > 0`).

- lambda:

  Numeric vector of rate parameters (mean of the Poisson distribution).

- log.p:

  Logical; if `TRUE`, returns log probabilities instead of
  probabilities.

## Value

A numeric vector of probabilities (or log-probabilities if
`log.p = TRUE`).

## Details

The zero-truncated Poisson probability for an observation \\y \> 0\\ is:
\$\$ P(Y = y \mid Y \> 0) = \frac{P(Y = y)}{1 - P(Y = 0)} \$\$ where
\\P(Y = y)\\ and \\P(Y = 0)\\ are evaluated using the standard Poisson
PMF and CDF, respectively. The function uses
[`dpois`](https://rdrr.io/r/stats/Poisson.html) and
[`ppois`](https://rdrr.io/r/stats/Poisson.html) internally.

This function automatically vectorizes inputs so that each probability
corresponds elementwise to the provided parameter values.

## See also

[`dpois`](https://rdrr.io/r/stats/Poisson.html),
[`ppois`](https://rdrr.io/r/stats/Poisson.html), and
[`cdf.tp`](https://tiw150.github.io/Zresidual/reference/cdf.tp.md) for
the corresponding cumulative function.

## Examples

``` r
# Example: Zero-truncated Poisson probabilities
y <- 1:5
lambda <- 2
pdf.tp(y, lambda)
#> [1] 0.31303529 0.31303529 0.20869019 0.10434510 0.04173804

# Log probabilities
pdf.tp(y, lambda, log.p = TRUE)
#> [1] -1.161439 -1.161439 -1.566904 -2.260052 -3.176342
```
