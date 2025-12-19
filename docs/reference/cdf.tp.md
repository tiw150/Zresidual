# Cumulative Distribution Function (CDF) of Zero-Truncated Poisson Distribution

Computes the cumulative distribution function (CDF) for a zero-truncated
Poisson (TP) distribution. This function modifies the standard Poisson
CDF to account for truncation at zero (i.e., only supports \\y \> 0\\).

## Usage

``` r
cdf.tp(y, lambda, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of quantiles (count values) for which to compute the
  CDF.

- lambda:

  Numeric vector or scalar giving the mean parameter (\\\lambda\\) of
  the Poisson distribution.

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
the Poisson distribution: \$\$P(Y \le y \mid Y \> 0) = \frac{P(Y \le
y) - P(Y = 0)}{1 - P(Y = 0)}.\$\$

Internally, it uses the
[`ppois()`](https://rdrr.io/r/stats/Poisson.html) function from base R
for computing Poisson probabilities, and performs calculations on the
log scale for improved numerical stability.

When `lower.tail = FALSE`, it returns the upper-tail probabilities \\P(Y
\> y \mid Y \> 0)\\ instead.

## See also

[`ppois`](https://rdrr.io/r/stats/Poisson.html) for the standard Poisson
CDF.

## Examples

``` r
# Example: Compute the zero-truncated Poisson CDF for y = 1:5
lambda <- 2
cdf.tp(1:5, lambda)
#> [1] 0.68696471 0.37392943 0.16523924 0.06089414 0.01915611

# Compute upper-tail probabilities
cdf.tp(1:5, lambda, lower.tail = FALSE)
#> [1] 0.68696471 0.37392943 0.16523924 0.06089414 0.01915611

# Compute log-CDF values
cdf.tp(1:5, lambda, log.p = TRUE)
#> [1] -0.3754723 -0.9836882 -1.8003609 -2.7986183 -3.9551338
```
