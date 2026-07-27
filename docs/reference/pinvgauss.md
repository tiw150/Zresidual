# Inverse Gaussian CDF and survival function

Distribution function for the inverse Gaussian distribution using the
same parameterization as brms' inverse.gaussian family: mean parameter
`mu` and shape parameter `shape`.

## Usage

``` r
pinvgauss(q, mu, shape, lower.tail = TRUE, log.p = FALSE)
```

## Arguments

- q:

  Numeric vector of quantiles.

- mu:

  Positive mean parameter.

- shape:

  Positive shape parameter.

- lower.tail:

  Logical; if `TRUE`, return `P(Y <= q)`, otherwise return `P(Y > q)`.

- log.p:

  Logical; if `TRUE`, return probabilities on the log scale.

## Value

A numeric vector of CDF/survival values or log-CDF/log-survival values.

## Details

This function is written to be stable on the log scale, especially for
right tail probabilities used by Z-residual calculations.
