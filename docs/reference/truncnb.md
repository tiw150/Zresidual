# Zero-truncated negative binomial distribution

Density and distribution functions for the zero-truncated negative
binomial distribution.

## Usage

``` r
dtruncnb(y, mu, size, log.p = FALSE)

ptruncnb(y, mu, size, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of counts.

- mu:

  Mean parameter of the negative binomial distribution.

- size:

  Dispersion parameter of the negative binomial distribution.

- log.p:

  Logical; if `TRUE`, return values on the log scale. For `dtruncnb()`
  this returns log-densities; for `ptruncnb()` this returns
  log-probabilities.

- lower.tail:

  Logical; used only by `ptruncnb()`. If `TRUE`, return \\P(Y \le y)\\;
  otherwise return \\P(Y \> y)\\.

## Value

A numeric vector of densities, distribution values, or their logarithms.

## Examples

``` r
dtruncnb(1:4, mu = 2, size = 1.5)
#> [1] 0.3342697 0.2387641 0.1591760 0.1023275
ptruncnb(1:4, mu = 2, size = 1.5)
#> [1] 0.6657303 0.4269663 0.2677902 0.1654628
```
