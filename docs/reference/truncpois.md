# Zero-truncated Poisson distribution

Density and distribution functions for the zero-truncated Poisson
distribution.

## Usage

``` r
dtruncpois(y, lambda, log.p = FALSE)

ptruncpois(y, lambda, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of counts.

- lambda:

  Mean parameter of the Poisson distribution.

- log.p:

  Logical; if `TRUE`, return values on the log scale. For `dtruncpois()`
  this returns log-densities; for `ptruncpois()` this returns
  log-probabilities.

- lower.tail:

  Logical; used only by `ptruncpois()`. If `TRUE`, return \\P(Y \le
  y)\\; otherwise return \\P(Y \> y)\\.

## Value

A numeric vector of densities, distribution values, or their logarithms.

## Examples

``` r
dtruncpois(1:4, lambda = 2)
#> [1] 0.3130353 0.3130353 0.2086902 0.1043451
ptruncpois(1:4, lambda = 2)
#> [1] 0.68696471 0.37392943 0.16523924 0.06089414
```
