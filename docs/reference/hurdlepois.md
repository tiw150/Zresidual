# Hurdle Poisson distribution

Density and distribution functions for the hurdle Poisson distribution.

## Usage

``` r
dhurdlepois(y, lambda, pi, log = FALSE)

phurdlepois(y, lambda, pi, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of counts.

- lambda:

  Mean parameter of the Poisson component.

- pi:

  Probability of a structural zero.

- log:

  Logical; used only by `dhurdlepois()`. If `TRUE`, return
  log-densities.

- lower.tail:

  Logical; used only by `phurdlepois()`. If `TRUE`, return \\P(Y \le
  y)\\; otherwise return \\P(Y \> y)\\.

- log.p:

  Logical; used only by `phurdlepois()`. If `TRUE`, return probabilities
  on the log scale.

## Value

A numeric vector of densities, distribution values, or their logarithms.

## Examples

``` r
dhurdlepois(0:3, lambda = 2, pi = 0.3)
#> [1] 0.3000000 0.2191247 0.2191247 0.1460831
phurdlepois(0:3, lambda = 2, pi = 0.3)
#> [1] 0.7000000 0.4808753 0.2617506 0.1156675
```
