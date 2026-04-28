# Hurdle negative binomial distribution

Density and distribution functions for the hurdle negative binomial
distribution.

## Usage

``` r
dhurdlenb(y, mu, size, pi, log = FALSE)

phurdlenb(y, mu, size, pi, lower.tail = FALSE, log.p = FALSE)
```

## Arguments

- y:

  Numeric vector of counts.

- mu:

  Mean parameter of the negative binomial component.

- size:

  Dispersion parameter of the negative binomial component.

- pi:

  Probability of a structural zero.

- log:

  Logical; used only by `dhurdlenb()`. If `TRUE`, return log-densities.

- lower.tail:

  Logical; used only by `phurdlenb()`. If `TRUE`, return \\P(Y \le y)\\;
  otherwise return \\P(Y \> y)\\.

- log.p:

  Logical; used only by `phurdlenb()`. If `TRUE`, return probabilities
  on the log scale.

## Value

A numeric vector of densities, distribution values, or their logarithms.

## Examples

``` r
dhurdlenb(0:3, mu = 2, size = 1.5, pi = 0.3)
#> [1] 0.3000000 0.2339888 0.1671348 0.1114232
phurdlenb(0:3, mu = 2, size = 1.5, pi = 0.3)
#> [1] 0.7000000 0.4660112 0.2988764 0.1874532
```
