# Probability Mass Function of the Hurdle Poisson Distribution

Computes the probability density (or log-density) for the Poisson hurdle
distribution. This distribution combines a point mass at zero with a
truncated-at-zero Poisson distribution for positive counts.

## Usage

``` r
dhurdle.pois(y, lambda, pi, log = FALSE)
```

## Arguments

- y:

  Numeric vector of observed counts.

- lambda:

  Numeric vector of Poisson mean parameters (must be positive).

- pi:

  Numeric vector of hurdle probabilities (probability of structural
  zeros), where each value must be between 0 and 1.

- log:

  Logical; if `TRUE`, probabilities are returned on the log scale.

## Value

A numeric vector of the same length as `y`, giving the density (or
log-density) of the Poisson hurdle distribution.

## Details

The hurdle Poisson distribution assumes: \$\$ P(Y = 0) = \pi \$\$ and
for \\y \> 0\\: \$\$ P(Y = y) = (1 - \pi) \frac{P\_{\text{Pois}}(Y =
y)}{1 - P\_{\text{Pois}}(Y = 0)} \$\$ where \\P\_{\text{Pois}}(Y = y)\\
is the standard Poisson probability mass function.

The function is vectorized over all parameters.

## See also

[`dpois`](https://rdrr.io/r/stats/Poisson.html),
[`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html),
[`dhurdle.nb`](https://tiw150.github.io/Zresidual/reference/dhurdle.nb.md)

## Examples

``` r
# Example usage:
y <- 0:5
lambda <- 2
pi <- 0.3
dhurdle.pois(y, lambda, pi)
#> [1] 0.30000000 0.21912470 0.21912470 0.14608313 0.07304157 0.02921663
dhurdle.pois(y, lambda, pi, log = TRUE)
#> [1] -1.203973 -1.518114 -1.518114 -1.923579 -2.616727 -3.533017
```
