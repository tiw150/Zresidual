# Probability Mass Function of the Hurdle Negative Binomial Distribution

Computes the probability (or log-probability) of observing a count \\y\\
under a hurdle model with a Bernoulli process for zeros and a truncated
Negative Binomial distribution for positive counts (\\y \> 0\\).

## Usage

``` r
dhurdlenb(y, mu, size, pi, log = FALSE)
```

## Arguments

- y:

  Numeric vector of observed count values.

- mu:

  Mean parameter (\\\mu\\) of the Negative Binomial component.

- size:

  Dispersion (shape) parameter (\\r\\) of the Negative Binomial
  component.

- pi:

  Probability of observing a structural zero (from the hurdle
  component).

- log:

  Logical; if `TRUE`, probabilities \\p\\ are returned on the log scale.

## Value

A numeric vector of:

- Probabilities \\P(Y = y)\\ if `log = FALSE`.

- Log-probabilities \\\log P(Y = y)\\ if `log = TRUE`.

## Details

The hurdle Negative Binomial (HNB) distribution models counts \\Y\\ as:

\$\$ P(Y = 0) = \pi, \$\$ \$\$ P(Y = y) = (1 - \pi)
\frac{f\_{NB}(y)}{1 - f\_{NB}(0)}, \quad \text{for } y \> 0, \$\$

where \\f\_{NB}(y)\\ is the Negative Binomial probability mass function
with parameters \\\mu\\ and \\r\\ (dispersion).

Internally, computations are performed on the log scale for numerical
stability. The function handles vectorized inputs and returns a numeric
vector of the same length as `y`.

## See also

[`dnbinom`](https://rdrr.io/r/stats/NegBinomial.html) for the Negative
Binomial PMF, [`pnbinom`](https://rdrr.io/r/stats/NegBinomial.html) for
its CDF.

## Examples

``` r
# Example parameters
mu <- 2
size <- 1
pi <- 0.3

# Compute hurdle NB probabilities for counts 0–5
dhurdlenb(0:5, mu, size, pi)
#> [1] 0.30000000 0.23333333 0.15555556 0.10370370 0.06913580 0.04609053

# Log-scale probabilities
dhurdlenb(0:5, mu, size, pi, log = TRUE)
#> [1] -1.203973 -1.455287 -1.860752 -2.266217 -2.671683 -3.077148

# Compare structural zero vs positive counts
dhurdlenb(0, mu, size, pi)      # P(Y = 0)
#> [1] 0.3
dhurdlenb(1, mu, size, pi)      # P(Y = 1)
#> [1] 0.2333333
```
