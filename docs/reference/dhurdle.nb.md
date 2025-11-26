# Probability Mass Function of the Hurdle Negative Binomial Distribution

Computes the probability (or log-probability) of observing a count \\y\\
under a hurdle model with a Bernoulli process for zeros and a truncated
Negative Binomial distribution for positive counts (\\y \> 0\\).

## Usage

``` r
dhurdle.nb(y, mu, size, pi, log = FALSE)
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
dhurdle.nb(0:5, mu, size, pi)
#> Error in dhurdle.nb(0:5, mu, size, pi): could not find function "dhurdle.nb"

# Log-scale probabilities
dhurdle.nb(0:5, mu, size, pi, log = TRUE)
#> Error in dhurdle.nb(0:5, mu, size, pi, log = TRUE): could not find function "dhurdle.nb"

# Compare structural zero vs positive counts
dhurdle.nb(0, mu, size, pi)      # P(Y = 0)
#> Error in dhurdle.nb(0, mu, size, pi): could not find function "dhurdle.nb"
dhurdle.nb(1, mu, size, pi)      # P(Y = 1)
#> Error in dhurdle.nb(1, mu, size, pi): could not find function "dhurdle.nb"
```
