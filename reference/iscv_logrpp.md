# Compute Log Randomized Predictive p-values using ISCV Method

Calculates log randomized predictive p-values (log-RPPs) based on the
Importance Sampling Cross-Validation (ISCV) method, using precomputed
log cumulative distribution function (log-CDF) and log probability mass
function (log-PMF) values.

## Usage

``` r
iscv_logrpp(log_cdf, log_pmf)
```

## Arguments

- log_cdf:

  A numeric matrix of log-CDF values with dimensions \\M \times N\\,
  where \\M\\ is the number of posterior draws and \\N\\ is the number
  of observations.

- log_pmf:

  A numeric matrix of log-PMF values with the same dimensions as
  `log_cdf`.

## Value

A numeric vector of length \\N\\, giving the log randomized predictive
p-values for each observation.

## Details

This function implements the Importance Sampling Cross-Validation (ISCV)
version of the log randomized predictive p-value (log-RPP) computation.
It uses numerically stable log-sum-exp operations to avoid overflow or
underflow when summing in log space.

Randomized predictive p-values (\\rpp\\) are computed using uniform
random draws \\u_i \sim \text{Uniform}(0,1)\\, combined with the
model-predicted \\\text{CDF}\\ and \\\text{PMF}\\ values for each
observation and posterior sample. The result is returned in log scale
for numerical stability.

Internally, this function:

- Generates uniform random values \\u_i\\ for each observation.

- Computes \\\log(u_i)\\ and adds it to \\\log(\text{PMF})\\.

- Applies numerically stable log-sum-exp operations column-wise.

Values numerically equal to 0 or 1 are replaced by \\\log(1e-5)\\ and
\\\log(9e-5)\\ respectively to keep calculations stable.

## See also

`colLogSumExps`, `log_sum_exp`, `loocv_logrpp`, `posterior_logrpp`

## Examples

``` r
# Example matrices (toy example)
log_cdf <- matrix(log(pnorm(rnorm(20))), 4000, 5)
log_pmf <- matrix(log(dnorm(rnorm(20))), 4000, 5)
result <- iscv_logrpp(log_cdf, log_pmf)
#> Error in iscv_logrpp(log_cdf, log_pmf): could not find function "iscv_logrpp"
head(result)
#> Error: object 'result' not found
```
