# Z-residuals for Bernoulli models fitted with brms

`Zresidual.bernoulli.brms()` is the S3 method for
[`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)
when applied to Bernoulli (binary) regression models fitted with
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) and
`family = bernoulli()`. Objects are dispatched here when the fitted
object is a `"brmsfit"` with family `"bernoulli"` and has been
internally tagged with the class `"bernoulli.brms"` by
[`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

In most cases users should call
[`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)
directly on the `brmsfit` object, e.g. `Zresidual(fit)`, rather than
calling `Zresidual.bernoulli.brms()` explicitly.

## Usage

``` r
# S3 method for class 'bernoulli.brms'
Zresidual(object, nrep = 1, data = NULL, type = NULL, method = "iscv", ...)
```

## Arguments

- object:

  A `brmsfit` object with `brms::family(object)$family == "bernoulli"`.

- nrep:

  Integer; number of replicated Z-residual sets to generate. Defaults to
  `1`.

- data:

  Optional data frame used for prediction. If `NULL`, the data stored
  inside the `brmsfit` object are used.

- type:

  Optional character string controlling the residual type, interpreted
  by the underlying implementation (if used).

- method:

  Character string specifying the residual calculation method: `"iscv"`
  for importance-sampled cross-validated randomized predictive p-values,
  `"rpost"` for randomized posterior predictive p-values, or `"mpost"`
  for middle-value posterior predictive p-values. Default is `"iscv"`.

- ...:

  Further arguments passed to the underlying implementation

## Value

A numeric matrix of Z-residuals with one column per replication, as
returned by
[`Zresidual_bernoulli_brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual_bernoulli_brms.md),
but with the class `"zresid"` added to its class vector.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(brms)
  fit_bern <- brm(y ~ x1 + x2, data = df, family = bernoulli())

  ## ISCV-based Z-residuals
  z1 <- Zresidual(fit_bern, method = "iscv", nrep = 2)

  ## Posterior predictive Z-residuals
  z2 <- Zresidual(fit_bern, method = "rpost")
} # }
```
