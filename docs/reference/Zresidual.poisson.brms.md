# Z-residuals for Poisson models fitted with brms

\`Zresidual.poisson.brms()\` is the S3 method for \[Zresidual()\] when
applied to Poisson models fitted with \[brms::brm()\] and \`family =
poisson()\`. Objects are dispatched here when the fitted object is a
\`"brmsfit"\` with \`brms::family(object)\$family == "poisson"\` and has
been internally tagged with the class \`"poisson.brms"\` by
\[Zresidual()\].

In normal use, users should call \[Zresidual()\] directly on the
\`brmsfit\` object (for example \`Zresidual(fit)\`), rather than calling
\`Zresidual.poisson.brms()\` explicitly.

## Usage

``` r
# S3 method for class 'poisson.brms'
Zresidual(object, nrep = 1, data, type, method = "iscv", ...)
```

## Arguments

- object:

  A \`brmsfit\` object with Poisson family
  (\`brms::family(object)\$family == "poisson"\`).

- nrep:

  Integer; the number of replicated Z-residual sets to generate. Default
  is \`1\`.

- data:

  Optional data frame used for prediction or residual computation. If
  `NULL` (default), the data stored inside the `brmsfit` object are
  used.

- type:

  Optional character string controlling the residual type, interpreted
  by the underlying implementation (if used).

- method:

  Character string specifying the residual calculation method:
  \`"iscv"\` for importance-sampled cross-validated randomized
  predictive p-values, \`"rpost"\` for posterior predictive p-values, or
  \`"mpost"\` for marginal posterior predictive p-values. Default is
  \`"iscv"\`.

- ...:

  Further arguments passed from \[Zresidual()\]. They are ignored by
  this method but are accepted for consistency with the generic.

## Value

A numeric matrix of Z-residuals (one column per replication) as returned
by \[Zresidual_poisson_brms()\], with the class \`"zresid"\` added to
its class vector.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(brms)
  fit_pois <- brm(y ~ x1 + x2, data = df,
                  family = poisson())

  ## ISCV-based Z-residuals
  z_pois <- Zresidual(fit_pois, method = "iscv")

  ## Posterior predictive Z-residuals with 2 replicates
  z_pois_post <- Zresidual(fit_pois, method = "rpost", nrep = 2)
} # }
```
