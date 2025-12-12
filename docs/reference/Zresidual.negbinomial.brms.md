# Z-residuals for negative binomial models fitted with brms

\`Zresidual.negbinomial.brms()\` is the S3 method for \[Zresidual()\]
when applied to negative binomial models fitted with \[brms::brm()\] and
\`family = negbinomial()\`. Objects are dispatched here when the fitted
object is a \`"brmsfit"\` with \`brms::family(object)\$family ==
"negbinomial"\` and has been internally tagged with the class
\`"negbinomial.brms"\` by \[Zresidual()\].

In normal use, users should call \[Zresidual()\] directly on the
\`brmsfit\` object (for example \`Zresidual(fit)\`), rather than calling
\`Zresidual.negbinomial.brms()\` explicitly.

## Usage

``` r
# S3 method for class 'negbinomial.brms'
Zresidual(object, nrep = 1, data, type, method = "iscv", ...)
```

## Arguments

- object:

  A \`brmsfit\` object with negative binomial family
  (\`brms::family(object)\$family == "negbinomial"\`).

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

  Further arguments passed to the underlying implementation function
  \[Zresidual_negbinomial_brms()\].

## Value

A numeric matrix of Z-residuals (one column per replication) as returned
by \[Zresidual_negbinomial_brms()\], with the class \`"zresid"\` added
to its class vector.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(brms)
  fit_nb <- brm(y ~ x1 + x2, data = df,
                family = negbinomial())

  ## ISCV-based Z-residuals
  z_nb <- Zresidual(fit_nb, method = "iscv")

  ## Posterior predictive Z-residuals with 2 replicates
  z_nb_post <- Zresidual(fit_nb, method = "rpost", nrep = 2)
} # }
```
