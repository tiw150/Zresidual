# Z-residuals for hurdle negative binomial models fitted with brms

\`Zresidual.hurdle_negbinomial.brms()\` is the S3 method for
\[Zresidual()\] when applied to hurdle negative binomial models fitted
with \[brms::brm()\] and \`family = hurdle_negbinomial()\`. Objects are
dispatched here when the fitted object is a \`"brmsfit"\` with
\`brms::family(object)\$family == "hurdle_negbinomial"\` and has been
internally tagged with the class \`"hurdle_negbinomial.brms"\` by
\[Zresidual()\].

In normal use, users should call \[Zresidual()\] directly on the
\`brmsfit\` object (for example \`Zresidual(fit)\`), rather than calling
\`Zresidual.hurdle_negbinomial.brms()\` explicitly.

## Usage

``` r
# S3 method for class 'hurdle_negbinomial.brms'
Zresidual(
  object,
  nrep = 1,
  data,
  type = c("hurdle", "count", "zero"),
  method = "iscv",
  ...
)
```

## Arguments

- object:

  A \`brmsfit\` object with hurdle negative binomial family
  (\`brms::family(object)\$family == "hurdle_negbinomial"\`).

- nrep:

  Integer; number of replicated Z-residual sets to generate. Default is
  \`1\`.

- data:

  Optional data frame used for prediction or residual computation. If
  `NULL`, the data stored inside the `brmsfit` object are used.

- type:

  Character string specifying which part of the model to compute
  Z-residuals for:

  - \`"zero"\` — the hurdle/zero part;

  - \`"count"\` — the truncated negative binomial count part;

  - \`"hurdle"\` — the full hurdle negative binomial model.

  The default is \`"hurdle"\`.

- method:

  Character string specifying the residual calculation method:
  \`"iscv"\` for importance-sampled cross-validated randomized
  predictive p-values, \`"rpost"\` for randomized posterior predictive
  p-values, or \`"mpost"\` for middle-value posterior predictive
  p-values. Default is \`"iscv"\`.

- ...:

  Further arguments passed to the underlying implementation function
  \[Zresidual_hurdle_negbinomial_brms()\].

## Value

A numeric matrix of Z-residuals (one column per replication) as returned
by \[Zresidual_hurdle_negbinomial_brms()\], with the class \`"zresid"\`
added to its class vector.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(brms)
  fit_hnb <- brm(y ~ x1 + x2, data = df,
                 family = hurdle_negbinomial())

  ## Counts part only
  z_count <- Zresidual(fit_hnb, type = "count", method = "iscv")

  ## Full hurdle model with 2 replicates
  z_hurdle <- Zresidual(fit_hnb, type = "hurdle", nrep = 2)
} # }
```
