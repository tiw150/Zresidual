# Z-residuals for parametric survival regression models (survival)

S3 method for \[Zresidual()\] when the fitted model is a
\[survival::survreg()\] object (internally tagged as
\`"survreg.survival"\`). This is a thin wrapper around the existing
\`Zresidual.survreg()\` core implementation: it simply passes the fitted
object and optional \`data\` to the core function and then adds the
\`"zresid"\` class.

Users are expected to call \[Zresidual()\] directly rather than calling
\`Zresidual.survreg.survival()\`.

## Usage

``` r
# S3 method for class 'survreg.survival'
Zresidual(object, nrep = 1, data = NULL, ...)
```

## Arguments

- object:

  A fitted \[survival::survreg()\] model.

- nrep:

  Integer; number of randomized Z-residual replicates to generate.
  Defaults to \`1\`.

- data:

  Optional data frame containing the survival response and covariates;
  if \`NULL\`, the original model frame is used.

- ...:

  Currently ignored; included for method compatibility.

## Value

A numeric matrix of dimension \\n \times\\ `nrep`, with additional
attributes as produced by \`Zresidual_survreg_survival()\`. The returned
object is given class \`"zresid"\` in addition to any existing classes.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)
  fit_surv <- survreg(Surv(time, status) ~ age + sex,
                      data = lung, dist = "weibull")
  z_surv <- Zresidual(fit_surv, nrep = 5, data = lung)
} # }
```
