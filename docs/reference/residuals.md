# Residuals for supported survival models

High-level wrapper to compute residuals for supported survival models,
with a unified interface. Depending on the model type and presence of a
shared frailty term, this function dispatches to:

- [`residual.coxph()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.md)
  for standard Cox proportional hazards models,

- [`residual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.frailty.md)
  for shared frailty Cox models, and

- [`residual.survreg()`](https://tiw150.github.io/Zresidual/reference/residual.survreg.md)
  for parametric survival models fitted with
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

## Usage

``` r
residuals(
  fit.object,
  data,
  residual.type = c("censored Z-residual", "Cox-Snell", "martingale", "deviance")
)
```

## Arguments

- fit.object:

  A fitted model object. Currently, only
  [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) and
  [`survreg`](https://rdrr.io/pkg/survival/man/survreg.html) objects are
  supported. Cox models may optionally include a shared frailty term
  (e.g., `frailty(group)`).

- data:

  A `data.frame` containing the variables needed to evaluate residuals
  for `fit.object`. It must contain the survival response and all
  covariates (and, for shared frailty models, the grouping factor)
  appearing in the model formula.

- residual.type:

  Character string specifying the type of residual to compute. Valid
  values depend on the underlying model and residual function, but for
  Cox models typically include: `"censored Z-residual"`, `"Cox-Snell"`,
  `"martingale"`, and `"deviance"`. The default is the full vector
  `c("censored Z-residual", "Cox-Snell", "martingale", "deviance")`,
  which is typically resolved via
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html) in the
  underlying residual functions.

## Value

An object containing residuals, as returned by the corresponding
lower-level function:

- [`residual.coxph()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.md)
  or
  [`residual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.frailty.md)
  for Cox / shared frailty Cox models, or

- [`residual.survreg()`](https://tiw150.github.io/Zresidual/reference/residual.survreg.md)
  for parametric survival models.

The basic structure (numeric vector or matrix) and attributes are
preserved from the underlying function, but the returned object is
assigned an additional class `"zresid"` on top of its original classes.

## Details

This function is intended as a convenience entry point for users of the
package: it automatically routes to the appropriate residual computation
for the given model type and adds the class `"zresid"` to the result.
This extra class can be used by downstream plotting or diagnostic
functions (e.g.,
[`plot.zresid()`](https://tiw150.github.io/Zresidual/reference/plot.zresid.md)).

## See also

[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`survreg`](https://rdrr.io/pkg/survival/man/survreg.html),
`residual.coxph`, `residual.coxph.frailty`, `residual.survreg`

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  ## Cox PH model
  fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
  r1 <- residuals(fit_cox, data = lung,
                  residual.type = "censored Z-residual")

  ## Shared frailty Cox model
  lung$inst <- factor(lung$inst)
  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
                     data = lung)
  r2 <- residuals(fit_frail, data = lung,
                  residual.type = "Cox-Snell")

  ## Parametric survival model (Weibull)
  fit_weib <- survreg(Surv(time, status) ~ age + sex, data = lung,
                      dist = "weibull")
  r3 <- residuals(fit_weib, data = lung,
                  residual.type = "censored Z-residual")
} # }
```
