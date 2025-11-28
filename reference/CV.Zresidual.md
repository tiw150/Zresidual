# Cross-Validated Z-Residual Diagnostics for Survival Models

Computes cross-validated Z-residuals for fitted survival models,
including Cox proportional hazards models (\`coxph\`) with or without
frailty terms, and parametric survival regression models (\`survreg\`).
The function automatically detects the model type from the fitted model
object and applies the appropriate Z-residual cross-validation method.

## Usage

``` r
CV.Zresidual(fit.object, nfolds, foldlist = NULL, data = NULL, nrep = 1)
```

## Arguments

- fit.object:

  A fitted survival model object of class \`coxph\` or \`survreg\`.

- nfolds:

  Integer. Number of folds for cross-validation.

- foldlist:

  Optional list specifying custom fold assignments. If \`NULL\`, folds
  are generated internally.

- data:

  Optional dataset used to refit the model during cross-validation.
  Required when \`foldlist\` is provided or when the original model call
  does not contain the data explicitly.

- nrep:

  Integer. Number of repeated cross-validations to perform. Default is
  1.

## Value

An object of class \`"cvzresid"\` containing the cross-validated
Z-residual results and any model-specific diagnostic information.

## Details

The function identifies whether the fitted model is: - a Cox model
(\`coxph\`) with frailty terms, - a Cox model without frailty, - or a
parametric survival model (\`survreg\`), and dispatches to the
appropriate internal cross-validation function:
\`CV.Zresidual.coxph.frailty()\`, \`CV.Zresidual.coxph()\`, or
\`CV.Zresidual.survreg()\`.

All required packages are loaded via \`pacman::p_load()\`.

## See also

\`CV.Zresidual.coxph()\`, \`CV.Zresidual.coxph.frailty()\`,
\`CV.Zresidual.survreg()\`

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
out <- CV.Zresidual(fit, nfolds = 5)
} # }
```
