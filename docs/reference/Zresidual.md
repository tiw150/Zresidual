# Compute Z-residuals via an S3 generic

`Zresidual()` is an S3 generic for computing Z-residuals for a range of
survival and regression models. It currently supports Cox proportional
hazards models (`coxph`) with or without frailty terms, parametric
survival regression models (`survreg`), and Bayesian regression models
fitted with **brms** (`brmsfit`) for several common count and hurdle
families.

The function inspects the fitted object, assigns an internal model– and
package–specific class (e.g. `"coxph.survival"`, `"survreg.survival"`,
`"poisson.brms"`, `"hurdle_poisson.brms"`), and dispatches to the
corresponding `Zresidual.*()` method.

## Usage

``` r
Zresidual(object, nrep = 1, data = NULL, type = NULL, method = "iscv", ...)
```

## Arguments

- object:

  A fitted model object. Currently supported objects include:

  - `"coxph"` – Cox proportional hazards model (from **survival**)

  - `"survreg"` – Parametric survival regression (from **survival**)

  - `"brmsfit"` – Bayesian regression model (from **brms**) with one of
    the supported families listed under **Details**.

- nrep:

  Integer. Number of replicated Z-residual samples to compute. Default
  is `1`.

- data:

  Optional data frame used for generating model predictions and
  residuals. For some models (e.g. `coxph` with frailty terms) this
  should match the data used to fit the model.

- type:

  Optional character string. Residual type for hurdle and count models
  (mainly used for **brms** and other Bayesian/count-model interfaces).
  The interpretation of this argument is model-specific.

- method:

  Character string indicating which predictive p-value scheme to use
  when computing Z-residuals for Bayesian or simulation-based models.
  Common options include `"iscv"` (importance-sampling cross-validated),
  `"loocv"`, and `"posterior"`. The default is `"iscv"`.

- ...:

  Further arguments passed on to model-specific methods such as
  [`Zresidual.coxph.survival()`](https://tiw150.github.io/Zresidual/reference/Zresidual.coxph.survival.md),
  [`Zresidual.survreg.survival()`](https://tiw150.github.io/Zresidual/reference/Zresidual.survreg.survival.md),
  [`Zresidual.poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.poisson.brms.md),
  etc.

## Value

An object of class `"zresid"` with additional class information
inherited from the model-specific method. The object contains the
computed Z-residuals and any model-specific diagnostic quantities
returned by the underlying `Zresidual.*()` implementation.

## Details

Internally, `Zresidual()` modifies the class of `object` to encode both
the model type and the fitting package, and then uses S3 method
dispatch:

- For **survival** models

  - `coxph` objects receive the class `"coxph.survival"` and are handled
    by
    [`Zresidual.coxph.survival()`](https://tiw150.github.io/Zresidual/reference/Zresidual.coxph.survival.md),
    which further distinguishes between models with and without frailty
    terms.

  - `survreg` objects receive the class `"survreg.survival"` and are
    handled by
    [`Zresidual.survreg.survival()`](https://tiw150.github.io/Zresidual/reference/Zresidual.survreg.survival.md).

- For **brms** models

  - The family is obtained via `family(object)$family`.

  - Currently supported families include `"hurdle_negbinomial"`,
    `"hurdle_poisson"`, `"negbinomial"`, `"poisson"`, and `"bernoulli"`.

  - Each supported family is mapped to a class of the form
    `"<family>.brms"`, and dispatched to a method such as
    [`Zresidual.hurdle_poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.hurdle_poisson.brms.md),
    [`Zresidual.negbinomial.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.negbinomial.brms.md),
    [`Zresidual.poisson.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.poisson.brms.md),
    or
    [`Zresidual.bernoulli.brms()`](https://tiw150.github.io/Zresidual/reference/Zresidual.bernoulli.brms.md).

Additional model classes (e.g. `glmmTMB` fits) can be supported by
defining new S3 methods such as `Zresidual.poisson.glmmTMB()` and
ensuring that `Zresidual()` assigns the corresponding internal class
(for example `"poisson.glmmTMB"`).

If `object` does not match any supported model class, or if a `brmsfit`
object uses an unsupported family, `Zresidual()` raises an error.

## See also

Zresidual.coxph.survival(), Zresidual.survreg.survival(),
Zresidual.hurdle_poisson.brms(), Zresidual.hurdle_negbinomial.brms(),
Zresidual.negbinomial.brms(), Zresidual.poisson.brms(),
Zresidual.bernoulli.brms(), and related model-specific methods.

## Examples

``` r
if (FALSE) { # \dontrun{
## Cox proportional hazards model
library(survival)
fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
z_cox <- Zresidual(fit_cox, nrep = 10, data = lung)

## Parametric survival regression
fit_surv <- survreg(Surv(time, status) ~ age + sex,
                    data = lung, dist = "weibull")
z_surv <- Zresidual(fit_surv, nrep = 5, data = lung)

## Bayesian Poisson regression with brms
library(brms)
fit_brms <- brm(count ~ x, data = df, family = poisson)
z_brms <- Zresidual(fit_brms, method = "posterior")
} # }
```
