# Z-residuals for Cox proportional hazards models (survival package)

\`Zresidual.coxph.survival()\` computes randomized Z-residuals for Cox
proportional hazards models fitted with \[survival::coxph()\], with or
without shared frailty terms. It is the S3 method for \[Zresidual()\]
when the fitted object is internally tagged with the class
\`"coxph.survival"\`, and is normally called via \[Zresidual()\] rather
than directly.

The function automatically detects whether the model formula contains a
frailty term (e.g. \`frailty(group)\`) and dispatches to separate
internal implementations for standard Cox models and shared frailty Cox
models. In both cases, Z-residuals are computed for a single data set
(typically the data used to fit the model) and are intended for
in-sample diagnostics.

## Usage

``` r
# S3 method for class 'coxph.survival'
Zresidual(object, nrep = 1, data = NULL, ...)
```

## Arguments

- object:

  A fitted \[survival::coxph()\] model. The function supports both
  standard Cox models and shared frailty Cox models specified with a
  term such as `frailty(group)` in the formula.

- nrep:

  Integer; number of independent randomized Z-residual replicates to
  generate. Defaults to `1`. Each replicate corresponds to a different
  randomization of censored observations.

- data:

  Optional `data.frame` containing the survival response and covariates
  used in `object$terms`. When `NULL` (default), the model frame is
  reconstructed from `object` and residuals are computed on the original
  data.

- ...:

  Further arguments passed to the underlying implementation functions.
  Currently unused.

## Value

A numeric matrix of dimension \\n \times\\ `nrep`, where \\n\\ is the
number of observations in the data set on which residuals are evaluated
(either `data` if supplied, or the original model frame). Each column
corresponds to one set of Z-residuals. The returned matrix has several
attributes attached, including:

- `Survival.Prob`: vector of survival probabilities \\S_i(t_i)\\ (or
  \\S\_{ij}(t_i)\\ in the frailty case).

- `linear.pred`: vector of linear predictors \\\eta_i\\ (fixed effects;
  excluding the frailty term in shared frailty models).

- `covariates`: data frame of covariates (model frame without the
  survival response and grouping factor).

- `censored.status`: event indicator (1 = event, 0 = censored).

- `object.model.frame`: the `model.frame` used to compute the residuals.

- `type`: character string, typically `"survival"`.

## Details

There are two main usage patterns:

- **Standard Cox models (no frailty term)**:

  When the model formula does not contain a frailty term, the method
  computes Z-residuals for a fitted \`coxph\` object using the survival
  response and covariates in `data`. If `data` is `NULL`, the model
  frame is reconstructed from the fitted object and residuals are
  computed on the original data used to fit the model.

  Internally, this branch delegates to an implementation function such
  as `Zresidual_coxph_survival()`.

- **Shared frailty Cox models**:

  If the model includes a multiplicative frailty term (e.g.
  `frailty(group)`), the method computes randomized Z-residuals that
  account for the cluster-level frailty. The same data set (either
  `data` if supplied, or the original model frame) is used both to
  reconstruct the baseline hazard / frailty effects and to evaluate
  residuals.

  Internally, this branch calls a dedicated frailty implementation, such
  as `Zresidual_coxph_frailty_survival()`.

## See also

Zresidual(), Zresidual_coxph_survival(),
Zresidual_coxph_frailty_survival()

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  ## Standard Cox model (no frailty term)
  fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
  z_cox <- Zresidual(fit_cox, nrep = 10, data = lung)

  ## Shared frailty Cox model (in-sample residuals)
  lung$inst <- factor(lung$inst)
  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
                     data = lung)
  z_in <- Zresidual(fit_frail, nrep = 5)
} # }
```
