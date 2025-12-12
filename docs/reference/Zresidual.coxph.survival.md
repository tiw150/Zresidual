# Z-residuals for Cox proportional hazards models (survival package)

`Zresidual.coxph.survival()` computes randomized Z-residuals for Cox
proportional hazards models fitted with
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html), supporting both
standard and shared frailty models. This S3 method is designed to be
called via the generic function
[`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

The function automatically detects the presence of a frailty term (e.g.,
`frailty(group)`) and dispatches the calculation to one of two internal
implementations. These residuals are intended for **in-sample**
diagnostics and model assessment.

## Usage

``` r
# S3 method for class 'coxph.survival'
Zresidual(object, nrep = 1, data = NULL, type = NULL, method = NULL, ...)
```

## Arguments

- object:

  A fitted [`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model.
  Supports both standard Cox models and shared frailty models.

- nrep:

  Integer; number of independent randomized Z-residual replicates to
  generate. Defaults to `1`.

- data:

  Optional `data.frame` containing the survival response and covariates.
  When `NULL` (default), the residuals are computed on the data used to
  fit the `object`. This parameter is often aliased as `newdata` in the
  internal worker functions.

- type:

  Optional character string controlling the residual type. Set
  internally to `"survival"` for Cox models.

- method:

  Character string specifying the residual calculation method. Currently
  unused.

- ...:

  Further arguments passed to the underlying implementation functions.

## Value

A numeric matrix of class `"zresid"` with dimension \\N \times nrep\\.
Each column is an independent set of Z-residuals. The following
diagnostic attributes are attached:

- `Survival.Prob`: Vector of predicted survival probabilities
  \\S_i(t_i)\\.

- `linear.pred`: Vector of linear predictors \\\eta_i = \mathbf{x}\_i
  \mathbf{\hat{\beta}}\\.

- `covariates`: Data frame of covariates used in the model.

- `censored.status`: Event indicator (1 = event, 0 = censored).

- `object.model.frame`: The `model.frame` used for computation.

- `type`: Character string, always `"survival"`.

## Details

This method dispatches work based on the model formula:

- **Standard Cox models (No Frailty)**:

  The function calls
  [`Zresidual_coxph_survival`](https://tiw150.github.io/Zresidual/reference/Zresidual_coxph_survival.md)
  to compute Z-residuals using the fixed effects
  (\\\mathbf{x}\hat{\mathbf{\beta}}\\) and the estimated baseline
  cumulative hazard function \\\hat{H}\_0(t)\\.

- **Shared Frailty Cox models**:

  The function calls
  [`Zresidual_coxph_frailty_survival`](https://tiw150.github.io/Zresidual/reference/Zresidual_coxph_frailty_survival.md).
  This implementation computes residuals accounting for the
  cluster-level frailty effect (\\\hat{z}\_{\text{group}}\\). It
  requires the data used for fitting (`traindata`) to reconstruct the
  baseline hazard and estimate the frailty term.

**Randomization for Censored Observations:** Since the true survival
probability for a censored observation \$i\$ is known only to be greater
than \\S_i(t_i)\\, the Z-residual uses a randomized survival
probability: \\S\_{i, \text{rand}}(t_i) = S_i(t_i) \cdot U\\, where \\U
\sim \text{Unif}(0, 1)\\. This randomization is repeated `nrep` times.

## See also

[`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md),
[`Zresidual_coxph_survival`](https://tiw150.github.io/Zresidual/reference/Zresidual_coxph_survival.md),
[`Zresidual_coxph_frailty_survival`](https://tiw150.github.io/Zresidual/reference/Zresidual_coxph_frailty_survival.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(survival)

  ## Standard Cox model (no frailty term)
  fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
  # Note: The internal class 'coxph.survival' is usually added by a wrapper,
  # but Zresidual() handles dispatch automatically.
  z_cox <- Zresidual(fit_cox, nrep = 10, data = lung)

  ## Shared frailty Cox model (in-sample residuals)
  # Note: 'inst' must be a grouping factor.
  lung$inst_f <- factor(lung$inst)
  fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst_f),
                     data = lung)
  z_in <- Zresidual(fit_frail, nrep = 5)
} # }
```
