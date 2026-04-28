# Extract aligned model metadata for Z-residual diagnostics

Extracts model-aligned response information, response type, covariates,
and linear predictors for use in Z-residual diagnostics and plotting
functions. The exact components returned depend on the fitted model
class and the value of `detail`.

## Usage

``` r
Zcov(fit, data, detail = c("basic", "full"), type = NULL, ...)

# S3 method for class 'brmsfit'
Zcov(fit, data, detail = c("basic", "full"), type = NULL, ...)

# S3 method for class 'coxph'
Zcov(fit, data, detail = c("basic", "full"), type = NULL, ...)

# S3 method for class 'survreg'
Zcov(fit, data, detail = c("basic", "full"), type = NULL, ...)

# S3 method for class 'glm'
Zcov(fit, data, detail = c("basic", "full"), type = NULL, ...)
```

## Arguments

- fit:

  A fitted model object, such as `brmsfit`,
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html), or
  [`survival::survreg`](https://rdrr.io/pkg/survival/man/survreg.html).

- data:

  Data used to extract aligned metadata. Must be provided.

- detail:

  Character string; either `"basic"` or `"full"`.

- type:

  Optional component selector. For example, in `brms` hurdle models, use
  `type = "zero"` to isolate the binary hurdle component,
  `type = "count"` for the positive-count component, and
  `type = "hurdle"` for the full hurdle model.

- ...:

  Additional arguments passed to class-specific methods.

## Value

A named list. Depending on the model and `detail`, this typically
includes aligned response information, response-type metadata,
covariates, linear predictors, and optional model-frame details.

## Methods (by class)

- `Zcov(brmsfit)`: Method for `brmsfit` objects. Currently, the
  top-level `type` argument is implemented **only** for hurdle or
  zero-inflated models. For these models, setting `type = "zero"`
  isolates the binary zero-inflation/hurdle process, `type = "count"`
  isolates the truncated count process, and `type = "hurdle"` evaluates
  the entire model. For all other standard families (e.g., standard
  logistic, Poisson, Gaussian), the `type` argument is unnecessary and
  will be ignored with a warning. For these standard models, the
  function simply returns the standard linear predictors, covariates,
  and response.

- `Zcov(coxph)`: Method for `coxph` objects from the survival package.
  Extracts the `Surv` response, linear predictors, and censoring
  metadata.

- `Zcov(survreg)`: Method for `survreg` objects from the survival
  package. Extracts the `Surv` response, linear predictors, distribution
  scale, and censoring metadata.

- `Zcov(glm)`: Method for `glm` objects. Standard generalized linear
  models are single-process models. Therefore, the `type` argument is
  unnecessary and will be ignored with a warning if provided. The
  function returns the standard linear predictors on the link scale,
  covariates, and response.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n <- 20
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.2 * x))
  t_cens  <- rexp(n, rate = 0.5)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)

  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
  info <- Zcov(fit, data = dat)
  names(info)
}
#>  [1] "type"          "family"        "response_name" "response"     
#>  [5] "covariates"    "linear_pred"   "obs_id"        "y_type"       
#>  [9] "y_type_kind"   "y_type_levels" "extra"        
```
