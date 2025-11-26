# Z-residuals for Accelerated failure time model models

This function calculates Z-residuals based on a fitted \`survreg\`
object from the \`survival\` package.

## Usage

``` r
Zresidual.survreg(fit_survreg, newdata, n.rep = 1)
```

## Arguments

- fit_survreg:

  A fitted object from the \`survreg\` function in the \`survival\`
  package.

- newdata:

  Optional `data.frame` containing the variables used in
  `fit_coxph$formula`. If `NULL` (the default), residuals are computed
  on the original data used to fit the model. If supplied, `newdata`
  must contain the survival response and all covariates appearing in the
  original model formula.

- n.rep:

  An integer specifying the number of random draws to use for
  calculating the Z-residuals for censored observations. Defaults to
  \`nrep\` (which should be defined in your environment, a common choice
  is 100).

## Value

A numeric matrix of dimension \\n \times\\ `n.rep`, where \\n\\ is the
number of observations in the (new) model frame. Each column corresponds
to one set of Z-residuals. The returned matrix has the following
attributes attached:

- `Survival.Prob`: vector of survival probabilities \\S_i(t_i)\\.

- `linear.pred`: vector of linear predictors \\\eta_i\\.

- `covariates`: data frame of covariates (model frame without the
  response).

- `censored.status`: event indicator (1 = event,0 = censored).

- `object.model.frame`: the `model.frame` used to compute the residuals.

- `type`: character string `"survival"`.

## Examples

``` r
library(survival)

# Fit a Weibull survival regression model
fit_weibull <- survreg(Surv(time, status) ~ x, data = lung, dist = "weibull")
#> Error in eval(predvars, data, env): object 'x' not found

# Calculate Z-residuals for the fitted model
z_residuals <- Zresidual.survreg(fit_weibull)
#> Error in Zresidual.survreg(fit_weibull): argument "newdata" is missing, with no default
head(z_residuals)
#> Error: object 'z_residuals' not found

# Calculate Z-residuals for new data
new_data <- data.frame(x = c(1, 2, 0.5))
Z_residuals_new <- Zresidual.survreg(fit_weibull, newdata = new_data)
#> Error: object 'fit_weibull' not found
head(z_residuals_new)
#> Error: object 'z_residuals_new' not found

```
