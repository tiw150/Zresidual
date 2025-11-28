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
data(cancer, package="survival")
fit_weibull <- survreg(Surv(time, status) ~ age, data = lung, dist = "weibull")

# Calculate Z-residuals for the fitted model
z_residuals <- Zresidual.survreg(fit_weibull,newdata=lung)
head(z_residuals)
#>      Z-residual 1
#> [1,]    0.1427079
#> [2,]    0.5487250
#> [3,]    2.7479618
#> [4,]   -0.5074260
#> [5,]    1.4224828
#> [6,]    2.3325122

```
