# Calculate Z-residuals for a fitted survival regression model

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

  An optional data frame containing new observations for which to
  calculate the Z-residuals. If \`NULL\` (default), residuals are
  calculated for the data used to fit the model.

- n.rep:

  An integer specifying the number of random draws to use for
  calculating the Z-residuals for censored observations. Defaults to
  \`nrep\` (which should be defined in your environment, a common choice
  is 100).

## Value

- Survival.Prob: The estimated survival probabilities.

- linear.pred: The linear predictors from the survival regression model.

- covariates: The covariate values used in the model.

- censored.status: The censoring status (0 for censored, 1 for event).

- object.model.frame: The model frame used for the analysis.

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
