# Predictive quantities for survival::coxph with frailty

Predictive quantities for survival::coxph with frailty

## Usage

``` r
log_pointpred_survival_coxph.penal(fit, data, traindata, ...)
```

## Arguments

- fit:

  A fitted coxph object, usually inheriting from coxph.penal.

- data:

  New data to evaluate.

- traindata:

  Training data used to fit the model. This is required for
  reconstructing the baseline cumulative hazard.

- ...:

  Additional arguments passed from Zresidual.

## Value

A list containing log-survival probabilities, log-likelihood
contributions, discreteness indicators, and linear predictors.
