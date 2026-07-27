# Predictive quantities for survival::coxph with frailty

Compute point-wise predictive quantities for a Cox proportional hazards
model with a shared frailty term fitted by survival::coxph.

## Usage

``` r
log_pointpred_survival_coxph_frailty(fit, data, traindata, ...)
```

## Arguments

- fit:

  A fitted survival::coxph object with a frailty term.

- data:

  The new data to evaluate.

- traindata:

  The training data used to fit the model. This is required for
  constructing the fold-specific baseline cumulative hazard in
  cross-validation.

- ...:

  Additional arguments passed from Zresidual().

## Value

A named list containing log_surv, log_like, and is_discrete.
