# Predictive quantities for survival::coxph with frailty

Predictive quantities for survival::coxph with frailty

## Usage

``` r
log_pointpred_survival_coxph.penal(fit, data, traindata, ...)
```

## Arguments

- fit:

  A fitted coxph object (usually class coxph.penal).

- data:

  The new data to evaluate.

- traindata:

  The original data used to fit the model (required for baseline
  hazard).

- ...:

  Additional arguments passed from Zresidual.
