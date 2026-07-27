# Cross-validated Z-residuals for Cox models

Internal Cox backend used by
[`Zresidual_CV()`](https://tiw150.github.io/Zresidual/reference/Zresidual_CV.md).
The function supports ordinary `coxph` models and `coxph.penal` frailty
models. Each fold is refitted on its training observations and evaluated
on its held-out observations.

## Usage

``` r
Zresidual_CV_survival_coxph(
  object,
  data,
  nfolds = NULL,
  foldlist = NULL,
  nrep = 1L,
  randomized = TRUE,
  type = NULL,
  log_pointpred = NULL,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A fitted
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html)
  object.

- data:

  Data containing all variables used in `object`.

- nfolds:

  Number of folds when `foldlist` is not supplied.

- foldlist:

  Optional list of held-out indices relative to the internally retained
  complete-case data.

- nrep:

  Number of randomized residual replicates.

- randomized:

  Logical; generate randomized residuals for censored observations when
  `TRUE`.

- type:

  Optional model subtype passed to
  [`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

- log_pointpred:

  Optional predictive backend function or function name.

- seed:

  Optional integer seed.

- ...:

  Additional arguments passed to
  [`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

## Value

A `cvzresid` matrix.

## Details

For frailty models, cross-validation is currently within-group: every
group represented in a test fold must also be represented in that fold's
training data.
