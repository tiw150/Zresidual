# Cross-validated Z-residuals

Generic function for computing cross-validated Z-residuals.

## Usage

``` r
CV.Zresidual(
  object,
  nfolds = NULL,
  foldlist = NULL,
  data = NULL,
  nrep = 1,
  log_pointpred = NULL,
  mcmc_summarize = c("post", "iscv"),
  type = NULL,
  randomized = TRUE,
  seed = NULL,
  ...
)

# S3 method for class 'survreg'
CV.Zresidual(
  object,
  nfolds = NULL,
  foldlist = NULL,
  data = NULL,
  nrep = 1,
  log_pointpred = NULL,
  mcmc_summarize = c("post", "iscv"),
  type = NULL,
  randomized = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model object.

- nfolds:

  Integer number of folds. If `NULL` and `foldlist` is also `NULL`, a
  default value is chosen.

- foldlist:

  Optional custom list of test-set indices, one element per fold.

- data:

  Data used for cross-validation refitting. Must be supplied.

- nrep:

  Integer number of Z-residual replicates.

- log_pointpred:

  Optional function or function name passed to
  [`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

- mcmc_summarize:

  Posterior summarization method for Bayesian fits.

- type:

  Optional component selector used by
  [`Zresidual()`](https://tiw150.github.io/Zresidual/reference/Zresidual.md)
  and [`Zcov()`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- randomized:

  Logical; whether to generate randomized Z-residuals.

- seed:

  Optional integer seed.

- ...:

  Further arguments passed to
  [`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md).

## Value

A numeric matrix of class `"cvzresid"` and `"zresid"`, with one row per
observation and one column per replicate.
