# Extract point-level quantities for predictive checks

Thin wrapper around
[`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md) that
requests the standardized point-level quantities required by predictive
checks.

## Usage

``` r
point_Zcov(fit, data, type = NULL, ndraws = NULL, draw_ids = NULL, ...)
```

## Arguments

- fit:

  A fitted model object.

- data:

  Evaluation data.

- type:

  Optional component selector.

- ndraws:

  Optional number of posterior draws.

- draw_ids:

  Optional posterior draw indices. When supplied, this is used to
  synchronize point quantities with `log_pointpred` and
  `postpred_simulate`.

- ...:

  Additional arguments passed to `Zcov`.

## Value

A standardized list containing at least `response`, `mean`, `variance`,
`lp`, and `covariate` when supported by the model backend.
