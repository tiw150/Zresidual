# Simulate posterior predictive draws for predictive checks

Generic bridge used by `predcheck_chisquare`. Methods should return an M
x N matrix, where rows are posterior draws or replications and columns
are observations.

## Usage

``` r
postpred_simulate(fit, data, type = NULL, ndraws = NULL, draw_ids = NULL, ...)
```

## Arguments

- fit:

  A fitted model object.

- data:

  Evaluation data.

- type:

  Optional component selector.

- ndraws:

  Optional number of draws.

- draw_ids:

  Optional posterior draw indices.

- ...:

  Additional arguments passed to methods.
