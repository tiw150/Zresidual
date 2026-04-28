# Holdout predictive checks for count models

Computes holdout predictive check summaries for a fitted count model on
a new dataset. The returned summaries use the same discrepancy measures
as [`ppc()`](https://tiw150.github.io/Zresidual/reference/ppc.md), but
are evaluated on `newdata`.

## Usage

``` r
hpc(
  fit,
  newdata,
  predcheck_pointpred = NULL,
  x = NULL,
  ndraws = NULL,
  seed = NULL,
  k_anova = 10,
  ...
)
```

## Arguments

- fit:

  A fitted model object.

- newdata:

  A data frame used for holdout predictive checking.

- predcheck_pointpred:

  Optional low-level backend function. This function must accept at
  least `fit`, `data`, and `draw_ids`, and return a list with components
  `support`, `family`, `y`, `n`, `ndraws`, `pmf`, `tail`, `rng`, and
  `moments`.

- x:

  Optional grouping variable for the ANOVA-style diagnostic. This can be
  either a column name in `newdata` or a vector of length
  `nrow(newdata)`.

- ndraws:

  Optional number of posterior draws to use. If `NULL`, all available
  draws are used.

- seed:

  Optional random seed used for draw subsampling and randomized residual
  generation.

- k_anova:

  Maximum number of bins used when `x` is numeric.

- ...:

  Additional arguments passed to `predcheck_pointpred`.

## Value

An object of class `"predcheck"` containing the predictive-check summary
statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- hpc(fit = fit_nb, newdata = dat_hold, x = "depth", ndraws = 500, seed = 1)
print(res)
} # }
```
