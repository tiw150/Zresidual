# Posterior predictive checks for count models

Computes posterior predictive check summaries for a fitted count model
using randomized and mid-point probability integral transform residuals,
a chi-squared discrepancy, and an optional ANOVA-style residual
diagnostic.

## Usage

``` r
ppc(
  fit,
  data,
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

- data:

  A data frame used for posterior predictive checking.

- predcheck_pointpred:

  Optional low-level backend function. This function must accept at
  least `fit`, `data`, and `draw_ids`, and return a list with components
  `support`, `family`, `y`, `n`, `ndraws`, `pmf`, `tail`, `rng`, and
  `moments`.

- x:

  Optional grouping variable for the ANOVA-style diagnostic. This can be
  either a column name in `data` or a vector of length `nrow(data)`.

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
res <- ppc(fit = fit_nb, data = dat, x = "depth", ndraws = 500, seed = 1)
print(res)
} # }
```
