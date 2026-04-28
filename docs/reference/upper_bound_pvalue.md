# Distribution-free upper-bound p-value from replicated p-values

Combines replicated p-values into a single valid upper-bound p-value
using an order-statistic correction.

## Usage

``` r
upper_bound_pvalue(p_rep, r = ceiling(length(p_rep)/2), H = NULL, ...)
```

## Arguments

- p_rep:

  Numeric vector of replicated p-values in `[0, 1]`.

- r:

  Integer order-statistic rank. The default is
  `ceiling(length(p_rep) / 2)`.

- H:

  Optional precomputed correction constant.

- ...:

  Additional arguments passed to the internal computation of
  \\H\_{J,r}\\ when `H = NULL`.

## Value

A single numeric p-value in `[0, 1]`.

## Examples

``` r
p <- c(0.03, 0.08, 0.12, 0.05, 0.09)
upper_bound_pvalue(p)
#> [1] 0.09582172
upper_bound_pvalue(p, r = 2)
#> [1] 0.08294784
```
