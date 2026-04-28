# ANOVA Test Core (internal)

Performs an ANOVA test for Z-residuals against a covariate. Numeric
covariates are binned; sparse bins (\<=2 obs) are removed.

## Usage

``` r
test.nl.aov(Zresidual, fitted.value, k.anova = 10)
```

## Arguments

- Zresidual:

  Numeric vector.

- fitted.value:

  Numeric or factor covariate.

- k.anova:

  Max bins for numeric covariates.
