# A function to draw qq plot of a z-residual

A function to draw qq plot of a z-residual

## Usage

``` r
# S3 method for class 'zresid'
qqnorm(
  Zresidual,
  irep = 1,
  diagnosis.test = "SW",
  main.title = ifelse(is.null(attr(Zresidual, "type")), "Normal Q-Q Plot",
    paste("Normal Q-Q Plot -", attr(Zresidual, "type"))),
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles",
  outlier.return = TRUE,
  outlier.value = 3.5,
  outlier.set = list(),
  my.mar = c(5, 4, 4, 6) + 0.1,
  legend.settings = list(),
  ...
)
```

## Arguments

- Zresidual:

  A z-residual.

- diagnosis.test:

  diagnosis test

- legend.settings:

  A list of parameters available in legend().

- point.clr:

  A vector of colors of points.

- point.type:

  A vector of type of points.

- outlier.settings:

  A list of parameters available in symbols() and text().

- X.anova:

  X variable to run if ANOVA is selected in diagnosis.test. Possible
  c("lp", "covariate")

- k.anova:

  Number of bins if X.anova is 'covariates'.
