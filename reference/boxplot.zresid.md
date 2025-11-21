# A function to draw box plot of a z-residual

A function to draw box plot of a z-residual

## Usage

``` r
# S3 method for class 'zresid'
boxplot(
  Zresidual,
  irep = 1,
  X = c("lp", "covariate"),
  num.bin = 10,
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  main.title = paste("Z-residual Boxplot -", attr(Zresidual, "type")),
  outlier.return = FALSE,
  outlier.value = 3.5,
  ...
)
```

## Arguments

- Zresidual:

  A z-residual.
