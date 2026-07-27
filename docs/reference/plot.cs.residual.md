# Cox-Snell residual plot for survival models

Draws a ggplot2 Cox-Snell residual diagnostic based on the
Fleming-Harrington estimate of the cumulative hazard. Under a correctly
specified model, the curve should be close to the 45-degree reference
line.

## Usage

``` r
# S3 method for class 'cs.residual'
plot(
  x,
  ylab = "Cumulative Hazard Function",
  main.title = "Cox-Snell Residuals Scatterplot",
  outlier.return = TRUE,
  point.args = list(),
  curve.args = list(),
  reference.args = list(),
  theme = ggplot2::theme_bw(),
  ...
)
```

## Arguments

- x:

  Numeric Cox-Snell residual vector or one-column matrix carrying a
  length-n `censored.status` event-indicator attribute.

- ylab:

  Y-axis label.

- main.title:

  Main plot title.

- outlier.return:

  If `TRUE`, mark and return observations satisfying the original
  `abs(x) > 3.5` rule.

- point.args:

  Named arguments for the KM-point `geom_point()` layer.

- curve.args:

  Named arguments for the cumulative-hazard `geom_step()`.

- reference.args:

  Named arguments for the `y = x` reference segment.

- theme:

  A ggplot2 theme.

- ...:

  Additional plotting arguments. Common base arguments `xlab`, `ylab`,
  `main`, `xlim`, `ylim`, `col`, `pch`, and `cex` are supported.

## Value

Invisibly returns `NULL`, or `list(outliers = ...)` when
`outlier.return = TRUE`, matching the original method.
