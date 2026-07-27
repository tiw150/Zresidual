# Plot martingale residuals for survival models

Draws a ggplot2 diagnostic plot against observation index, linear
predictor, or a selected covariate. The original LOWESS and external
`is.outlier` behavior are retained.

## Usage

``` r
# S3 method for class 'martg.resid'
plot(
  x,
  ylab = "Martingale Residual",
  x_axis_var = c("index", "lp", "covariate"),
  main.title = "Martingale Residual Plot",
  outlier.return = TRUE,
  point.args = list(),
  smooth.args = list(),
  theme = ggplot2::theme_bw(),
  ...
)
```

## Arguments

- x:

  A martingale residual object with `censored.status`, `linear.pred`,
  and `covariates` attributes as required by the selected x-axis.

- ylab:

  Y-axis label.

- x_axis_var:

  `"index"`, `"lp"`, `"covariate"`, or a covariate name.

- main.title:

  Plot title.

- outlier.return:

  Whether to use the calling environment's logical `is.outlier` vector
  and return its indices.

- point.args:

  Named arguments for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).

- smooth.args:

  Named arguments for the LOWESS `geom_line()`.

- theme:

  A ggplot2 theme.

- ...:

  Additional plotting arguments.

## Value

Invisibly returns `NULL`, or `list(outliers = ...)` when
`outlier.return = TRUE`, matching the original method.
