# Scatterplot diagnostics for Z-residuals

Draws ggplot2 scatterplots of Z-residuals against observation index, a
linear predictor, a model covariate, or a user-supplied x-axis variable.
Multiple residual replications are displayed in facets.

## Usage

``` r
# S3 method for class 'zresid'
plot(
  x,
  zcov = NULL,
  irep = 1L,
  ylab = "Z-Residual",
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10L,
  x_axis_var = "index",
  main.title = NULL,
  outlier.return = TRUE,
  outlier.value = 3.5,
  category = NULL,
  outlier.set = list(),
  xlab = NULL,
  my.mar = c(5, 4, 4, 6) + 0.1,
  add_lowess = FALSE,
  reference.lines = c(1.96, 3),
  facet = TRUE,
  x_scale = c("identity", "log10"),
  interactive = FALSE,
  point.args = list(),
  smooth.args = list(),
  theme = ggplot2::theme_bw(),
  ...
)
```

## Arguments

- x:

  A numeric vector or matrix of Z-residuals. Matrix columns are residual
  replications.

- zcov:

  Optional metadata returned by
  [`Zcov()`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- irep:

  Integer vector selecting residual columns.

- ylab:

  Label for the y-axis.

- normality.test:

  Any of `"SW"`, `"AOV"`, and `"BL"`.

- k.test:

  Number of groups used by AOV and Bartlett diagnostics.

- x_axis_var:

  `"index"`, `"lp"`, `"covariate"`, a covariate name, a length-n vector,
  or a function of `x` and `zcov`.

- main.title:

  Plot title.

- outlier.return:

  Whether to mark and label outliers.

- outlier.value:

  Absolute residual threshold for outliers.

- category:

  Optional grouping variable of length n.

- outlier.set:

  Named list controlling outlier colour, size, label size, and whether
  labels are shown. Supported names are `colour`, `size`,
  `label_colour`, `label_size`, and `label`.

- xlab:

  Optional x-axis label.

- my.mar:

  Retained for source compatibility; margins are controlled by the
  ggplot theme and this argument is ignored.

- add_lowess:

  Whether to add the original
  [`stats::lowess()`](https://rdrr.io/r/stats/lowess.html) smooth.

- reference.lines:

  Numeric horizontal reference lines. Their negatives are also drawn.

- facet:

  Whether multiple selected replications are shown in facets. When
  `FALSE`, replications are overlaid and distinguished by linetype only
  for the LOWESS curves; faceting is recommended.

- x_scale:

  Either `"identity"` or `"log10"`.

- interactive:

  If `TRUE`, convert the ggplot with
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).

- point.args:

  Named list of arguments passed to `geom_point()`.

- smooth.args:

  Named list of arguments passed to the LOWESS `geom_line()` layer.

- theme:

  A complete or partial ggplot2 theme.

- ...:

  Additional named arguments for the main `geom_point()` layer.

## Value

A ggplot object, invisibly. With `interactive = TRUE`, returns a plotly
htmlwidget. Diagnostics and outlier indices are stored in the
`zresid_diagnostics` and `zresid_outliers` attributes of the ggplot.
