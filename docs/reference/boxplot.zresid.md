# Boxplot diagnostics for Z-residuals

Produces ggplot2 boxplots of Z-residuals grouped by binned x-axis
values. The original binning rules and SW, AOV, and Bartlett diagnostics
are retained.

## Usage

``` r
# S3 method for class 'zresid'
boxplot(
  x,
  zcov = NULL,
  info = NULL,
  irep = 1,
  x_axis_var = "lp",
  num.bin = 10,
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  main.title = ifelse(is.null(attr(x, "type")), "Z-residual Boxplot",
    paste("Z-residual Boxplot -", attr(x, "type"))),
  outlier.return = FALSE,
  outlier.value = 3.5,
  interactive = FALSE,
  theme = ggplot2::theme_bw(),
  boxplot.args = list(),
  ...
)
```

## Arguments

- x:

  A numeric vector or matrix of Z-residuals, with one column per
  residual replication.

- zcov:

  Optional metadata returned by
  [`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- info:

  Legacy alias for `zcov`.

- irep:

  Integer vector selecting residual columns.

- x_axis_var:

  `"lp"`, `"covariate"`, a covariate name, a length-n vector, or a
  function returning such a vector.

- num.bin:

  Number of bins for a numeric x-axis variable.

- normality.test:

  Any of `"SW"`, `"AOV"`, and `"BL"`.

- k.test:

  Grouping parameter passed to the diagnostic functions.

- main.title:

  Main plot title.

- outlier.return:

  Whether to report and return observations satisfying the absolute
  Z-residual outlier rule.

- outlier.value:

  Absolute Z-residual outlier threshold.

- interactive:

  Convert each plot with
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).

- theme:

  A ggplot2 theme.

- boxplot.args:

  Named list passed to
  [`ggplot2::geom_boxplot()`](https://ggplot2.tidyverse.org/reference/geom_boxplot.html).

- ...:

  Additional named boxplot aesthetics. Common base names `col`, `lwd`,
  and `cex` are translated where possible.

## Value

Invisibly returns a list containing `outliers` and `plots`. Multiple
selected replications are printed sequentially.
