# Normal Q-Q plot for Z-residuals

Produces a ggplot2 normal Q-Q plot for one or more columns of a
`"zresid"` object. The Shapiro-Wilk diagnostic, outlier annotation,
extreme-value compression, reference lines, and axis-break marks retain
the behavior of the original plotting method.

## Usage

``` r
# S3 method for class 'zresid'
qqnorm(
  y,
  zcov = NULL,
  info = NULL,
  irep = 1,
  diagnosis.test = "SW",
  main.title = ifelse(is.null(attr(y, "type")), "Normal Q-Q Plot",
    paste("Normal Q-Q Plot -", attr(y, "type"))),
  xlab = "Theoretical Quantiles",
  ylab = "Sample Quantiles",
  outlier.return = TRUE,
  outlier.value = 3.5,
  outlier.set = list(),
  outlier.label.xpad = 0.5,
  my.mar = c(5, 4, 4, 6) + 0.1,
  legend.settings = list(),
  clip.extreme = TRUE,
  clip.threshold = 6,
  interactive = FALSE,
  theme = ggplot2::theme_bw(),
  ...
)
```

## Arguments

- y:

  A numeric matrix of Z-residuals, typically returned by
  [`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md),
  with one column per residual replicate.

- zcov:

  Optional metadata, typically returned by
  [`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- info:

  Legacy alias for `zcov`.

- irep:

  Integer vector specifying which column(s) of `y` to plot.

- diagnosis.test:

  Character string specifying the normality diagnostic. Currently `"SW"`
  is supported.

- main.title:

  Main title of the plot.

- xlab:

  Label for the x-axis.

- ylab:

  Label for the y-axis.

- outlier.return:

  Logical; mark and return observations satisfying the outlier rule when
  `TRUE`.

- outlier.value:

  Numeric absolute-residual outlier threshold.

- outlier.set:

  Named list controlling outlier annotation. Supported ggplot settings
  are `colour`, `size`, `label_colour`, `label_size`, and `label`.

- outlier.label.xpad:

  Numeric padding added to the right x-axis limit.

- my.mar:

  Retained for source compatibility; ggplot themes control plot margins.

- legend.settings:

  Named list modifying the line legend. Supported names include
  `position`, `text_size`, and `title_size`.

- clip.extreme:

  Logical; visually compress extreme residuals.

- clip.threshold:

  Numeric clipping threshold.

- interactive:

  Logical; convert each ggplot using
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html)
  when `TRUE`.

- theme:

  A ggplot2 theme.

- ...:

  Additional named arguments passed to the QQ-point `geom_point()`
  layer.

## Value

Invisibly returns a list containing `outliers` and `plots`. A single
selected replication produces one ggplot in `plots`; multiple selected
replications are printed sequentially, matching the original method.
