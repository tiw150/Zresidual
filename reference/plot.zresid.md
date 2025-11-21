# A function to draw scatter plot of a Z-residual

A function to draw scatter plot of a Z-residual

## Usage

``` r
# S3 method for class 'zresid'
plot(
  Zresidual,
  irep = 1:ncol(Zresidual),
  ylab = "Z-Residual",
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  X = c("index", "covariate", "lp"),
  main.title = ifelse(is.null(attr(Zresidual, "type")), "Z-residual Scatterplot",
    paste("Z-residual Scatterplot -", attr(Zresidual, "type"))),
  outlier.return = TRUE,
  outlier.value = 3.5,
  category = NULL,
  outlier.set = list(),
  xlab = NULL,
  my.mar = c(5, 4, 4, 6) + 0.1,
  ...
)
```

## Arguments

- Zresidual:

  A Z-residual.

- outlier.set:

  A list of parameters available in symbols() and text().

- xlab:

  Optional label for the x-axis. If NULL (default), the function
  determines the label based on 'X'. Accepts LaTeX code enclosed in
  'tex()' for enhanced formatting using the latex2exp package.
