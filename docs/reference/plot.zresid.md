# Scatterplot diagnostics for Z-residuals

Produces scatterplots of Z-residuals against observation index, linear
predictors, model covariates, or a user-supplied x-axis variable.
Optional metadata supplied through `info` (or stored in attributes of
`x`) are used to fill legacy plotting attributes such as `"covariates"`,
`"linear.pred"`, `"zero_id"`, and `"censored.status"`.

## Usage

``` r
# S3 method for class 'zresid'
plot(
  x,
  zcov = NULL,
  info = NULL,
  irep = 1,
  ylab = "Z-Residual",
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  x_axis_var = "index",
  main.title = ifelse(is.null(attr(x, "type")), "Z-residual Scatterplot",
    paste("Z-residual Scatterplot -", attr(x, "type"))),
  outlier.return = TRUE,
  outlier.value = 3.5,
  category = NULL,
  outlier.set = list(),
  xlab = NULL,
  my.mar = c(5, 4, 4, 6) + 0.1,
  add_lowess = FALSE,
  ...
)
```

## Arguments

- x:

  A numeric matrix of Z-residuals, typically returned by
  [`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md),
  with one column per residual replicate.

- zcov:

  Optional metadata, typically returned by
  [`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md).

- info:

  Legacy alias for `zcov`.

- irep:

  Integer vector specifying which column(s) of `x` to plot.

- ylab:

  Label for the y-axis.

- normality.test:

  Character vector specifying which diagnostic p-values to display.
  Supported values are `"SW"`, `"AOV"`, and `"BL"`.

- k.test:

  Integer controlling grouping used by the diagnostic tests.

- x_axis_var:

  Variable used on the x-axis. It may be one of `"index"`, `"lp"`,
  `"covariate"`, a covariate name stored in `attr(x, "covariates")`, a
  length-\\n\\ vector, or a function returning such a vector.

- main.title:

  Main title of the plot. If omitted, a default title is constructed
  from `attr(x, "type")`, when available.

- outlier.return:

  Logical; if `TRUE`, mark observations with `|Z| > outlier.value` (and
  non-finite residuals) and invisibly return their indices.

- outlier.value:

  Numeric threshold used to define outliers.

- category:

  Optional grouping variable of length \\n\\ used to modify point
  appearance.

- outlier.set:

  A named list of graphical arguments passed to
  [`symbols`](https://rdrr.io/r/graphics/symbols.html) and
  [`text`](https://rdrr.io/r/graphics/text.html) when annotating
  outliers.

- xlab:

  Label for the x-axis. If `NULL`, an automatic label is used.

- my.mar:

  Numeric vector passed to
  [`par`](https://rdrr.io/r/graphics/par.html)`(mar = ...)`.

- add_lowess:

  Logical; if `TRUE`, add a LOWESS smooth when the x-axis is numeric.

- ...:

  Additional graphical arguments passed to plotting functions.

## Value

Invisibly returns a list with component `outliers`, containing the
indices of observations flagged as outliers for the plotted replicate.
The main effect of the function is the diagnostic scatterplot.

## Details

Depending on the metadata available, the plot can distinguish zero
versus positive observations for hurdle-type models, or censored versus
uncensored observations for survival models.

## See also

[`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md),
[`Zcov`](https://tiw150.github.io/Zresidual/reference/Zcov.md)

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  set.seed(1)
  n <- 30
  x <- rnorm(n)
  t_event <- rexp(n, rate = exp(0.3 * x))
  t_cens  <- rexp(n, rate = 0.4)
  status  <- as.integer(t_event <= t_cens)
  time    <- pmin(t_event, t_cens)
  dat <- data.frame(time = time, status = status, x = x)

  fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
  z <- Zresidual(fit, data=dat, nrep = 1, seed = 1)
  info <- Zcov(fit, data = dat)

  plot(z, info = info, x_axis_var = "index")
  plot(z, info = info, x_axis_var = "lp")
}
#> Warning: NaNs produced


```
