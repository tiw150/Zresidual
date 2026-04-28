# Boxplot diagnostics for Z-residuals

Produces boxplots of Z-residuals grouped by binned x-axis values. The
x-axis variable can be a linear predictor, a model covariate, or a
user-supplied vector. Optional metadata supplied through `info` (or
stored in attributes of `x`) are used to fill legacy plotting attributes
such as `"covariates"`, `"linear.pred"`, and `"type"`.

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

- x_axis_var:

  Variable used for grouping on the x-axis. It may be `"lp"`,
  `"covariate"`, a covariate name stored in `attr(x, "covariates")`, a
  length-\\n\\ vector, or a function returning such a vector.

- num.bin:

  Integer giving the number of bins used when the x-axis variable is
  numeric.

- normality.test:

  Character vector specifying which diagnostic p-values to display.
  Supported values are `"SW"`, `"AOV"`, and `"BL"`.

- k.test:

  Integer controlling grouping used by the diagnostic tests.

- main.title:

  Main title of the plot. If omitted, a default title is constructed
  from `attr(x, "type")`, when available.

- outlier.return:

  Logical; if `TRUE`, invisibly return the indices of observations with
  `|Z| > outlier.value`.

- outlier.value:

  Numeric threshold used to define outliers.

- ...:

  Additional graphical arguments passed to plotting functions.

## Value

Invisibly returns a list with component `outliers`, containing the
indices of observations flagged as outliers for the plotted replicate.
The main effect of the function is the boxplot.

## Details

This plot is useful for checking whether the distribution of Z-residuals
changes systematically across fitted values or covariates.

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
  z <- Zresidual(fit,data=dat, nrep = 1, seed = 1)
  info <- Zcov(fit, data = dat)

  boxplot(z, info = info, x_axis_var = "lp")
}
#> Warning: NaNs produced

```
