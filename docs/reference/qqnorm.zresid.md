# Normal Q-Q plot for Z-residuals

Produces a normal Q-Q plot for one or more columns of a `"zresid"`
object. Optional metadata supplied through `info` (or stored in
attributes of `y`) are used only to fill legacy plotting attributes such
as `"type"`.

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
  my.mar = c(5, 4, 4, 6) + 0.1,
  legend.settings = list(),
  clip.extreme = TRUE,
  clip.threshold = 6,
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

  Legacy alias for `zcov`. When provided, it is used to fill missing
  legacy attributes such as `"type"`.

- irep:

  Integer vector specifying which column(s) of `y` to plot.

- diagnosis.test:

  Character string specifying the normality diagnostic to display.
  Currently `"SW"` is supported.

- main.title:

  Main title of the plot. If omitted, a default title is constructed
  from `attr(y, "type")`, when available.

- xlab:

  Label for the x-axis.

- ylab:

  Label for the y-axis.

- outlier.return:

  Logical; if `TRUE`, mark observations with `|Z| > outlier.value` and
  invisibly return their indices.

- outlier.value:

  Numeric threshold used to define outliers.

- outlier.set:

  A named list of graphical arguments passed to
  [`symbols`](https://rdrr.io/r/graphics/symbols.html) and
  [`text`](https://rdrr.io/r/graphics/text.html) when annotating
  outliers.

- my.mar:

  Numeric vector passed to
  [`par`](https://rdrr.io/r/graphics/par.html)`(mar = ...)`.

- legend.settings:

  Optional named list used to modify the default legend settings.

- clip.extreme:

  Logical; if `TRUE`, very large residuals are visually clipped to
  improve readability.

- clip.threshold:

  Numeric threshold used when `clip.extreme = TRUE`.

- ...:

  Additional graphical arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns a list with component `outliers`, containing the
indices of observations flagged as outliers for the plotted replicate.
The main effect of the function is the Q-Q plot.

## Details

The method can optionally report a Shapiro-Wilk normality diagnostic,
mark observations with large absolute residuals, and visually compress
extreme values when they would otherwise dominate the plot.

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
  qqnorm(z)
}
#> Warning: NaNs produced

```
