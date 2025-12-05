# Plot martingale residuals for survival models

Produce diagnostic plots for martingale residuals from survival models,
with the option to plot against the observation index, the linear
predictor, or a selected covariate. The function expects a vector (or
one-column matrix) of martingale residuals with attributes attached by
the residual computation functions in this package.

## Usage

``` r
# S3 method for class 'martg.resid'
plot(
  x,
  ylab = "Martingale Residual",
  x_axis_var = c("index", "lp", "covariate"),
  main.title = "Martingale Residual Plot",
  outlier.return = FALSE,
  ...
)
```

## Arguments

- x:

  Numeric vector (or one-column matrix) of martingale residuals,
  typically returned by one of the residual functions in this package
  with `residual.type = "martingale"`. It must carry the attributes
  `"censored.status"`, `"linear.pred"`, and `"covariates"` as described
  above.

- ylab:

  Character string for the y-axis label. Default is
  `"Martingale Residual"`.

- x_axis_var:

  Character string controlling the x-axis. Must be one of `"index"`,
  `"lp"`, `"covariate"`, or the name of a covariate contained in
  `attr(x, "covariates")`. The default is effectively `"lp"` if
  `x_axis_var` is not supplied.

- main.title:

  Character string for the main plot title. Default is
  `"Martingale Residual Plot"`.

- outlier.return:

  Logical; if `TRUE`, attempted outliers (as indicated by an external
  logical vector `is.outlier` in the calling environment) are
  highlighted in the plot and their indices are returned invisibly. If
  `FALSE` (default), no outlier indices are returned. Note that this
  function does not compute outliers internally: it assumes that a
  logical vector `is.outlier` of the same length as `x` is available if
  outlier highlighting is desired.

- ...:

  Additional arguments passed to the underlying plotting functions.

## Value

The function is primarily called for its side-effect of producing a
plot. If `outlier.return = TRUE`, it prints the indices of outlying
points to the console and invisibly returns a list with component
`outliers`, containing the indices where `is.outlier` is `TRUE`.
Otherwise, it returns `NULL` invisibly.

## Details

The input `x` is typically obtained from
[`residual.coxph()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.md),
[`residual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.frailty.md),
or
[`residual.survreg()`](https://tiw150.github.io/Zresidual/reference/residual.survreg.md)
with `residual.type = "martingale"`.

The `x_axis_var` argument controls the x-axis:

- `"index"`: plot martingale residuals against observation index.

- `"lp"`: plot martingale residuals against the linear predictor
  (attribute `"linear.pred"`).

- `"covariate"`: prompt the user and print the available covariate names
  to the console.

- a character string matching one of the covariate names in
  `attr(x, "covariates")`: plot martingale residuals against that
  covariate.

In the `"lp"` and covariate cases, a LOWESS smooth is added to the plot
to highlight systematic patterns in the residuals.

Non-finite martingale residuals are detected and truncated to lie
slightly beyond the largest finite residual, with a warning message
printed to alert the user that there may be problems with the model fit.
Censored and uncensored observations are distinguished by color and
plotting symbol in all display modes.

## See also

`residual.coxph`, `residual.coxph.frailty`, `residual.survreg`,
[`Surv`](https://rdrr.io/pkg/survival/man/Surv.html),
[`survfit`](https://rdrr.io/pkg/survival/man/survfit.html)

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)

data(lung)
fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
r_m <- residual.coxph(fit, newdata = lung,
                      residual.type = "martingale")

## Basic plot vs. index
plot(r_m, x_axis_var = "index")

## Plot vs. linear predictor
plot(r_m, x_axis_var = "lp")

## Plot vs. a specific covariate, e.g. "age"
plot(r_m, x_axis_var = "age")
} # }
```
