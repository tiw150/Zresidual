# Cox–Snell residual plot for survival models

Produce a Cox–Snell residual diagnostic plot for survival models, based
on the cumulative hazard of the Cox–Snell residuals. Under a correctly
specified model, the Cox–Snell residuals should follow an exponential
distribution with mean 1, so the nonparametric estimate of the
cumulative hazard should lie close to the 45-degree line.

## Usage

``` r
# S3 method for class 'cs.residual'
plot(
  cs.residual,
  ylab = "Cumulative Hazard Function",
  main.title = "Cox-Snell Residuals Scatterplot",
  outlier.return = FALSE,
  ...
)
```

## Arguments

- cs.residual:

  Numeric vector (or one-column matrix) of Cox–Snell residuals,
  typically returned by one of the residual functions in this package
  with `residual.type = "Cox-Snell"`. It must carry the attribute
  `"censored.status"`, giving the event indicator (1 = event, 0 =
  censored).

- ylab:

  Character string for the y-axis label. Default is
  `"Cumulative Hazard Function"`.

- main.title:

  Character string for the main plot title. Default is
  `"Cox-Snell Residuals Scatterplot"`.

- outlier.return:

  Logical; if `TRUE`, potential outliers are identified using a simple
  cutoff `|cs.residual| > 3.5`. Their indices are printed to the console
  and returned invisibly. If `FALSE` (default), no outlier indices are
  returned. Note that the current implementation attempts to highlight
  outliers using additional plotting calls and assumes access to objects
  named `Zresidual` and `j` in the calling environment; users may wish
  to adapt this part of the code for their own workflows.

- ...:

  Additional arguments passed to
  [`plot.survfit`](https://rdrr.io/pkg/survival/man/plot.survfit.html)
  or to the underlying base graphics functions.

## Value

The function is primarily called for its side-effect of producing a
plot. If `outlier.return = TRUE`, it prints the indices of points
flagged as outliers (`|cs.residual| > 3.5`) and invisibly returns a list
with component `outliers`, containing these indices. Otherwise, it
returns `NULL` invisibly.

## Details

The input `cs.residual` is typically obtained from the residual
functions in this package (e.g.,
[`residual.coxph()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.md),
[`residual.coxph.frailty()`](https://tiw150.github.io/Zresidual/reference/residual.coxph.frailty.md),
or
[`residual.survreg()`](https://tiw150.github.io/Zresidual/reference/residual.survreg.md))
with `residual.type = "Cox-Snell"`.

Non-finite Cox–Snell residuals are detected and truncated to lie
slightly beyond the largest finite residual, with a message printed to
alert the user that there may be problems with the model fit. The
cumulative hazard is drawn using `fun = "cumhaz"` and compared visually
to the exponential(1) reference line \\H(t) = t\\.

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
  cs_resid <- residual.coxph(fit, newdata = lung,
                             residual.type = "Cox-Snell")

  ## Cox–Snell residual plot
  plot.cs.residual(cs_resid)

  ## Return indices of large residuals
  out <- plot.cs.residual(cs_resid, outlier.return = TRUE)
  out$outliers
} # }
```
