# Boxplot of Z-Residuals

Produces a boxplot of Z-residuals grouped by binned fitted values or a
selected covariate. This diagnostic plot supports other count-data
models (e.g. Bayesian hurdle, zero-truncated) by visualizing residual
distribution, detecting outliers, and evaluating normality assumptions
using Shapiro-Wilk, ANOVA, or Bartlett-type tests for Z-residuals.

## Usage

``` r
# S3 method for class 'zresid'
boxplot(
  x,
  irep = 1,
  x_axis_var = c("lp", "covariate"),
  num.bin = 10,
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  main.title = paste("Z-residual Boxplot -", attr(x, "type")),
  outlier.return = FALSE,
  outlier.value = 3.5,
  ...
)
```

## Arguments

- x:

  A matrix of Z-residuals with one column per MCMC iteration. Must
  contain attributes:

  `"type"`

  :   Model type used to generate residuals (e.g., hurdle, truncated).

  `"fitted.value"`

  :   Vector of fitted values for the model.

  `"covariates"`

  :   Optional data frame of covariates.

  `"zero_id"`

  :   Indices of zero observations (if applicable).

- irep:

  Integer vector indicating which columns of `Zresidual` to plot.
  Default is `1`.

- x_axis_var:

  Character string specifying the x-axis variable:

  - `"fitted.value"` (default): Bin fitted values.

  - `"covariate"`: Display a list of covariate names.

  - A specific covariate name present in
    `attr(Zresidual, "covariates")`.

- num.bin:

  Integer. Number of bins for grouping fitted values or selected
  covariate. Defaults to `10`.

- normality.test:

  Character vector specifying which normality tests to report:

  - `"SW"` - Shapiro-Wilk test for Z-residuals.

  - `"AOV"` - ANOVA-based test for variance/mean structure.

  - `"BL"` - Bartlett-type test for variance homogeneity.

  Defaults to `c("SW","AOV","BL")`.

- k.test:

  Integer. Number of groups to use for ANOVA/Bartlett-type tests.
  Default is `10`.

- main.title:

  Character. Main title of the plot. Default includes the model type
  automatically.

- outlier.return:

  Logical. If `TRUE`, returns the index of Z-residual values exceeding
  the threshold defined by `outlier.value`. Default is `FALSE`.

- outlier.value:

  Numeric. Threshold for defining outliers based on absolute Z-residual
  magnitude. Default is `3.5`.

- ...:

  Additional graphical parameters passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) or
  [`legend`](https://rdrr.io/r/graphics/legend.html).

## Value

If `outlier.return = TRUE`, returns a list containing:

- `outliers` - vector of indices where `|Zresidual| > outlier.value`.

Otherwise, the function returns `NULL` (invisible) and produces a
diagnostic plot.

## Details

The function generates boxplots of Z-residuals across binned fitted
values (or a selected covariate), which helps detect lack of fit,
heteroscedasticity, and model misspecification.

Infinite and non-finite residuals are automatically replaced with a
maximal finite value (with preserved sign), and a warning message is
displayed.

When `x_axis_var="covariate"`, users may supply any covariate name
available in the `"covariates"` attribute. If a covariate contains too
few unique bins, fitted values are transformed using
[`log()`](https://rdrr.io/r/base/Log.html) to stabilize binning, with a
message provided.

Normality diagnostics are displayed in the plot legend. Internally, the
function calls:
[`sw.test.zresid()`](https://tiw150.github.io/Zresidual/reference/sw.test.zresid.md),
[`aov.test.zresid()`](https://tiw150.github.io/Zresidual/reference/aov.test.zresid.md),
and
[`bartlett.test.zresid()`](https://tiw150.github.io/Zresidual/reference/bartlett.test.zresid.md).

## See also

[`plot`](https://rdrr.io/r/graphics/plot.default.html),
[`boxplot`](https://rdrr.io/r/graphics/boxplot.html),
[`legend`](https://rdrr.io/r/graphics/legend.html), `sw.test.zresid`,
`aov.test.zresid`, `bartlett.test.zresid`

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'zres' is a Z-residual matrix produced by a model-fitting function:
boxplot.zresid(zres)

# Plot against a specific covariate
boxplot.zresid(zres, x_axis_var = "age")

# Return outliers
box.out <- boxplot.zresid(zres, outlier.return = TRUE)
} # }
```
