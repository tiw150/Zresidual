# Plot Z-Residuals for Bayesian and Frequentist Count / Hurdle / Zero / Survival Models

Produces diagnostic scatterplots of Z-residuals for a wide range of
model types, including hurdle models, zero models, count models, and
survival models. This function is designed to be compatible with
Z-residual matrices generated from Bayesian models (e.g., **brms**,
Stan-based models) as well as classical models.

The function supports:

- Plotting residuals against index, covariates, linear predictors, or
  user-specified vectors.

- Visual outlier detection with customizable coloring, emphasis, and
  labels.

- Automatic handling of censored/un-censored and hurdle count/zero
  classifications.

- Multiple normality tests (Shapiro-Wilk, ANOVA, Bartlett) with
  per-iteration reporting.

- Extensive customization through graphical parameters.

## Usage

``` r
# S3 method for class 'zresid'
plot(
  x,
  irep = 1:ncol(x),
  ylab = "Z-Residual",
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  x_axis_var = c("index", "covariate", "lp"),
  main.title = ifelse(is.null(attr(x, "type")),
                      "Z-residual Scatterplot",
                      paste("Z-residual Scatterplot -", attr(x, "type"))),
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

- x:

  A numeric matrix of Z-residuals (dimensions: *n* × *m*) of class
  `"zresid"`. Attributes used by the function include: `"type"`,
  `"zero_id"`, `"censored.status"`, `"covariates"`, and `"linear.pred"`.

- irep:

  A vector specifying the residual column(s) (iterations) to plot.
  Defaults to all columns of **`x`**.

- ylab:

  Label for the y-axis.

- normality.test:

  A character vector specifying normality tests applied per iteration:
  `"SW"`, `"AOV"`, or `"BL"`. Helper functions (e.g., `sw.test.zresid`)
  must exist and accept inputs (`x`, `x_axis_var`, `k.test`).

- k.test:

  Bin size for normality tests (used for grouping continuous
  predictors).

- x_axis_var:

  Specifies the x-axis values. Options include: `"index"`,
  `"covariate"`, `"lp"`, a covariate name present in
  `attr(x, "covariates")`, or a numeric vector of length *n*.

- main.title:

  Title of the plot. Defaults to a type-based informative title.

- outlier.return:

  If `TRUE`, outliers are printed to console and returned invisibly.

- outlier.value:

  Threshold above which a residual is flagged as an outlier. Defaults to
  `3.5`.

- category:

  Optional vector categorizing observations (length *n*). Used for
  coloring and shaping points in scatterplots.

- outlier.set:

  A named list of arguments passed to
  [`symbols()`](https://rdrr.io/r/graphics/symbols.html) and
  [`text()`](https://rdrr.io/r/graphics/text.html) for marking and
  labeling outliers. Overrides defaults.

- xlab:

  Label for the x-axis. May include LaTeX syntax using the form
  `tex("...")`, which will be interpreted via **latex2exp** (if
  installed).

- my.mar:

  A numeric vector passed to `par(mar=...)` to adjust plot margins.

- ...:

  Additional graphical arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`legend()`](https://rdrr.io/r/graphics/legend.html),
  [`symbols()`](https://rdrr.io/r/graphics/symbols.html), and
  [`text()`](https://rdrr.io/r/graphics/text.html).

## Value

Invisibly returns (when `outlier.return = TRUE`) a list:

- outliers:

  Vector of outlier indices

Otherwise returns `NULL`. Always produces a scatterplot as its primary
output.

## Details

This function employs S3 method dispatch. The input `x` must be an
object of class `"zresid"`.

### **Model-Type-Specific Behavior**

#### **Hurdle models**

- Zeros and counts are colored differently (`red` vs `blue`).

#### **Survival models**

- Censored and uncensored observations are visually separated.

### **Outlier Detection**

Outliers are defined as:

\$\$\|Z\| \> \mbox{outlier.value} \hspace{1em} \mbox{or non-finite
values}\$\$

They are marked, labeled, and returned to the user if
`outlier.return = TRUE`.

## Note

This function modifies graphical parameters (`par(mar=...)`) during
execution and resets them at the end.
