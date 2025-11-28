# Plot Z-Residuals for Bayesian and Frequentist Count / Hurdle / Zero / Survival Models

Produces diagnostic scatterplots of Z-residuals for a wide range of
model types, including hurdle models, zero models, count models, and
survival models. This function is designed to be compatible with
Z-residual matrices generated from Bayesian models (e.g., \*\*brms\*\*,
Stan-based models) as well as classical models.

The function supports: - Plotting residuals against index, covariates,
linear predictors, or user-specified vectors. - Visual outlier detection
with customizable coloring, emphasis, and labels. - Automatic handling
of censored/un-censored and hurdle count/zero classifications. -
Multiple normality tests (Shapiro-Wilk, ANOVA, Bartlett) with
per-iteration reporting. - Extensive customization through graphical
parameters.

It is intended to support model diagnostics, specifically for detecting
violations in distributional assumptions in hierarchical and Bayesian
modeling.

## Usage

``` r
plot.zresid(
  Zresidual,
  irep = 1:ncol(Zresidual),
  ylab = "Z-Residual",
  normality.test = c("SW", "AOV", "BL"),
  k.test = 10,
  X = c("index", "covariate", "lp"),
  main.title = ifelse(is.null(attr(Zresidual, "type")),
                      "Z-residual Scatterplot",
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

  A numeric matrix of Z-residuals (dimensions: \*n\* × \*m\*), typically
  generated from Bayesian models. Attributes used by the function
  include: - \`"type"\`: character string (\`"hurdle"\`, \`"count"\`,
  \`"zero"\`, \`"survival"\`, or \`NULL\`) - \`"zero_id"\`: indices of
  zero model observations (for hurdle models) - \`"censored.status"\`:
  binary vector indicating censoring status for survival models -
  \`"covariates"\`: model matrix of covariates - \`"linear.pred"\`:
  linear predictor values

- irep:

  A vector specifying the residual column(s) (iterations) to plot.
  Defaults to all columns.

- ylab:

  Label for the y-axis.

- normality.test:

  A character vector specifying normality tests applied per iteration: -
  \`"SW"\`: Shapiro-Wilk test using \`sw.test.zresid()\` - \`"AOV"\`:
  ANOVA-type normality check using \`aov.test.zresid()\` - \`"BL"\`:
  Bartlett test using \`bartlett.test.zresid()\`

  These functions must exist in the environment and accept inputs
  \`(Zresidual, X, k.test)\`.

- k.test:

  Bin size for normality tests (used for grouping continuous
  predictors).

- X:

  Specifies the x-axis values. Options include: - \`"index"\`:
  observation index - \`"covariate"\`: first covariate from
  \`attr(Zresidual, "covariates")\` - \`"lp"\`: linear predictor - A
  covariate name present in the covariate matrix attribute - A numeric
  vector of length \*n\* (user-specified)

- main.title:

  Title of the plot. Defaults to a type-based informative title.

- outlier.return:

  If \`TRUE\`, outliers are printed to console and returned invisibly.

- outlier.value:

  Threshold above which a residual is flagged as an outlier. Defaults to
  \`3.5\`.

- category:

  Optional vector categorizing observations (length \*n\*). Used for
  coloring and shaping points in scatterplots.

- outlier.set:

  A named list of arguments passed to \`symbols()\` and \`text()\` for
  marking and labeling outliers. Overrides defaults.

- xlab:

  Label for the x-axis. May include LaTeX syntax using the form
  \`tex("...")\`, which will be interpreted via \*\*latex2exp\*\* (if
  installed).

- my.mar:

  A numeric vector passed to \`par(mar=...)\` to adjust plot margins.

- ...:

  Additional graphical arguments passed to \`plot()\`, \`legend()\`,
  \`symbols()\`, and \`text()\`.

## Value

Invisibly returns (when \`outlier.return = TRUE\`) a list:

- outliers:

  Vector of outlier indices

Otherwise returns \`NULL\`.

Always produces a scatterplot as its primary output.

## Details

\## \*\*Model-Type-Specific Behavior\*\*

\### \*\*Hurdle models\*\* - Zeros and counts are colored differently
(\`red\` vs \`blue\`). - A legend is automatically created.

\### \*\*Survival models\*\* - Censored and uncensored observations are
visually separated.

\### \*\*Category-Based Styling\*\* If \`category\` is provided: -
Unique categories receive distinct colors and point shapes. - Users can
override palettes via \`col=\` or \`pch=\` in \`...\`.

\## \*\*Normality Testing and Legends\*\* - Normality tests are
performed \*per iteration\*. - P-values are displayed as a right-side
legend.

\## \*\*Outlier Detection\*\* Outliers are defined as:

\$\$\|Z\| \> \text{outlier.value} \quad \text{or non-finite values}\$\$

They are: - Marked using \`symbols()\` - Labeled using \`text()\` -
Returned to the user if \`outlier.return = TRUE\`

\## \*\*Handling Infinite or NaN Residuals\*\* - Infinite values are
replaced by signed \`(max + 0.1)\` for plotting. - A warning message is
issued.

## Note

This function modifies graphical parameters (\`par(mar=...)\`) during
execution and resets them at the end.

## References

Dunn, P.K., & Smyth, G.K. (1996). Randomized Quantile Residuals.
\*Journal of Computational and Graphical Statistics\*.

Gelman, A., Carlin, J.B., Stern, H.S., Dunson, D., Vehtari, A., & Rubin,
D.B. (2013). \*Bayesian Data Analysis\*.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose Z is a matrix of Z-residuals from a Bayesian hurdle model

attr(Z, "type")      <- "hurdle"
attr(Z, "zero_id")   <- which(y == 0)
attr(Z, "covariates") <- model.matrix(~ x1 + x2)
attr(Z, "linear.pred") <- fitted_values

# Basic plot against index
plot.zresid(Z)

# Plot against a covariate
plot.zresid(Z, X = "x1")

# User-defined x vector
plot.zresid(Z, X = fitted_values, normality.test = "SW")

# Custom outlier appearance
plot.zresid(Z, outlier.set = list(col = "red", cex = 1.4))

# With categories
plot.zresid(Z, category = group)
} # }
```
