# Resolve or Identify Required point-wise predictive Functions

Determines which native point-wise predictive function should be
utilized for evaluating predictive tail probabilities. It dynamically
searches the calling environment, global environment, and package
namespace.

This function evaluates multiple naming conventions by extracting the
model's package namespace, S3 class, and statistical family. If no
suitable function is found (or if `show_names` is `TRUE`), it compiles a
descriptive message detailing the exact function names expected by the
package.

## Usage

``` r
required_log_pointpred(
  fit = NULL,
  pred_method = c("analytic", "simulation"),
  type = NULL,
  show_names = (sys.nframe() == 1),
  ...
)
```

## Arguments

- fit:

  A fitted model object. Used to automatically extract the model
  package, class, and family.

- pred_method:

  Character string indicating the prediction method. Must be one of
  `"analytic"` or `"simulation"`.

- type:

  Optional character string used as a component selector or variant tag.

- show_names:

  Logical; if `TRUE`, bypasses the environment search, prints the
  compiled string of expected function names to the console, and returns
  it invisibly. Defaults to `FALSE`.

- ...:

  Additional arguments (currently ignored).

## Value

The resolved function object. If no suitable function is found, it
returns a character string containing the compiled diagnostic message.
If `show_names = TRUE`, it prints the message and returns the string
invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Standard GLM
fit_glm <- glm(vs ~ mpg, data = mtcars, family = binomial())
required_log_pointpred(fit_glm, show_names = TRUE)

# Example 2: Survival Model (coxph)
library(survival)
fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
required_log_pointpred(fit_cox, show_names = TRUE)

# Example 3: Bayesian Model (brms) using pre-compiled package template
template_path <- system.file("extdata", "brms_template.rds", package = "Zresidual")
if (file.exists(template_path)) {
  fit_brm <- readRDS(template_path)
  required_log_pointpred(fit_brm, show_names = TRUE)
}
} # }
```
