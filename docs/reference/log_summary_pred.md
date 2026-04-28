# Extract and summarize MCMC predictive tail probabilities

Evaluates observation-level predictive log-survival probabilities and
log-probability masses (for discrete data) over the MCMC draws, and
aggregates them into a single summary vector per observation based on
the chosen `mcmc_summarize` method. This function isolates the
computationally expensive matrix operations (such as log-sum-exp) from
the randomization step.

## Usage

``` r
log_summary_pred(
  fit,
  data,
  log_pointpred = NULL,
  pred_method = c("analytic", "simulation"),
  mcmc_summarize = c("post", "iscv"),
  type = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object.

- data:

  A data frame or list containing the data used to evaluate predictive
  quantities. Must be provided.

- log_pointpred:

  Optional function (or function name as a character string) returning
  predictive pointwise quantities. If `NULL`, the function is
  automatically resolved based on the model backend, family, and
  `pred_method` using
  [`required_log_pointpred`](https://tiw150.github.io/Zresidual/reference/required_log_pointpred.md).
  The resolved function must return a list containing `log_surv`,
  `log_like`, and `is_discrete`.

- pred_method:

  Character string indicating the prediction method. Must be one of
  `"analytic"` or `"simulation"`.

- mcmc_summarize:

  Character string indicating how posterior predictive quantities are
  summarized. Must be one of `"post"` (posterior mean) or `"iscv"`
  (importance sampling cross-validation).

- type:

  Optional character string used as a component selector or variant tag.

- ...:

  Additional arguments passed to the underlying prediction backend.

## Value

A list containing three elements:

- `log_surv_hat`:

  A numeric vector of the summarized log-survival probabilities.

- `log_pmf_hat`:

  A numeric vector of the summarized log-probability masses (or `NA` for
  continuous observations).

- `is_discrete`:

  An integer vector indicating whether each observation is treated as
  discrete (`1`) or continuous (`0`).

## See also

[`Zresidual`](https://tiw150.github.io/Zresidual/reference/Zresidual.md),
[`required_log_pointpred`](https://tiw150.github.io/Zresidual/reference/required_log_pointpred.md)
