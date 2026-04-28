# Compute Z-residuals from predictive tail probabilities

Computes Z-residuals by transforming predictive upper-tail probabilities
to the standard normal scale. When point-wise predictive mass at the
observed value is available, randomized Z-residuals are constructed by
adding a uniform perturbation within that mass. This allows the same
interface to handle both continuous-type and discrete-type outcomes.

## Usage

``` r
Zresidual(
  fit = NULL,
  data = NULL,
  log_pointpred = NULL,
  pred_method = c("analytic", "simulation"),
  mcmc_summarize = c("post", "iscv"),
  type = NULL,
  randomized = TRUE,
  nrep = 30,
  eps = 1e-12,
  seed = NULL,
  log_sumpred = NULL,
  ...
)
```

## Arguments

- fit:

  A fitted model object. Optional if `log_sumpred` is provided.

- data:

  A data frame or list containing the data used to evaluate predictive
  quantities. Optional if `log_sumpred` is provided.

- log_pointpred:

  Optional function (or function name as a character string) returning
  predictive pointwise quantities. If `NULL`, the function name is
  automatically constructed based on the model backend, family, and
  `pred_method`.

- pred_method:

  Character string indicating the prediction method. Must be one of
  `"analytic"` or `"simulation"`.

- mcmc_summarize:

  Character string indicating how posterior predictive quantities are
  summarized. Must be one of `"post"` (posterior mean) or `"iscv"`
  (importance sampling cross-validation).

- type:

  Optional character string used as a component selector or variant tag.

- randomized:

  Logical; if `TRUE`, generate randomized Z-residuals when point mass or
  discrete-type contribution is available.

- nrep:

  Integer giving the number of residual replicates to generate. If
  `randomized = FALSE`, this is forced to `1L`.

- eps:

  Small positive numeric constant used to keep probability values
  strictly bounded away from exactly 0 and 1.

- seed:

  Optional integer seed for reproducible residual randomization.

- log_sumpred:

  Optional precomputed output list from
  [`log_summary_pred`](https://tiw150.github.io/Zresidual/reference/log_summary_pred.md)
  to bypass MCMC summarization.

- ...:

  Additional arguments passed to the underlying prediction backend.

## Value

A numeric matrix with one row per observation and `nrep` columns
containing the computed Z-residual replicates. The returned object has
class `"zresid"`. Attribute `"rsp"` stores the probability-scale
residuals used to form the Z-residuals. Other contextual attributes
(e.g., covariates, linear predictors) may be attached if `Zcov`
successfully extracts them.

## Details

To avoid repeating expensive MCMC summaries across multiple
randomizations, users can precompute the predictive summaries using
[`log_summary_pred`](https://tiw150.github.io/Zresidual/reference/log_summary_pred.md)
and pass the result to the `log_sumpred` argument.

## See also

[`log_summary_pred`](https://tiw150.github.io/Zresidual/reference/log_summary_pred.md)
