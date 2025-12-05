# Extract Posterior Predicted Parameters from a Hurdle or Count Model

Computes posterior predicted values for a specified distributional
parameter (e.g., mean, shape, or hurdle probability) from a fitted
Bayesian count or hurdle model. The function supports extracting
parameters for positive counts only or for all observations.

## Usage

``` r
posterior.pred(fit, dpar, count.only = TRUE)
```

## Arguments

- fit:

  A fitted brms model object containing data, formula, and MCMC
  posterior draws.

- dpar:

  Character string specifying the distributional parameter to extract:
  `"mu"` (mean), `"shape"` (dispersion), or `"zero"` (hurdle
  probability).

- count.only:

  Logical; if `TRUE` (default), computes predicted parameters only for
  positive counts (`y > 0`); otherwise, includes all observations.

## Value

A numeric matrix of predicted parameter values for each observation
(columns) and posterior draw (rows).

## Details

The function performs the following steps:

1.  Builds the model matrix for the chosen parameter and observation
    subset.

2.  Extracts the corresponding posterior MCMC draws from the fitted
    model.

3.  Computes the linear predictor via matrix multiplication of draws and
    model matrix.

4.  Applies the link function associated with the parameter (e.g.,
    logit, log) to obtain the predicted parameter values on their
    natural scale.

## See also

[`log.pred.dist.HNB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.HNB.md),
[`log.pred.dist.NB`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.NB.md),
[`log.pred.dist.HP`](https://tiw150.github.io/Zresidual/reference/log.pred.dist.HP.md),
[`posterior_linpred`](https://mc-stan.org/rstantools/reference/posterior_linpred.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Extract posterior predicted mean (mu) for all observations
mu_pred <- posterior.pred(fit, dpar = "mu", count.only = TRUE)

# Extract hurdle probabilities
pi_pred <- posterior.pred(fit, dpar = "zero", count.only = FALSE)
} # }
```
