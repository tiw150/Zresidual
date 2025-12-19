# Z-residuals for Standard Cox Models (Internal Worker)

Internal function to compute randomized Z-residuals for a **standard**
Cox proportional hazards model (without frailty).

## Usage

``` r
Zresidual_coxph_survival(fit_coxph, newdata, n.rep = 1, ...)
```

## Arguments

- fit_coxph:

  A fitted
  [`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html) model
  object.

- newdata:

  Optional data frame on which to compute the residuals. If `NULL`, the
  original model frame is used.

- n.rep:

  Integer. Number of randomized residual replicates to generate. Default
  is 1.

- ...:

  Additional arguments (currently unused).

## Value

A matrix containing the Z-residuals (\\N \times nrep\\) with diagnostic
attributes: `Survival.Prob`, `linear.pred`, `covariates`,
`censored.status`, `object.model.frame`, and `type = "survival"`.

## Details

The Z-residual for an observation \$i\$ is calculated as \$\$Z_i =
-\Phi^{-1}(\hat{S}\_i(t_i, \text{rand}))\$\$ where \\\Phi^{-1}\\ is the
inverse standard normal CDF, and \\\hat{S}\_i(t_i, \text{rand})\\ is the
predicted survival probability at time \\t_i\\.

For uncensored observations (\\t_i\\ is an event time),
\\\hat{S}\_i(t_i, \text{rand}) = \hat{S}\_i(t_i)\\. For censored
observations, \\\hat{S}\_i(t_i, \text{rand}) = \hat{S}\_i(t_i) \cdot
U\\, where \\U \sim \text{Unif}(0, 1)\\.

The predicted survival is calculated as \$\$\hat{S}\_i(t_i) =
\exp(-\exp(\mathbf{x}\_i \mathbf{\hat{\beta}}) \hat{H}\_0(t_i))\$\$.
