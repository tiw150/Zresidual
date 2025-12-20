# An Animation to Explain Z-residuals

## Z-residuals

The **Zresidual** package implements diagnostic residuals based on the
**predictive distribution** of each observation. By utilizing the full
probabilistic information of the model, the package generates residuals
that are approximately normally distributed, enabling further standard
diagnostics as using Pearson’s residuals for OLS.

For a given observation y_i, the randomized survival probability (also
known as randomized probability integral transform) represents the value
of the predictive survival function (CDF) at y_i, with a randomization
term to ensure continuity for discrete outcomes. It is defined as:

RSP_i(y_i \mid \theta) = S_i(y_i \mid \theta) + U_i \\ p_i(y_i \mid
\theta)

where S_i and p_i represent the survival function and probability mass
function, respectively, derived from the predictive distribution of y_i
given covariates \mathbf{x}\_i and parameters \theta. U_i \sim
\text{Unif}(0,1) is a random uniform variable used to handle discrete
outcomes. Under a correctly specified model where the observed data y_i
arises from the assumed predictive distribution, the RSP follows a
\text{Unif}(0,1) distribution.

The **Z-residual** (also known as randomized quantile residual) is then
derived by transforming the RSP via the inverse standard normal
cumulative distribution function:

z_i^{RSP}(y_i \mid \theta) = -\Phi^{-1}\[RSP_i(y_i \mid \theta)\]

Under correct model specification, these residuals follow a standard
normal N(0,1) distribution. This mapping allows researchers to assess
the quality of the predictive distribution—including its mean, variance,
and shape—using standard diagnostics for Gaussian OLS.
