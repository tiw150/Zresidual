# An Animation to Understand Z-residuals

## 1 Z-residuals

The **Zresidual** package implements diagnostic residuals based on the
**predictive distribution** of each observation. By utilizing the full
probabilistic information of the model, the package generates residuals
that are approximately normally distributed. This allows for standard
diagnostic techniques similar to using Pearson residuals in OLS
regression.

For a given observation y_i, the **Randomized Survival Probability
(RSP)**—also known as the randomized probability integral
transform—represents the value of the predictive survival function
(upper-tail probability) at y_i, adjusted with a randomization term to
ensure continuity for discrete outcomes. It is defined as:

RSP_i(y_i \mid \theta) = S_i(y_i \mid \theta) + U_i \\ p_i(y_i \mid
\theta) \tag{1}

where S_i and p_i represent the survival function and probability mass
function, respectively, derived from the predictive distribution of y_i
given covariates \mathbf{x}\_i and parameters \theta. U_i \sim
\text{Unif}(0,1) is a random uniform variable used to smooth discrete
outcomes. Under a correctly specified model—where the observed data y_i
arises from the assumed predictive distribution—the RSP follows a
\text{Unif}(0,1) distribution.

The **Z-residual** (or randomized quantile residual) is derived by
transforming the RSP via the inverse standard normal cumulative
distribution function (CDF):

z_i(y_i \mid \theta) = -\Phi^{-1}\left(RSP_i(y_i \mid \theta)\right)=
\Psi^{-1}\left(RSP_i(y_i \mid \theta)\right) \tag{2}

where \Psi^{-1}(p) = -\Phi^{-1}(p) = \Phi^{-1}(1-p) represents the
inverse survival function (upper-tailed quantile) of the standard normal
distribution N(0,1).

In simpler terms, the Z-residual of y_i is the N(0,1) quantile that
shares the same **upper-tail probability** as y_i in its predictive
distribution. For discrete y_i, this position is discrete (a step
function), so randomization is applied to fill the gap. This
transformation effectively maps any predictive distribution to N(0,1).
If the model is correctly specified, the resulting Z-residuals will
follow a standard normal distribution.

The animation in **[Figure 1](#fig-aniz)** demonstrates this
transformation.

- **Top Plot:** A value y_i is randomly drawn from a Negative Binomial
  distribution (\mu=5, dispersion parameter size =5). The **dark blue**
  area represents the survival probability S(y_i), while the **light
  blue** area represents the random jitter U_i p(y_i). The combined
  shaded area equals the RSP_i.
- **Bottom Plot:** This total area (probability mass) is mapped to the
  upper tail of the standard normal distribution. The resulting z_i is
  marked by a red dot.

As more points are simulated, the distribution of z_i values (shown by
the accumulating points) aligns with the theoretical N(0,1) density
curve.

![](z_residual_anim.gif)

Figure 1: Animation to Demonstrate the Z-residuals given the true
parameter.
