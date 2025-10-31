#' Cumulative Distribution Function of the Hurdle Negative Binomial Distribution
#'
#' Computes the cumulative distribution function (CDF) or its logarithm for the
#' hurdle negative binomial (HNB) distribution. The hurdle model combines a point
#' mass at zero with a truncated negative binomial distribution for positive counts.
#'
#' @param y Numeric vector of observed count values.
#' @param mu Numeric vector of mean parameters of the negative binomial distribution.
#' @param size Numeric vector of shape (dispersion) parameters of the negative binomial distribution.
#' @param pi Numeric vector of hurdle probabilities (probability of structural zeros).
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le y)};
#'   otherwise, they are \eqn{P(Y > y)}.
#' @param log.p Logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @details
#' The hurdle negative binomial model assumes:
#' \deqn{
#' P(Y = 0) = \pi, \quad
#' P(Y = y \mid Y > 0) = (1 - \pi) \frac{F_{NB}(y)-F_{NB}(0)}{1 - F_{NB}(0)}, \quad y > 0
#' }
#' where \eqn{F_{NB}(y)} is CDF of the standard
#' negative binomial distribution.
#'
#' The function computes the upper or lower tail probabilities for both zeros and
#' positive counts using the logarithmic form for numerical stability. Internal helper
#' functions (\code{log_diff_exp}, \code{log_sum_exp}) are used to handle differences
#' and sums of log-scale probabilities safely.
#'
#' @return
#' A numeric vector of cumulative probabilities (or log-probabilities if \code{log.p = TRUE}).
#'
#' @examples
#' # Example: Hurdle Negative Binomial CDF
#' y <- 0:5
#' mu <- 2
#' size <- 1.5
#' pi <- 0.3
#' phurdle.nb(y, mu, size, pi)
#'
#' # Upper tail probabilities on log scale
#' phurdle.nb(y, mu, size, pi, lower.tail = FALSE, log.p = TRUE)
#'
#' @seealso
#' \code{\link[stats]{pnbinom}}, \code{\link[stats]{dnbinom}},
#' \code{\link{pdf.tnb}} for the zero-truncated negative binomial PMF.
#'
phurdle.nb <- function(y, mu, size, pi, lower.tail = FALSE, log.p = FALSE)
{
  log_diff_exp <- function (la, lb) # compute the logarithm of the difference between the exponentials of two log values.
  {
    ifelse(la >= lb, la + log1mexp(la-lb), NaN)
  }

  # compute the logarithm of the sum of exponential values.
  log_sum_exp <- function (lx)
  {
    mlx <- max (lx)
    log (sum (exp (lx - mlx))) + mlx
  }

  n.y <- max(length(y), length(mu),length(size),length(pi))
  y <- rep(y,length=n.y)
  mu <- rep(mu, length=n.y)
  size <- rep(size, length=n.y)
  pi <- rep(pi, length=n.y)

  # this code computes log upper tail probability
  log_upper_prob <- rep (0, n.y)
  i.py <- which (y>0) ## for y < 0, upper tail probability == 1 (log=0)
  i.zy <- which (y==0)

  ## computing upper tail probabilities for zero y only
  log_upper_prob [i.zy] <- log_diff_exp(0, log(pi[i.zy]))

  ## computing upper tail probabilities for positive y only
  log_upper_prob [i.py] <-
    log_diff_exp(0, log(pi[i.py])) +
    pnbinom(q=y[i.py],mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE) +
    - pnbinom(q=0,mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE)

  #- log_diff_exp(0, dnbinom(x=0,mu=mu[i.py],size=size[i.py],log=TRUE))
  upper_prob <- exp(log_upper_prob)
  log_lower_prob <- log_diff_exp(0, log_upper_prob)
  lower_prob <- exp (log_lower_prob)
  if (lower.tail==FALSE)
  {
    if(log.p) log_upper_prob
    else upper_prob
  } else
  {
    if(log.p) log_lower_prob
    else lower_prob
  }
}
