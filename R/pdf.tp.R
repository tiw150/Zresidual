#' Probability Mass Function of the Zero-Truncated Poisson Distribution
#'
#' Computes the probability mass function (PMF) or log-PMF for the zero-truncated
#' Poisson (TP) distribution. This version excludes zeros and rescales
#' the probabilities so that they sum to one over positive counts only.
#'
#' @param y Numeric vector of observed count values (\code{y > 0}).
#' @param lambda Numeric vector of rate parameters (mean of the Poisson distribution).
#' @param log.p Logical; if \code{TRUE}, returns log probabilities instead of probabilities.
#'
#' @details
#' The zero-truncated Poisson probability for an observation \eqn{y > 0} is:
#' \deqn{
#' P(Y = y \mid Y > 0) = \frac{P(Y = y)}{1 - P(Y = 0)}
#' }
#' where \eqn{P(Y = y)} and \eqn{P(Y = 0)} are evaluated using the standard
#' Poisson PMF and CDF, respectively. The function uses
#' \code{\link[stats]{dpois}} and \code{\link[stats]{ppois}} internally.
#'
#' This function automatically vectorizes inputs so that each probability
#' corresponds elementwise to the provided parameter values.
#'
#' @return A numeric vector of probabilities (or log-probabilities if \code{log.p = TRUE}).
#'
#' @examples
#' # Example: Zero-truncated Poisson probabilities
#' y <- 1:5
#' lambda <- 2
#' pdf.tp(y, lambda)
#'
#' # Log probabilities
#' pdf.tp(y, lambda, log.p = TRUE)
#'
#' @seealso
#' \code{\link[stats]{dpois}}, \code{\link[stats]{ppois}},
#' and \code{\link{cdf.tp}} for the corresponding cumulative function.
#'
pdf.tp <- function(y, lambda, log.p = FALSE)
{

  n.y <- max(length(y), length(lambda))
  y <- rep(y,length=n.y)
  lambda <- rep(lambda, length=n.y)

  log_prob <- rep (-Inf, n.y)
  i.py <- which (y>0) ## for y > 0

  ## computing log probabilities for positive y only (y>0)
  log_prob [i.py] <- dpois(x=y[i.py],lambda=lambda[i.py], log = TRUE) -
    ppois(q=0,lambda=lambda[i.py], lower.tail = FALSE, log.p = TRUE)

  prob <- exp(log_prob)

  if(log.p) log_prob
  else prob

}
