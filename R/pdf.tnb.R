#' Probability Mass Function of the Zero-Truncated Negative Binomial Distribution
#'
#' Computes the probability mass function (PMF) or log-PMF for the zero-truncated
#' negative binomial (TNB) distribution. This version excludes zeros and rescales
#' the probabilities accordingly so that they sum to one over positive counts only.
#'
#' @param y Numeric vector of observed count values (\code{y > 0}).
#' @param mu Numeric vector of mean parameters for the negative binomial distribution.
#' @param size Numeric vector of shape (dispersion) parameters.
#' @param log.p Logical; if \code{TRUE}, returns log probabilities instead of probabilities.
#'
#' @details
#' The zero-truncated negative binomial probability for an observation \eqn{y > 0} is:
#' \deqn{
#' P(Y = y \mid Y > 0) = \frac{P(Y = y)}{1 - P(Y = 0)}
#' }
#' where \eqn{P(Y = y)} and \eqn{P(Y = 0)} are evaluated using the standard
#' negative binomial PMF and CDF, respectively. The implementation uses
#' \code{\link[stats]{dnbinom}} and \code{\link[stats]{pnbinom}} for computation.
#'
#' The function automatically vectorizes inputs, ensuring that the output
#' corresponds elementwise to each set of parameters.
#'
#' @return A numeric vector of probabilities (or log-probabilities if \code{log.p = TRUE}).
#'
#' @examples
#' # Example: Zero-truncated negative binomial probabilities
#' y <- 1:5
#' mu <- 2
#' size <- 1.5
#' pdf.tnb(y, mu, size)
#'
#' # Log probabilities
#' pdf.tnb(y, mu, size, log.p = TRUE)
#'
#' @seealso
#' \code{\link[stats]{dnbinom}}, \code{\link[stats]{pnbinom}},
#' and \code{\link{cdf.tnb}} for the corresponding cumulative function.
#'
pdf.tnb <- function(y, mu, size, log.p = FALSE)
{

  n.y <- max(length(y), length(mu),length(size),length(pi))
  y <- rep(y,length=n.y)
  mu <- rep(mu, length=n.y)
  size <- rep(size, length=n.y)

  log_prob <- rep (-Inf, n.y)
  #if(count_only) i.py <- which(y>0) else i.py <- 1:n
  i.py <- which (y>0) ## for y > 0

  ## computing log probabilities for positive y only (y>0)
  log_prob [i.py] <- dnbinom(x=y[i.py],mu=mu[i.py],size=size[i.py], log = TRUE) -
    pnbinom(q=0,mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE)

  prob <- exp(log_prob)

  if(log.p) log_prob
  else prob

}
