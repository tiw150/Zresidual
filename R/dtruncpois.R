#' Zero-truncated Poisson distribution
#'
#' Density and distribution functions for the zero-truncated Poisson
#' distribution.
#'
#' @param y Numeric vector of counts.
#' @param lambda Mean parameter of the Poisson distribution.
#' @param log.p Logical; if \code{TRUE}, return values on the log scale. For
#'   \code{dtruncpois()} this returns log-densities; for \code{ptruncpois()} this
#'   returns log-probabilities.
#' @param lower.tail Logical; used only by \code{ptruncpois()}. If \code{TRUE},
#'   return \eqn{P(Y \le y)}; otherwise return \eqn{P(Y > y)}.
#'
#' @return A numeric vector of densities, distribution values, or their
#' logarithms.
#'
#' @examples
#' dtruncpois(1:4, lambda = 2)
#' ptruncpois(1:4, lambda = 2)
#'
#' @name truncpois
#' @export
dtruncpois <- function(y, lambda, log.p = FALSE)
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
