#' Zero-truncated negative binomial distribution
#'
#' Density and distribution functions for the zero-truncated negative binomial
#' distribution.
#'
#' @param y Numeric vector of counts.
#' @param mu Mean parameter of the negative binomial distribution.
#' @param size Dispersion parameter of the negative binomial distribution.
#' @param log.p Logical; if \code{TRUE}, return values on the log scale. For
#'   \code{dtruncnb()} this returns log-densities; for \code{ptruncnb()} this
#'   returns log-probabilities.
#' @param lower.tail Logical; used only by \code{ptruncnb()}. If \code{TRUE},
#'   return \eqn{P(Y \le y)}; otherwise return \eqn{P(Y > y)}.
#'
#' @return A numeric vector of densities, distribution values, or their
#' logarithms.
#'
#' @examples
#' dtruncnb(1:4, mu = 2, size = 1.5)
#' ptruncnb(1:4, mu = 2, size = 1.5)
#'
#' @name truncnb
#' @export
dtruncnb <- function(y, mu, size, log.p = FALSE) {
  
  n.y <- max(length(y), length(mu), length(size))
  y <- rep(y, length.out = n.y)
  mu <- rep(mu, length.out = n.y)
  size <- rep(size, length.out = n.y)
  
  log_prob <- rep(-Inf, n.y)
  i.py <- which(y > 0)
  
  log_prob[i.py] <- dnbinom(
    x = y[i.py],
    mu = mu[i.py],
    size = size[i.py],
    log = TRUE
  ) - pnbinom(
    q = 0,
    mu = mu[i.py],
    size = size[i.py],
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (log.p) log_prob else exp(log_prob)
}