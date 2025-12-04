#' Probability Mass Function of the Hurdle Poisson Distribution
#'
#' Computes the probability density (or log-density) for the Poisson hurdle distribution.
#' This distribution combines a point mass at zero with a truncated-at-zero Poisson
#' distribution for positive counts.
#'
#' @param y Numeric vector of observed counts.
#' @param lambda Numeric vector of Poisson mean parameters (must be positive).
#' @param pi Numeric vector of hurdle probabilities (probability of structural zeros),
#'   where each value must be between 0 and 1.
#' @param log Logical; if \code{TRUE}, probabilities are returned on the log scale.
#'
#' @details
#' The hurdle Poisson distribution assumes:
#' \deqn{
#' P(Y = 0) = \pi
#' }
#' and for \eqn{y > 0}:
#' \deqn{
#' P(Y = y) = (1 - \pi) \frac{P_{\text{Pois}}(Y = y)}{1 - P_{\text{Pois}}(Y = 0)}
#' }
#' where \eqn{P_{\text{Pois}}(Y = y)} is the standard Poisson probability mass function.
#'
#' The function is vectorized over all parameters.
#'
#' @return
#' A numeric vector of the same length as \code{y}, giving the density (or log-density)
#' of the Poisson hurdle distribution.
#'
#' @examples
#' # Example usage:
#' y <- 0:5
#' lambda <- 2
#' pi <- 0.3
#' dhurdle.pois(y, lambda, pi)
#' dhurdle.pois(y, lambda, pi, log = TRUE)
#'
#' @seealso
#' \code{\link{dpois}}, \code{\link{pnbinom}}, \code{\link{dhurdle.nb}}
#' @export
dhurdle.pois <- function(y, lambda, pi, log = FALSE) {
  log1mexp <- function(x) {
    ifelse(x <= log(2),
           log(-expm1(-x)),
           log1p(-exp(-x)))
  }

  # log(1 - pi)
  log1m_pi <- function(log_pi) {
    log_diff_exp(0, log_pi)
  }

  # Log difference of exponentials
  log_diff_exp <- function(la, lb) {
    ifelse(la >= lb, la + log1mexp(la - lb), NaN)
  }

  n <- max(length(y), length(lambda), length(pi))
  y <- rep(y, length.out = n)
  lambda <- rep(lambda, length.out = n)
  pi <- rep(pi, length.out = n)

  log_pi <- log(pi)
  log_dens <- rep(NA_real_, n)

  # y == 0
  i_zero <- which(y == 0)
  if (length(i_zero) > 0) {
    log_dens[i_zero] <- log_pi[i_zero]
  }

  # y > 0
  i_pos <- which(y > 0)
  if (length(i_pos) > 0) {
    log_pois <- dpois(y[i_pos], lambda = lambda[i_pos], log = TRUE)
    log_pnz <- ppois(q = 0, lambda = lambda[i_pos], lower.tail = FALSE, log.p = TRUE)  # P(Y > 0)
    log1m_pi_vals <- log1m_pi(log_pi[i_pos])

    log_dens[i_pos] <- log1m_pi_vals + log_pois - log_pnz
  }

  if (log) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}
