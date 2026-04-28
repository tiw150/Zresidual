#' Hurdle Poisson distribution
#'
#' Density and distribution functions for the hurdle Poisson distribution.
#'
#' @param y Numeric vector of counts.
#' @param lambda Mean parameter of the Poisson component.
#' @param pi Probability of a structural zero.
#' @param log Logical; used only by \code{dhurdlepois()}. If \code{TRUE}, return
#'   log-densities.
#' @param lower.tail Logical; used only by \code{phurdlepois()}. If \code{TRUE},
#'   return \eqn{P(Y \le y)}; otherwise return \eqn{P(Y > y)}.
#' @param log.p Logical; used only by \code{phurdlepois()}. If \code{TRUE}, return
#'   probabilities on the log scale.
#'
#' @return A numeric vector of densities, distribution values, or their
#' logarithms.
#'
#' @examples
#' dhurdlepois(0:3, lambda = 2, pi = 0.3)
#' phurdlepois(0:3, lambda = 2, pi = 0.3)
#'
#' @name hurdlepois
#' @export
dhurdlepois <- function(y, lambda, pi, log = FALSE) {
  log1m_pi <- function(log_pi) {
    log_diff_exp(0, log_pi)
  }
  
  n <- max(length(y), length(lambda), length(pi))
  y <- rep(y, length.out = n)
  lambda <- rep(lambda, length.out = n)
  pi <- rep(pi, length.out = n)
  
  log_pi <- log(pi)
  log_dens <- rep(NA_real_, n)
  
  i_zero <- which(y == 0)
  if (length(i_zero) > 0) {
    log_dens[i_zero] <- log_pi[i_zero]
  }
  
  i_pos <- which(y > 0)
  if (length(i_pos) > 0) {
    log_pois <- dpois(y[i_pos], lambda = lambda[i_pos], log = TRUE)
    log_pnz <- ppois(
      q = 0,
      lambda = lambda[i_pos],
      lower.tail = FALSE,
      log.p = TRUE
    )
    log1m_pi_vals <- log1m_pi(log_pi[i_pos])
    
    log_dens[i_pos] <- log1m_pi_vals + log_pois - log_pnz
  }
  
  if (log) {
    log_dens
  } else {
    exp(log_dens)
  }
}