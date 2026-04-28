#' Hurdle negative binomial distribution
#'
#' Density and distribution functions for the hurdle negative binomial
#' distribution.
#'
#' @param y Numeric vector of counts.
#' @param mu Mean parameter of the negative binomial component.
#' @param size Dispersion parameter of the negative binomial component.
#' @param pi Probability of a structural zero.
#' @param log Logical; used only by \code{dhurdlenb()}. If \code{TRUE}, return
#'   log-densities.
#' @param lower.tail Logical; used only by \code{phurdlenb()}. If \code{TRUE},
#'   return \eqn{P(Y \le y)}; otherwise return \eqn{P(Y > y)}.
#' @param log.p Logical; used only by \code{phurdlenb()}. If \code{TRUE}, return
#'   probabilities on the log scale.
#'
#' @return A numeric vector of densities, distribution values, or their
#' logarithms.
#'
#' @examples
#' dhurdlenb(0:3, mu = 2, size = 1.5, pi = 0.3)
#' phurdlenb(0:3, mu = 2, size = 1.5, pi = 0.3)
#'
#' @name hurdlenb
#' @export
dhurdlenb <- function(y, mu, size, pi, log = FALSE) {
  log1m_pi <- function(log_pi) {
    log_diff_exp(0, log_pi)
  }
  
  n <- max(length(y), length(mu), length(size), length(pi))
  y <- rep(y, length.out = n)
  mu <- rep(mu, length.out = n)
  size <- rep(size, length.out = n)
  pi <- rep(pi, length.out = n)
  
  log_pi <- log(pi)
  log_dens <- rep(NA_real_, n)
  
  i_zero <- which(y == 0)
  if (length(i_zero) > 0) {
    log_dens[i_zero] <- log_pi[i_zero]
  }
  
  i_pos <- which(y > 0)
  if (length(i_pos) > 0) {
    log_nb <- dnbinom(y[i_pos], mu = mu[i_pos], size = size[i_pos], log = TRUE)
    log_pnz <- pnbinom(
      q = 0,
      mu = mu[i_pos],
      size = size[i_pos],
      lower.tail = FALSE,
      log.p = TRUE
    )
    log1m_pi_vals <- log1m_pi(log_pi[i_pos])
    
    log_dens[i_pos] <- log1m_pi_vals + log_nb - log_pnz
  }
  
  if (log) {
    log_dens
  } else {
    exp(log_dens)
  }
}