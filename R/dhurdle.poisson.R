#' A function to calculate pdf of hurdle poisson.
#'
#' @param y y values.
#' @param mu mu parameter.
#' @param pi hu parameter.
#'
dhurdle.poisson <- function(y, mu, pi, log = FALSE) {
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

  n <- max(length(y), length(mu), length(pi))
  y <- rep(y, length.out = n)
  mu <- rep(mu, length.out = n)
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
    log_pois <- dpois(y[i_pos], lambda = mu[i_pos], log = TRUE)
    log_pnz <- ppois(q = 0, lambda = mu[i_pos], lower.tail = FALSE, log.p = TRUE)  # P(Y > 0)
    log1m_pi_vals <- log1m_pi(log_pi[i_pos])

    log_dens[i_pos] <- log1m_pi_vals + log_pois - log_pnz
  }

  # Return result
  if (log) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}
