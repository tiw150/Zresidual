#' Probability Mass Function of the Hurdle Negative Binomial Distribution
#'
#' Computes the probability (or log-probability) of observing a count \eqn{y} under
#' a hurdle model with a Bernoulli process for zeros and a truncated Negative Binomial
#' distribution for positive counts (\eqn{y > 0}).
#'
#' @param y Numeric vector of observed count values.
#' @param mu Mean parameter (\eqn{\mu}) of the Negative Binomial component.
#' @param size Dispersion (shape) parameter (\eqn{r}) of the Negative Binomial component.
#' @param pi Probability of observing a structural zero (from the hurdle component).
#' @param log Logical; if \code{TRUE}, probabilities \eqn{p} are returned on the log scale.
#'
#' @details
#' The hurdle Negative Binomial (HNB) distribution models counts \eqn{Y} as:
#'
#' \deqn{
#' P(Y = 0) = \pi,
#' }
#' \deqn{
#' P(Y = y) = (1 - \pi) \frac{f_{NB}(y)}{1 - f_{NB}(0)}, \quad \text{for } y > 0,
#' }
#'
#' where \eqn{f_{NB}(y)} is the Negative Binomial probability mass function with parameters
#' \eqn{\mu} and \eqn{r} (dispersion).
#'
#' Internally, computations are performed on the log scale for numerical stability.
#' The function handles vectorized inputs and returns a numeric vector of the same length as \code{y}.
#'
#' @return
#' A numeric vector of:
#' \itemize{
#'   \item Probabilities \eqn{P(Y = y)} if \code{log = FALSE}.
#'   \item Log-probabilities \eqn{\log P(Y = y)} if \code{log = TRUE}.
#' }
#'
#' @seealso
#' \code{\link[stats]{dnbinom}} for the Negative Binomial PMF,
#' \code{\link[stats]{pnbinom}} for its CDF.
#'
#' @examples
#' # Example parameters
#' mu <- 2
#' size <- 1
#' pi <- 0.3
#'
#' # Compute hurdle NB probabilities for counts 0–5
#' dhurdle.nb(0:5, mu, size, pi)
#'
#' # Log-scale probabilities
#' dhurdle.nb(0:5, mu, size, pi, log = TRUE)
#'
#' # Compare structural zero vs positive counts
#' dhurdle.nb(0, mu, size, pi)      # P(Y = 0)
#' dhurdle.nb(1, mu, size, pi)      # P(Y = 1)
#' @export
dhurdle.nb <- function(y, mu, size, pi, log = FALSE) {
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

  n <- max(length(y), length(mu), length(size), length(pi))
  y <- rep(y, length.out = n)
  mu <- rep(mu, length.out = n)
  size <- rep(size, length.out = n)
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
    log_nb <- dnbinom(y[i_pos], mu = mu[i_pos], size = size[i_pos], log = TRUE)
    log_pnz <- pnbinom(q = 0, mu = mu[i_pos], size = size[i_pos], lower.tail = FALSE, log.p = TRUE)  # P(Y > 0)
    log1m_pi_vals <- log1m_pi(log_pi[i_pos])

    log_dens[i_pos] <- log1m_pi_vals + log_nb - log_pnz
  }

  # Return result
  if (log) {
    return(log_dens)
  } else {
    return(exp(log_dens))
  }
}
