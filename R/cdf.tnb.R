#' Cumulative Distribution Function (CDF) of Zero-Truncated Negative Binomial Distribution
#'
#' Computes the cumulative distribution function (CDF) for a zero-truncated negative binomial (TNB) distribution.
#' This function adjusts the standard negative binomial CDF to account for truncation at zero (i.e., only supports \eqn{y > 0}).
#'
#' @param y Numeric vector of quantiles (count values) for which to compute the CDF.
#' @param mu Mean parameter (\eqn{\mu}) of the negative binomial distribution.
#' @param size Dispersion (shape) parameter (\eqn{r}) of the negative binomial distribution.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le y)}.
#'   If \code{FALSE}, probabilities are \eqn{P(Y > y)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#'
#' @details
#' The function computes probabilities for the zero-truncated version of the negative binomial distribution:
#' \deqn{P(Y \le y \mid Y > 0) = \frac{P(Y \le y) - P(Y = 0)}{1 - P(Y = 0)}.}
#'
#' Internally, this is implemented using the log-scale for numerical stability.
#' When \code{lower.tail = FALSE}, it computes the upper-tail probabilities
#' \eqn{P(Y > y \mid Y > 0)} instead.
#'
#' @return
#' A numeric vector of the same length as the input, containing:
#' \itemize{
#'   \item CDF values (\eqn{P(Y \le y \mid Y > 0)}) if \code{lower.tail = TRUE}.
#'   \item Upper-tail probabilities (\eqn{P(Y > y \mid Y > 0)}) if \code{lower.tail = FALSE}.
#'   \item Log-probabilities if \code{log.p = TRUE}.
#' }
#'
#' @seealso
#' \code{\link[stats]{pnbinom}} for the standard negative binomial CDF.
#'
#' @examples
#' # Example: Compute the CDF for y = 1:5
#' mu <- 2
#' size <- 1
#' p_tnb(1:5, mu, size)
#'
#' # Compute the upper-tail probabilities
#' p_tnb(1:5, mu, size, lower.tail = FALSE)
#'
#' # Log probabilities
#' p_tnb(1:5, mu, size, log.p = TRUE)
#' @export

p_tnb <- function(y, mu, size, lower.tail = FALSE, log.p = FALSE) {

  n.y <- max(length(y), length(mu),length(size),length(pi))
  y <- rep(y,length=n.y)
  mu <- rep(mu, length=n.y)
  size <- rep(size, length=n.y)
  # pi <- rep(pi, length=n.y)

  log_upper_prob <- rep (0, n.y)
  #if(count_only) i.py <- which(y>=1) else i.py <- 1:n

  i.py <- which (y>=1) ## for y > 0

  ## computing upper tail probabilities for positive y only (y>0)
  log_upper_prob [i.py] <- pnbinom(q=y[i.py],mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE) -
    pnbinom(q=0,mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE)

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
