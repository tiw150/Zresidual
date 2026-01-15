#' Cumulative Distribution Function (CDF) of Zero-Truncated Poisson Distribution
#'
#' Computes the cumulative distribution function (CDF) for a zero-truncated Poisson (TP) distribution.
#' This function modifies the standard Poisson CDF to account for truncation at zero (i.e., only supports \eqn{y > 0}).
#'
#' @param y Numeric vector of quantiles (count values) for which to compute the CDF.
#' @param lambda Numeric vector or scalar giving the mean parameter (\eqn{\lambda}) of the Poisson distribution.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(Y \le y)}.
#'   If \code{FALSE}, probabilities are \eqn{P(Y > y)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \eqn{p} are given as \eqn{\log(p)}.
#'
#' @details
#' The function computes probabilities for the zero-truncated version of the Poisson distribution:
#' \deqn{P(Y \le y \mid Y > 0) = \frac{P(Y \le y) - P(Y = 0)}{1 - P(Y = 0)}.}
#'
#' Internally, it uses the \code{ppois()} function from base R for computing Poisson probabilities,
#' and performs calculations on the log scale for improved numerical stability.
#'
#' When \code{lower.tail = FALSE}, it returns the upper-tail probabilities
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
#' \code{\link[stats]{ppois}} for the standard Poisson CDF.
#'
#' @examples
#' # Example: Compute the zero-truncated Poisson CDF for y = 1:5
#' lambda <- 2
#' ptruncpois(1:5, lambda)
#'
#' # Compute upper-tail probabilities
#' ptruncpois(1:5, lambda, lower.tail = FALSE)
#'
#' # Compute log-CDF values
#' ptruncpois(1:5, lambda, log.p = TRUE)
#' @export

ptruncpois<- function(y, lambda, lower.tail = FALSE, log.p = FALSE) {

  n.y <- max(length(y), length(lambda))
  y <- rep(y,length=n.y)
  lambda <- rep(lambda, length=n.y)

  log_upper_prob <- rep (0, n.y)
  i.py <- which (y>=1) ## for y > 0

  ## computing upper tail probabilities for positive y only (y>0)
  log_upper_prob [i.py] <- ppois(q=y[i.py],lambda = lambda[i.py], lower.tail = FALSE, log.p = TRUE) -
    ppois(q=0,lambda = lambda[i.py], lower.tail = FALSE, log.p = TRUE)

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
