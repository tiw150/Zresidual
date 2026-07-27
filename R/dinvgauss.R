#' Inverse Gaussian density
#'
#' Density function for the inverse Gaussian distribution using the same
#' parameterization as brms' inverse.gaussian family: mean parameter `mu` and
#' shape parameter `shape`.
#'
#' @param x Numeric vector of positive values.
#' @param mu Positive mean parameter.
#' @param shape Positive shape parameter.
#' @param log Logical; if `TRUE`, return log-density.
#'
#' @return A numeric vector of density values or log-density values.
#'
#' @export
dinvgauss <- function(x, mu, shape, log = FALSE) {
  n <- max(length(x), length(mu), length(shape))
  x <- rep(x, length.out = n)
  mu <- rep(mu, length.out = n)
  shape <- rep(shape, length.out = n)

  log_dens <- rep(NA_real_, n)

  bad_param <- !is.finite(mu) | !is.finite(shape) | mu <= 0 | shape <= 0
  log_dens[bad_param] <- NaN

  valid_param <- !bad_param
  non_positive <- valid_param & is.finite(x) & x <= 0
  log_dens[non_positive] <- -Inf

  ok <- valid_param & is.finite(x) & x > 0

  if (any(ok)) {
    # log f(x) = 1/2 {log(shape) - log(2*pi) - 3 log(x)}
    #            - shape * (x - mu)^2 / (2 * mu^2 * x)
    # The exponent is written as -0.5 * shape * (x/mu - 1)^2 / x
    # to reduce unnecessary overflow in (x - mu)^2.
    r <- x[ok] / mu[ok] - 1
    log_dens[ok] <- 0.5 * (log(shape[ok]) - log(2 * pi) - 3 * log(x[ok])) -
      0.5 * shape[ok] * r * r / x[ok]
  }

  if (isTRUE(log)) log_dens else exp(log_dens)
}
