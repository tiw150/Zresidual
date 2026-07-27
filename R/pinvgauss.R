#' Inverse Gaussian CDF and survival function
#'
#' Distribution function for the inverse Gaussian distribution using the same
#' parameterization as brms' inverse.gaussian family: mean parameter `mu` and
#' shape parameter `shape`.
#'
#' This function is written to be stable on the log scale, especially for right
#' tail probabilities used by Z-residual calculations.
#'
#' @param q Numeric vector of quantiles.
#' @param mu Positive mean parameter.
#' @param shape Positive shape parameter.
#' @param lower.tail Logical; if `TRUE`, return `P(Y <= q)`, otherwise return
#'   `P(Y > q)`.
#' @param log.p Logical; if `TRUE`, return probabilities on the log scale.
#'
#' @return A numeric vector of CDF/survival values or log-CDF/log-survival values.
#'
#' @export
pinvgauss <- function(q, mu, shape, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(mu), length(shape))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  shape <- rep(shape, length.out = n)

  log_lower <- rep(NA_real_, n)
  log_upper <- rep(NA_real_, n)

  bad_param <- !is.finite(mu) | !is.finite(shape) | mu <= 0 | shape <= 0
  log_lower[bad_param] <- NaN
  log_upper[bad_param] <- NaN

  valid_param <- !bad_param

  le_zero <- valid_param & !is.na(q) & q <= 0
  log_lower[le_zero] <- -Inf
  log_upper[le_zero] <- 0

  pos_inf <- valid_param & is.infinite(q) & q > 0
  log_lower[pos_inf] <- 0
  log_upper[pos_inf] <- -Inf

  ok <- valid_param & is.finite(q) & q > 0

  if (any(ok)) {
    z <- sqrt(shape[ok] / q[ok])
    a <- z * (q[ok] / mu[ok] - 1)
    b <- -z * (q[ok] / mu[ok] + 1)

    # CDF formula:
    # F(q) = Phi(a) + exp(2 * shape / mu) * Phi(b)
    log_phi_a <- stats::pnorm(a, log.p = TRUE)
    log_second <- 2 * shape[ok] / mu[ok] + stats::pnorm(b, log.p = TRUE)
    ll <- .zres_logspace_add2(log_phi_a, log_second)
    ll <- pmin(ll, 0)

    # Survival formula:
    # S(q) = Phi(-a) - exp(2 * shape / mu) * Phi(b)
    log_phi_neg_a <- stats::pnorm(-a, log.p = TRUE)
    lu <- .zres_logspace_subtract(log_phi_neg_a, log_second)

    bad_u <- is.na(lu) | is.nan(lu) | lu > 0
    if (any(bad_u)) {
      # Safe fallback: use log(1 - F) from the log-CDF. This is less precise
      # in an extremely small right tail, but it avoids NA/NaN propagation.
      lu[bad_u] <- .zres_log1m_from_logp(ll[bad_u])
    }

    bad_l <- is.na(ll) | is.nan(ll) | ll > 0
    if (any(bad_l)) {
      ll[bad_l] <- .zres_log1m_from_logp(lu[bad_l])
    }

    log_lower[ok] <- pmin(ll, 0)
    log_upper[ok] <- pmin(lu, 0)
  }

  out <- if (isTRUE(lower.tail)) log_lower else log_upper
  if (isTRUE(log.p)) out else exp(out)
}

.zres_logspace_add2 <- function(a, b) {
  n <- max(length(a), length(b))
  a <- rep(a, length.out = n)
  b <- rep(b, length.out = n)

  out <- rep(NA_real_, n)
  both_neg_inf <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
  out[both_neg_inf] <- -Inf

  ok <- !both_neg_inf & !is.na(a) & !is.na(b)
  if (any(ok)) {
    m <- pmax(a[ok], b[ok])
    out[ok] <- m + log(exp(a[ok] - m) + exp(b[ok] - m))
  }

  out
}

.zres_logspace_subtract <- function(log_a, log_b) {
  n <- max(length(log_a), length(log_b))
  log_a <- rep(log_a, length.out = n)
  log_b <- rep(log_b, length.out = n)

  out <- rep(NaN, n)

  na_idx <- is.na(log_a) | is.na(log_b)
  out[na_idx] <- NA_real_

  b_neg_inf <- !na_idx & is.infinite(log_b) & log_b < 0
  out[b_neg_inf] <- log_a[b_neg_inf]

  a_neg_inf <- !na_idx & is.infinite(log_a) & log_a < 0 & !b_neg_inf
  out[a_neg_inf] <- NaN

  ok <- !na_idx & !b_neg_inf & !a_neg_inf & log_a >= log_b
  if (any(ok)) {
    out[ok] <- log_a[ok] + .zres_log1m_from_logp(log_b[ok] - log_a[ok])
  }

  out
}

.zres_log1m_from_logp <- function(logp) {
  out <- rep(NA_real_, length(logp))

  na_idx <- is.na(logp)
  out[na_idx] <- NA_real_

  too_high <- !na_idx & logp > 0
  out[too_high] <- NaN

  zero_idx <- !na_idx & logp == 0
  out[zero_idx] <- -Inf

  ok <- !na_idx & !too_high & !zero_idx
  if (any(ok)) {
    small <- ok & logp < log(0.5)
    large <- ok & !small

    if (any(small)) {
      out[small] <- log1p(-exp(logp[small]))
    }
    if (any(large)) {
      out[large] <- log(-expm1(logp[large]))
    }
  }

  out
}
