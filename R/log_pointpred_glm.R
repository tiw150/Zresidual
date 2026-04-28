# =========================================================================
# Zresidual: Analytic Predictive Functions for stats::glm
# =========================================================================

' Predictive quantities for glm Gaussian models
#'
#' Compute observation-wise predictive quantities for a fitted Gaussian
#' `stats::glm` model on a user-supplied dataset.
#'
#' @param fit A fitted `stats::glm` object with Gaussian family.
#' @param data A data frame on which predictive quantities are evaluated.
#'   This argument must be explicitly provided.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A list with components:
#' \describe{
#'   \item{log_like}{A 1 x n matrix of observation-wise log predictive densities.}
#'   \item{log_surv}{A 1 x n matrix of log survival probabilities,
#'   \eqn{\log P(Y > y_i)}.}
#'   \item{is_discrete}{An integer vector of length n, equal to 0 for all
#'   observations because the Gaussian outcome is continuous.}
#' }
#'
#' @export
log_pointpred_glm_gaussian <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_glm_gaussian: The 'data' argument must be explicitly provided.", call. = FALSE)
  }
  
  mf <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  y_obs <- as.numeric(stats::model.response(mf))
  n <- length(y_obs)
  
  mu_hat <- stats::predict(fit, newdata = data, type = "response")
  
  # In a Gaussian GLM, the dispersion parameter is the variance (sigma^2)
  dispersion <- summary(fit)$dispersion
  sigma_hat  <- sqrt(dispersion)
  
  lpmf <- matrix(stats::dnorm(y_obs, mean = mu_hat, sd = sigma_hat, log = TRUE), nrow = 1)
  lsf  <- matrix(stats::pnorm(y_obs, mean = mu_hat, sd = sigma_hat, lower.tail = FALSE, log.p = TRUE), nrow = 1)
  
  list(
    log_like    = lpmf,
    log_surv    = lsf,
    is_discrete = rep(0L, n)
  )
}

#' Predictive quantities for glm Gamma models
#'
#' Compute observation-wise predictive quantities for a fitted Gamma
#' `stats::glm` model on a user-supplied dataset.
#'
#' @param fit A fitted `stats::glm` object with Gamma family.
#' @param data A data frame on which predictive quantities are evaluated.
#'   This argument must be explicitly provided.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A list with components:
#' \describe{
#'   \item{log_like}{A 1 x n matrix of observation-wise log predictive densities.}
#'   \item{log_surv}{A 1 x n matrix of log survival probabilities,
#'   \eqn{\log P(Y > y_i)}.}
#'   \item{is_discrete}{An integer vector of length n, equal to 0 for all
#'   observations because the Gamma outcome is continuous.}
#' }
#'
#' @export
log_pointpred_glm_gamma <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_glm_Gamma: The 'data' argument must be explicitly provided.", call. = FALSE)
  }
  
  mf <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  y_obs <- as.numeric(stats::model.response(mf))
  n <- length(y_obs)
  
  mu_hat <- stats::predict(fit, newdata = data, type = "response")
  
  # For Gamma GLMs, Variance = dispersion * mu^2.
  # R's dgamma uses 'shape' (alpha) and 'rate' (beta).
  # Shape = 1 / dispersion. Rate = shape / mu.
  dispersion <- summary(fit)$dispersion
  shape_hat  <- 1 / dispersion
  rate_hat   <- shape_hat / mu_hat
  
  lpmf <- matrix(stats::dgamma(y_obs, shape = shape_hat, rate = rate_hat, log = TRUE), nrow = 1)
  lsf  <- matrix(stats::pgamma(y_obs, shape = shape_hat, rate = rate_hat, lower.tail = FALSE, log.p = TRUE), nrow = 1)
  
  list(
    log_like    = lpmf,
    log_surv    = lsf,
    is_discrete = rep(0L, n)
  )
}

#' Predictive quantities for glm Poisson models
#'
#' Compute observation-wise predictive quantities for a fitted Poisson
#' `stats::glm` model on a user-supplied dataset.
#'
#' @param fit A fitted `stats::glm` object with Poisson family.
#' @param data A data frame on which predictive quantities are evaluated.
#'   This argument must be explicitly provided.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A list with components:
#' \describe{
#'   \item{log_like}{A 1 x n matrix of observation-wise log predictive masses.}
#'   \item{log_surv}{A 1 x n matrix of log survival probabilities,
#'   \eqn{\log P(Y > y_i)}.}
#'   \item{is_discrete}{An integer vector of length n, equal to 1 for all
#'   observations because the Poisson outcome is discrete.}
#' }
#'
#' @export
log_pointpred_glm_poisson <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_glm_poisson: The 'data' argument must be explicitly provided.", call. = FALSE)
  }
  
  mf <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  y_obs <- as.numeric(stats::model.response(mf))
  n <- length(y_obs)
  
  lambda_hat <- stats::predict(fit, newdata = data, type = "response")
  
  lpmf <- matrix(stats::dpois(y_obs, lambda = lambda_hat, log = TRUE), nrow = 1)
  lsf  <- matrix(stats::ppois(y_obs, lambda = lambda_hat, lower.tail = FALSE, log.p = TRUE), nrow = 1)
  
  list(
    log_like    = lpmf,
    log_surv    = lsf,
    is_discrete = rep(1L, n)
  )
}


#' Predictive quantities for glm Binomial models
#'
#' Compute observation-wise predictive quantities for a fitted Binomial
#' `stats::glm` model on a user-supplied dataset.
#'
#' @param fit A fitted `stats::glm` object with Binomial family.
#' @param data A data frame on which predictive quantities are evaluated.
#'   This argument must be explicitly provided.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @return A list with components:
#' \describe{
#'   \item{log_like}{A 1 x n matrix of observation-wise log predictive masses.}
#'   \item{log_surv}{A 1 x n matrix of log survival probabilities,
#'   \eqn{\log P(Y > y_i)}.}
#'   \item{is_discrete}{An integer vector of length n, equal to 1 for all
#'   observations because the Binomial outcome is discrete.}
#' }
#'
#' @export
log_pointpred_glm_binomial <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_glm_binomial: The 'data' argument must be explicitly provided.", call. = FALSE)
  }
  
  # Utilize glm's internal parser for edge cases (cbind matrices, factors, proportions)
  mf <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  env <- new.env(parent = environment(stats::formula(fit)))
  env$y <- stats::model.response(mf)
  
  w <- stats::model.weights(mf)
  env$weights <- if (is.null(w)) rep(1, NROW(env$y)) else w
  env$nobs <- NROW(env$y)
  
  eval(stats::family(fit)$initialize, envir = env)
  
  size <- env$weights
  y_obs <- round(env$y * size)
  n <- length(y_obs)
  
  p_hat <- stats::predict(fit, newdata = data, type = "response")
  
  lpmf <- matrix(stats::dbinom(y_obs, size = size, prob = p_hat, log = TRUE), nrow = 1)
  lsf  <- matrix(stats::pbinom(y_obs, size = size, prob = p_hat, lower.tail = FALSE, log.p = TRUE), nrow = 1)
  
  list(
    log_surv    = lsf, 
    log_like    = lpmf,
    is_discrete = rep(1L, n)
  )
}
