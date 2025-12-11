#' Compute Log Predictive Distributions for a Poisson Model of a 'brms' Fit
#'
#' This function calculates the log predictive mass function (log-PMF) and the log
#' cumulative distribution function (log-CDF) for each observation
#' from a fitted Poisson model (fitted using \pkg{brms}).
#' The function extracts posterior samples for the model’s mean parameter and evaluates
#' the predictive distributions across all posterior draws.
#'
#' @param fit A fitted \pkg{brms} Poisson model object.
#'   The model must include the distributional parameter \code{mu} (mean parameter).
#'
#' @details
#' For each posterior draw and observation, the function computes:
#' \itemize{
#'   \item \code{lpmf_hat}: Log predictive mass function values using \code{dpois()}.
#'   \item \code{lcdf_hat}: Log cumulative distribution function values
#'   using \code{ppois()} with \code{lower.tail = FALSE}.
#' }
#'
#' The function also identifies indices of zero-valued observations.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{lpmf_hat}}{A matrix of log-PMF values (posterior samples × observations).}
#'   \item{\code{lcdf_hat}}{A matrix of log-CDF values (posterior samples × observations).}
#'   \item{\code{zero_id}}{Indices of observations with zero counts.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fit <- brm(bf(y ~ x1 + x2), family = poisson(), data = mydata)
#' pred_dist <- log_pred_dist_pois(fit)
#' str(pred_dist)
#' }
#'

log_pred_dist_pois <- function(fit){

  chains <- summary(fit)$chains
  iter <- summary(fit)$iter
  warmup <- summary(fit)$warmup
  mc_used <- chains*(iter - warmup)

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  id <- 1:n
  zero_id <- which(sim.y == 0)
  #if(count_only) id <- which(sim.y > 0) else id <- 1:n
  #if(!count_only) zero_id <- which(sim.y == 0)
  #data_id <- data[count_id,]
  #y_id <- sim.y[count_id]
  #y_id <- as.vector(sim.y)

  lambda <- posterior.pred(fit, dpar = "mu", count.only = F)
  #shape <- posterior.pred(fit, dpar = "shape", count.only = F)

  #mu <- mu[,count_id]
  #shape <- shape[,count_id]

  # ********* Validating parameters  *********
  # validate_mu <- posterior_linpred(fit, transform = F, dpar = "mu")[, count_id]
  # validate_shape <- posterior_linpred(fit, transform = F, dpar = "shape")[, count_id]
  # if (!all(mu == validate_mu)) warning("Validation failed for mu.")
  # if (!all(shape == validate_shape)) warning("Validation failed for shape.")
  # ******************************************

  lpmf_hat <- matrix(NA, mc_used, n)
  lcdf_hat <- matrix(NA, mc_used, n)

  for (i in id){
    lpmf_hat[,i] <- dpois(sim.y[i], lambda = lambda[,i], log = TRUE)
    lcdf_hat[,i] <- ppois(sim.y[i], lambda = lambda[,i], lower.tail = FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id)
  return(pred_dist)
}
