#' Compute Log Predictive Distributions for a Hurdle Negative Binomial Model  of a 'brms' Fit
#'
#' This function calculates the log probability mass function (log-PMF) and the log
#' cumulative distribution function (log-CDF) for each observation
#' from a fitted hurdle negative binomial model (fitted using \pkg{brms}).
#' The function extracts posterior samples for the model parameters and evaluates
#' the predictive distributions across all posterior draws.
#'
#' @param fit A fitted \pkg{brms} hurdle negative binomial model object.
#'   The model must include the distributional parameters \code{mu}, \code{shape}, and
#'   the hurdle probability \code{zero}.
#'
#' @details
#' For each posterior draw and each observation, the function computes:
#' \itemize{
#'   \item \code{lpmf_hat}: Log predictive mass function values using \code{dhurdle.nb()}.
#'   \item \code{lcdf_hat}: Log cumulative distribution function values
#'   using \code{phurdle.nb()} with \code{lower.tail = FALSE}.
#' }
#'
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{lpmf_hat}}{A matrix of log-PMF values (posterior samples × observations).}
#'   \item{\code{lcdf_hat}}{A matrix of log-CCDF values (posterior samples × observations).}
#'   \item{\code{zero_id}}{Indices of observations with zero counts.}
#'   \item{\code{count_id}}{Indices of observations with positive counts.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fit <- brm(bf(y ~ x1 + x2, hu ~ x1), family = hurdle_negbinomial(), data = mydata)
#' pred_dist <- log_pred_dist_HNB(fit)
#' str(pred_dist)
#' }
#'

log_pred_dist_HNB <- function(fit){

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

  count_id <- which(sim.y > 0)
  zero_id <- which(sim.y == 0)

  lpmf_hat <- matrix(0, mc_used, n)
  lcdf_hat <- matrix(0, mc_used, n)

  hu <- posterior.pred(fit, dpar = "zero")
  mu <- posterior.pred(fit, dpar = "mu")
  shape <- posterior.pred(fit, dpar = "shape")

  for (i in 1:n){
    #lpmf_hat[,i] <- dhurdle_negbinomial(sim.y[i], mu[,i], shape[,i], hu[,i], log = TRUE)
    lpmf_hat[,i] <- dhurdle.nb(y=sim.y[i], mu=mu[,i], size=shape[,i], pi=hu[,i], log = TRUE)
    lcdf_hat[,i] <- phurdle.nb(sim.y[i], mu[,i], shape[,i], hu[,i], lower.tail=FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id, count_id = count_id)
  return(pred_dist)
}
