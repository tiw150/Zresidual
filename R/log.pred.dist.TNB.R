#' Compute Log Predictive Distributions for a Truncated Negative Binomial Model of a 'brms' Fit
#'
#' This function calculates the log predictive mass function (log-PMF) and the log
#' cumulative distribution function (log-CDF) for each observation
#' from a fitted truncated negative binomial model (fitted using \pkg{brms}).
#' The function extracts posterior samples for the model’s distributional parameters
#' and evaluates the predictive distributions across all posterior draws.
#'
#' @param fit A fitted \pkg{brms} truncated negative binomial model object.
#'   The model must include the distributional parameters \code{mu} (mean parameter)
#'   and \code{shape} (dispersion parameter).
#'
#' @details
#' For each posterior draw and observation, the function computes:
#' \itemize{
#'   \item \code{lpmf_hat}: Log predictive mass function values using \code{pdf.tnb()}.
#'   \item \code{lcdf_hat}: Log cumulative distribution function values
#'   using \code{cdf.tnb()} with \code{lower.tail = FALSE}.
#' }
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{lpmf_hat}}{A matrix of log-PMF values (posterior samples × observations).}
#'   \item{\code{lcdf_hat}}{A matrix of log-CDF values (posterior samples × observations).}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' fit <- brm(bf(y | trunc(lb = 1) ~ x1 + x2), family = negbinomial(), data = mydata)
#' pred_dist <- log.pred.dist.TNB(fit)
#' str(pred_dist)
#' }
#'
#' @seealso
#' \code{\link{pdf.tnb}}, \code{\link{cdf.tnb}},
#' and \code{\link[brms]{posterior_predict}} for related computations.
#'
log.pred.dist.TNB <- function(fit){

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
  count.id <- which(sim.y > 0)
  zero_id <- which(sim.y == 0)
  #if(count_only) id <- which(sim.y > 0) else id <- 1:n
  #if(!count_only) zero_id <- which(sim.y == 0)
  #data_id <- data[count_id,]
  #y_id <- sim.y[count_id]
  #y_id <- as.vector(sim.y)

  mu <- posterior.pred(fit, dpar = "mu")
  shape <- posterior.pred(fit, dpar = "shape")

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

  for (i in count.id){
    lpmf_hat[,i] <- pdf.tnb(sim.y[i], mu = mu[,i], size = shape[,i], log = TRUE)
    lcdf_hat[,i] <- cdf.tnb(sim.y[i], mu = mu[,i], size = shape[,i], lower.tail = FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id)
  return(pred_dist)
}
