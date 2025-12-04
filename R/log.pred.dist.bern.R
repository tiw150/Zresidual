#' Compute Log Predictive Distributions for Logistic Regression of a 'brms' Fit
#'
#' Calculates the log probability mass function (log-PMF) and log
#' cumulative distribution function (log-CDF) for
#' a logistic model based on a fitted 'brms' model.
#'
#' @param fit A fitted model object from \code{brms}.
#'
#' @details
#' The function extracts the posterior predictions of the Bernoulli component
#' (structural zeros) using \code{posterior.pred()} for the mean parameter ("mu").
#'
#' For each observation:
#' \itemize{
#'   \item Computes the log-PMF of the observed binary outcome.
#'   \item Computes the log-CDF (upper-tail probability) of the observed outcome.
#' }
#' This produces matrices of size \eqn{M \times N}, where \eqn{M} is the number of posterior draws
#' and \eqn{N} is the number of observations.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{lpmf_hat}: numeric matrix of log-PMF values (posterior draws × observations).
#'   \item \code{lcdf_hat}: numeric matrix of log-CDF values (posterior draws × observations).
#' }
#'
#' @examples
#' # Assuming 'fit' is a fitted brms logistic model
#' # pred_dist <- log.pred.dist.bern(fit)
#' # lpmf_hat <- pred_dist$lpmf_hat
#' # lcdf_hat <- pred_dist$lcdf_hat
#'

log.pred.dist.bern <- function(fit){

  n <- dim(fit$data)[1]
  chains <- summary(fit)$chains
  iter <- summary(fit)$iter
  warmup <- summary(fit)$warmup
  mc_used <- chains*(iter - warmup)

  data <- fit$data
  data$hu <- 1

  hu <- posterior.pred(fit, dpar = "mu", count.only = F)

  y <- as.matrix(model.frame(fit$formula, data=data)[,1])
  O_i <- ifelse(y==0, 0, 1)

  lpmf_hat <- matrix(0, mc_used, n)
  lcdf_hat <- matrix(0, mc_used, n)

  for (i in 1:n){
    lpmf_hat[,i] <- dbern(O_i[i], hu[,i], log = TRUE) # 4000 x 100
    lcdf_hat[,i] <- pbern(O_i[i], hu[,i], lower.tail=FALSE, log.p = TRUE) # 4000 x 100
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat)
  return(pred_dist)
}
