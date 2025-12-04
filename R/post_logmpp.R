#' Compute Posterior Log Middle-value Predictive p-values (Log-MPP)
#'
#' Calculates the posterior log middle-value (0.5) predictive p-values (log-MPP) for a set
#' of observations given the log-PMF and log-CDF matrices obtained from posterior
#' predictive samples.
#'
#' @param log_cdf Numeric matrix of log cumulative probabilities for each observation
#'   (posterior samples × observations).
#' @param log_pmf Numeric matrix of log probability mass function values for each observation
#'   (posterior samples × observations).
#'
#' @details
#' The function computes a stabilized version of the posterior predictive p-values:
#' \enumerate{
#'   \item It sums the log-CDF and log-PMF across posterior samples using log-sum-exp for numerical stability.
#'   \item Probabilities exactly equal to 0 or 1 are replaced with small bounds (\eqn{10^{-5}} and \eqn{9 \cdot 10^{-5}}) to avoid numerical issues in further computations.
#' }
#'
#' @return A numeric vector of log posterior mid-value predictive probabilities for each observation.
#'
#' @examples NULL
#'

post_logmpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1]
  n <- dim(log_pmf)[2]
  log_u <- log(matrix(0.5, 1, n, byrow = TRUE))
  Log_Add_Exps <- function (log.cdf, log.pmf)
  {
    lM <- pmax(log.cdf,log.pmf)
    lM + log(exp (log.cdf-lM) + exp(log.pmf-lM))
  }

  log_sum_cdf <- colLogSumExps(log_cdf)
  log_sum_pmf <- colLogSumExps(log_pmf)+log_u
  log_pv<-Log_Add_Exps(log.cdf=log_sum_cdf, log.pmf=log_sum_pmf)
  logmpp_post <- as.vector(log_pv) - log(mc_used)
  # log_mean_exp <- function (lx)
  # {
  #   log_sum_exp (lx) - log(length (lx))
  # }
  # logmpp_post<- apply(log_pv, 2, log_mean_exp)

  id.large <- which(exp(logmpp_post) == 1.000000e+00)
  id.small <- which(exp(logmpp_post) == 0.000000e+00)
  logmpp_post[id.large] <- log(9e-5)
  logmpp_post[id.small] <- log(1e-5)

  logmpp_post
}
