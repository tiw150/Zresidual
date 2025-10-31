#' Compute Posterior Log Randomized Predictive p-values (Log-RPP)
#'
#' Calculates the posterior log randomized predictive p-values (log-RPP) for a set
#' of observations given the log-PMF and log-CDF matrices from posterior predictive samples.
#'
#' @param log_cdf Numeric matrix of log cumulative probabilities for each observation
#'   (posterior samples × observations).
#' @param log_pmf Numeric matrix of log probability mass function values for each observation
#'   (posterior samples × observations).
#'
#' @details
#' The function computes a randomized posterior predictive p-values as follows:
#' \enumerate{
#'   \item A uniform random number is generated for each observation to randomize the probability for zero counts.
#'   \item Probabilities exactly equal to 0 or 1 are replaced with small bounds (\eqn{10^{-5}} and \eqn{9 \cdot 10^{-5}}) to avoid numerical issues.
#' }
#'
#' @return A numeric vector of log randomized posterior predictive p-values for each observation.
#'
#' @examples
#' \dontrun{
#' # Assume log_cdf and log_pmf are matrices from posterior predictive draws
#' log_rpp <- post_logrpp(log_cdf, log_pmf)
#' head(exp(log_rpp))  # Randomized posterior predictive probabilities
#' }
#'
#' @seealso
#' \code{\link{post_logmpp}}, \code{\link{log.pred.dist.HNB}},
#' \code{\link{log.pred.dist.NB}}, \code{\link{log.pred.dist.TNB}}.
#'
post_logrpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1]
  n <- dim(log_pmf)[2]
  log_u <- log(matrix(runif(n), 1, n, byrow = TRUE))
  Log_Add_Exps <- function (log.cdf, log.pmf)
  {
    lM <- pmax(log.cdf,log.pmf)
    lM + log(exp (log.cdf-lM) + exp(log.pmf-lM))
  }

  #**** Edit: DS 05/10/2025
  log_sum_cdf <- colLogSumExps(log_cdf)
  log_sum_pmf <- colLogSumExps(log_pmf)+log_u
  log_pv<-Log_Add_Exps(log.cdf=log_sum_cdf, log.pmf=log_sum_pmf)
  logrpp_post <- as.vector(log_pv) - log(mc_used)
  # log_mean_exp <- function (lx)
  # {
  #   log_sum_exp (lx) - log(length (lx))
  # }
  # logrpp_post<- apply(log_pv, 2, log_mean_exp)

  id.large <- which(exp(logrpp_post) == 1.000000e+00)
  id.small <- which(exp(logrpp_post) == 0.000000e+00)
  logrpp_post[id.large] <- log(9e-5)
  logrpp_post[id.small] <- log(1e-5)

  logrpp_post
}
