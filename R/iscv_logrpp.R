#' Compute Log Randomized Predictive p-values using ISCV Method
#'
#' Calculates log randomized predictive p-values (log-RPPs) based on
#' the Importance Sampling Cross-Validation (ISCV) method, using
#' precomputed log cumulative distribution function (log-CDF) and
#' log probability mass function (log-PMF) values.
#'
#' @param log_cdf A numeric matrix of log-CDF values with dimensions
#'   \eqn{M \times N}, where \eqn{M} is the number of posterior draws
#'   and \eqn{N} is the number of observations.
#' @param log_pmf A numeric matrix of log-PMF values with the same dimensions
#'   as \code{log_cdf}.
#'
#' @details
#' This function implements the Importance Sampling Cross-Validation (ISCV)
#' version of the log randomized predictive p-value (log-RPP) computation.
#' It uses numerically stable log-sum-exp operations to avoid overflow or
#' underflow when summing in log space.
#'
#' Randomized predictive p-values (\eqn{rpp}) are computed using uniform
#' random draws \eqn{u_i \sim \text{Uniform}(0,1)}, combined with the
#' model-predicted \eqn{\text{CDF}} and \eqn{\text{PMF}} values for
#' each observation and posterior sample. The result is returned in
#' log scale for numerical stability.
#'
#' Internally, this function:
#' \itemize{
#'   \item Generates uniform random values \eqn{u_i} for each observation.
#'   \item Computes \eqn{\log(u_i)} and adds it to \eqn{\log(\text{PMF})}.
#'   \item Applies numerically stable log-sum-exp operations column-wise.
#' }
#'
#' Values numerically equal to 0 or 1 are replaced by \eqn{\log(1e-5)} and
#' \eqn{\log(9e-5)} respectively to keep calculations stable.
#'
#' @return
#' A numeric vector of length \eqn{N}, giving the log randomized predictive
#' p-values for each observation.
#'
#' @examples
#' # Example matrices (toy example)
#' log_cdf <- matrix(log(pnorm(rnorm(20))), 4000, 5)
#' log_pmf <- matrix(log(dnorm(rnorm(20))), 4000, 5)
#' result <- iscv_logrpp(log_cdf, log_pmf)
#' head(result)
#'
#' @seealso
#' \code{\link{colLogSumExps}}, \code{\link{log_sum_exp}},
#' \code{\link{loocv_logrpp}}, \code{\link{posterior_logrpp}}
#'
iscv_logrpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1] # 4000
  n <- dim(log_pmf)[2] # 100. No of observations
  log_u <- log(matrix(runif(n), mc_used, n, byrow = TRUE))

  Log_Add_Exps <- function (log.cdf, log.pmf) # Taking sum of 2 log values without losing precision by overflow/underflow.
  {
    lM <- pmax(log.cdf,log.pmf)
    lM + log(exp (log.cdf-lM) + exp(log.pmf-lM))
  }

  #**** Edit: DS 05/10/2025
  logrpp_iscv <- colLogSumExps(Log_Add_Exps(log.cdf = log_cdf, log.pmf = log_pmf + log_u) - log_pmf) -
    colLogSumExps(-log_pmf)
  # logrpp_iscv<-apply(log_pv - log_pmf, 2, log_sum_exp) -
  #   apply(-log_pmf, 2, log_sum_exp)


  id.large <- which(exp(logrpp_iscv) == 1.000000e+00)
  id.small <- which(exp(logrpp_iscv) == 0.000000e+00)
  logrpp_iscv[id.large] <- log(9e-5)
  logrpp_iscv[id.small] <- log(1e-5)

  logrpp_iscv
}
