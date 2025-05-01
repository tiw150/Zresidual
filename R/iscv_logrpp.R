#' A function to calculate important sampling cross validation log rpp
#'
#' @param log_cdf log cdf.
#' @param log_pmf log pmf.
#'
#'
iscv_logrpp <- function(log_cdf, log_pmf){
  mc_used <- dim(log_pmf)[1] # 4000
  n <- dim(log_pmf)[2] # 100. No of observations
  log_u <- log(matrix(runif(n), mc_used, n, byrow = TRUE)) # log uniform values 4000*100

  Log_Add_Exps <- function (log.cdf, log.pmf) # Taking sum of 2 log values without losing precision by overflow/underflow.
  {
    lM <- pmax(log.cdf,log.pmf) # finds the element-wise maximum between log.cdf and log.pmf
    lM + log(exp (log_cdf-lM) + exp(log.pmf-lM)) # computes the logarithm of the sum of the exponentials of the
    # shifted values (log.cdf - lM and log.pmf - lM), and adds back lM
    # to get the result.
  }

  log_pv<-Log_Add_Exps(log.cdf=log_cdf, log.pmf=log_pmf+log_u) #Calculating rpp

  logrpp_iscv<-apply(log_pv - log_pmf, 2, log_sum_exp) -
    apply(-log_pmf, 2, log_sum_exp)  # Why this?

  id.large <- which(exp(logrpp_iscv) == 1.000000e+00)
  id.small <- which(exp(logrpp_iscv) == 0.000000e+00)
  logrpp_iscv[id.large] <- log(9e-5)
  logrpp_iscv[id.small] <- log(1e-5)

  logrpp_iscv
}
