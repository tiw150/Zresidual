#' A function to calculate log predictive distribution (pmf and cdf) of logistic component of a 'brm' fit
#'
#' @param fit A `brm` fit.
#'
#' @import Rlab

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
