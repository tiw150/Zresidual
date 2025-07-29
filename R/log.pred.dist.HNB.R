#' A function to calculate predictive log hurdle negative binomial pmf and cdf of a 'brm' fit
#'
#' @param fit A `brm` fit.
#' @export log.pred.dist.HNB

log.pred.dist.HNB <- function(fit){

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
