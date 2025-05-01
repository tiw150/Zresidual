#' A function to calculate predictive log hurdle poisson pmf and cdf of 'brm' fit.
#'
#' @import brms
#' @param fit A `brm` fit.
#'

log.pred.dist.HP <- function(fit){

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

  hu <- posterior.pred.ds(fit, dpar = "zero")
  lambda <- posterior.pred.ds(fit, dpar = "mu")

  for (i in 1:n){
    lpmf_hat[,i] <- dhurdle_poisson(sim.y[i], lambda = lambda[,i], hu = hu[,i], log = TRUE)
    lcdf_hat[,i] <- phurdle.pois.li(sim.y[i], lambda = lambda[,i], pi = hu[,i], lower.tail=FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id, count_id = count_id)
  return(pred_dist)
}
