#' A function to calculate predictive log poisson pmf and cdf of 'brm' fit.
#'
#' @param fit A `brm` fit.


log.pred.dist.pois <- function(fit){

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
  id <- 1:n
  zero_id <- which(sim.y == 0)
  #if(count_only) id <- which(sim.y > 0) else id <- 1:n
  #if(!count_only) zero_id <- which(sim.y == 0)
  #data_id <- data[count_id,]
  #y_id <- sim.y[count_id]
  #y_id <- as.vector(sim.y)

  lambda <- posterior.pred.ds(fit, dpar = "mu", count.only = F)
  #shape <- posterior.pred.ds(fit, dpar = "shape", count.only = F)

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

  for (i in id){
    lpmf_hat[,i] <- dpois(sim.y[i], lambda = lambda[,i], log = TRUE)
    lcdf_hat[,i] <- ppois(sim.y[i], lambda = lambda[,i], lower.tail = FALSE, log.p = TRUE)
  }

  pred_dist <- list(lpmf_hat = lpmf_hat, lcdf_hat = lcdf_hat, zero_id = zero_id)
  return(pred_dist)
}
