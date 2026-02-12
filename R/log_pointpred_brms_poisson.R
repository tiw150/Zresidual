log_pointpred_brms_poisson <- function(fit, data = NULL, ...) {
  
  data_in <- if (is.null(data)) fit$data else data
  
  model.data <- model.frame(fit$formula, data = data_in)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  sim.y <- as.numeric(model.data[[response]])
  n <- length(sim.y)
  
 # y_type <- ifelse(sim.y == 0, 0L, 1L)
  
  posterior_pred_safe <- function(fit, dpar, newdata, ...) {
    f <- get("posterior.pred", mode = "function")
    if ("newdata" %in% names(formals(f))) return(f(fit, dpar = dpar, newdata = newdata, ...))
    f(fit, dpar = dpar, ...)
  }
  
  lambda <- posterior_pred_safe(fit, dpar = "mu", newdata = data_in, ...)
  mc_used <- nrow(lambda)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n)
  lsf_hat  <- matrix(NA_real_, mc_used, n)
  
  for (i in seq_len(n)) {
    lpmf_hat[, i] <- dpois(sim.y[i], lambda = lambda[, i], log = TRUE)
    lsf_hat[, i]  <- ppois(sim.y[i], lambda = lambda[, i],
                           lower.tail = FALSE, log.p = TRUE)
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = rep.int(1L, n)
  )
}
