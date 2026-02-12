log_pointpred_brms_hurdle_negbinomial <- function(fit, data = NULL, ...) {
  data_in <- if (is.null(data)) fit$data else data
  
  model.data <- model.frame(fit$formula, data = data_in)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  sim.y <- as.numeric(model.data[[response]])
  n <- length(sim.y)
  y_type <- ifelse(sim.y == 0, 0L, 1L)
  
  count_id <- which(sim.y > 0)
  zero_id  <- which(sim.y == 0)
  
  posterior_pred_safe <- function(fit, dpar, newdata) {
    f <- get("posterior.pred", mode = "function")
    if ("newdata" %in% names(formals(f))) return(f(fit, dpar = dpar, newdata = newdata, ...))
    f(fit, dpar = dpar, ...)
  }
  
  hu    <- posterior_pred_safe(fit, dpar = "zero",  newdata = data_in)
  mu    <- posterior_pred_safe(fit, dpar = "mu",    newdata = data_in)
  shape <- posterior_pred_safe(fit, dpar = "shape", newdata = data_in)
  
  mc_used <- nrow(mu)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n)
  lsf_hat  <- matrix(NA_real_, mc_used, n)
  
  for (i in seq_len(n)) {
    lpmf_hat[, i] <- dhurdlenb(y = sim.y[i], mu = mu[, i], size = shape[, i], pi = hu[, i], log = TRUE)
    lsf_hat[, i]  <- phurdlenb(sim.y[i], mu[, i], shape[, i], hu[, i], lower.tail = FALSE, log.p = TRUE)
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = y_type
  )
}
