log_pointpred_brms_hurdle_negbinomial_count <- function(fit, data = NULL, ...) {
  
  data_in <- if (is.null(data)) fit$data else data
  
  model.data <- model.frame(fit$formula, data = data_in)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  sim.y <- as.numeric(model.data[[response]])
  count_id <- which(sim.y > 0)
  
  posterior_pred_safe <- function(fit, dpar, newdata, ...) {
    f <- get("posterior.pred", mode = "function")
    if ("newdata" %in% names(formals(f))) return(f(fit, dpar = dpar, newdata = newdata, ...))
    f(fit, dpar = dpar, ...)
  }
  
  mu_all    <- posterior_pred_safe(fit, dpar = "mu",    newdata = data_in, ...)
  shape_all <- posterior_pred_safe(fit, dpar = "shape", newdata = data_in, ...)
  
  mu    <- if (length(count_id) > 0) mu_all[, count_id, drop = FALSE] else mu_all[, 0, drop = FALSE]
  shape <- if (length(count_id) > 0) shape_all[, count_id, drop = FALSE] else shape_all[, 0, drop = FALSE]
  
  y_sub <- sim.y[count_id]
  n_sub <- length(y_sub)
  mc_used <- nrow(mu)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n_sub)
  lsf_hat  <- matrix(NA_real_, mc_used, n_sub)
  
  if (n_sub > 0) {
    for (j in seq_len(n_sub)) {
      lpmf_hat[, j] <- dtruncnb(y_sub[j], mu = mu[, j], size = shape[, j], log.p = TRUE)
      lsf_hat[, j]  <- ptruncnb(y_sub[j], mu = mu[, j], size = shape[, j],
                                lower.tail = FALSE, log.p = TRUE)
    }
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = rep.int(1L, n_sub)  # truncated: only positive part
  )
}
