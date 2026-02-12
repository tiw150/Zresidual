log_pointpred_brms_bernoulli <- function(fit, data = NULL, ...) {
  
  data_in <- if (is.null(data)) fit$data else data
  
  model.data <- model.frame(fit$formula, data = data_in)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  y <- as.numeric(model.data[[response]])
  O_i <- ifelse(y == 0, 0L, 1L)  # 0/1
  n <- length(O_i)
  
 # y_type <- O_i
  
  posterior_pred_safe <- function(fit, dpar, newdata, ...) {
    f <- get("posterior.pred", mode = "function")
    if ("newdata" %in% names(formals(f))) return(f(fit, dpar = dpar, newdata = newdata, ...))
    f(fit, dpar = dpar, ...)
  }
  
  hu <- posterior_pred_safe(fit, dpar = "mu", newdata = data_in, ...)
  mc_used <- nrow(hu)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n)
  lsf_hat  <- matrix(NA_real_, mc_used, n)
  
  for (i in seq_len(n)) {
    lpmf_hat[, i] <- dbern(O_i[i], hu[, i], log = TRUE)
    lsf_hat[, i]  <- pbern(O_i[i], hu[, i],
                           lower.tail = FALSE, log.p = TRUE)
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = rep.int(1L, n)
  )
}
