log_pointpred_brms_hurdle_poisson_count <- function(fit, data = NULL, ...) {
  
  data_in <- if (is.null(data)) fit$data else data
  
  model.data <- model.frame(fit$formula, data = data_in)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  y_all <- as.numeric(model.data[[response]])
  count_id_orig <- which(y_all > 0)
  
  posterior_pred_safe <- function(fit, dpar, newdata, ...) {
    f <- get("posterior.pred", mode = "function")
    fmls <- names(formals(f))
    
    call_args <- list(fit, dpar = dpar)
    if (!is.null(newdata) && "newdata" %in% fmls) call_args$newdata <- newdata
    
    if ("..." %in% fmls) return(do.call(f, c(call_args, list(...))))
    do.call(f, call_args)
  }
  
  lambda_all <- posterior_pred_safe(fit, dpar = "mu", newdata = data_in, ...)
  mc_used <- nrow(lambda_all)
  
  y <- y_all[count_id_orig]
  lambda <- if (length(count_id_orig) > 0) lambda_all[, count_id_orig, drop = FALSE] else lambda_all[, 0, drop = FALSE]
  n <- length(y)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n)
  lsf_hat  <- matrix(NA_real_, mc_used, n)
  
  if (n > 0) {
    for (j in seq_len(n)) {
      lpmf_hat[, j] <- dtruncpois(y[j], lambda = lambda[, j], log.p = TRUE)
      lsf_hat[, j]  <- ptruncpois(y[j], lambda = lambda[, j], lower.tail = FALSE, log.p = TRUE)
    }
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = rep.int(1L, n)   # truncated: only positive part
  )
}
