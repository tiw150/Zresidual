log_pointpred_brms_hurdle_poisson<- function(fit, data = NULL, ...) {
  
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
  
  posterior_pred_safe <- function(fit, dpar, newdata, ...) {
    f <- get("posterior.pred", mode = "function")
    fmls <- names(formals(f))
    
    call_args <- list(fit, dpar = dpar)
    if (!is.null(newdata) && "newdata" %in% fmls) call_args$newdata <- newdata
    
    if ("..." %in% fmls) {
      return(do.call(f, c(call_args, list(...))))
    }
    do.call(f, call_args)
  }
  
  hu     <- posterior_pred_safe(fit, dpar = "zero", newdata = data_in, ...)
  lambda <- posterior_pred_safe(fit, dpar = "mu",   newdata = data_in, ...)
  
  mc_used <- nrow(lambda)
  
  lpmf_hat <- matrix(NA_real_, mc_used, n)
  lsf_hat  <- matrix(NA_real_, mc_used, n)
  
  for (i in seq_len(n)) {
    lpmf_hat[, i] <- dhurdlepois(sim.y[i], lambda = lambda[, i], pi = hu[, i], log = TRUE)
    lsf_hat[, i]  <- phurdlepois(sim.y[i], lambda = lambda[, i], pi = hu[, i],
                                 lower.tail = FALSE, log.p = TRUE)
  }
  
  list(
    lpmf_hat = lpmf_hat,
    lsf_hat  = lsf_hat,
    y_type   = y_type
  )
}
