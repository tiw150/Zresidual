log_pointpred_brmsfit_hurdle_poisson_count <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  response <- fit$formula$resp
  
  if (!(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: response variable not found in `data`.", call. = FALSE)
  }
  
  y_all <- as.numeric(model.data[[response]])
  n_all <- length(y_all)
  
  if (n_all < 1L) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: empty response vector.", call. = FALSE)
  }
  
  if (any(!is.finite(y_all))) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: response contains non-finite values.", call. = FALSE)
  }
  
  if (any(y_all < 0 | abs(y_all - round(y_all)) > 1e-8)) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: response must be non-negative integers.", call. = FALSE)
  }
  
  count_id <- which(y_all > 0)
  
  lambda_all <- posterior.pred(
    fit,
    dpar = "mu",
    data = data,
    count.only = FALSE,
    ...
  )
  
  if (!is.matrix(lambda_all)) {
    lambda_all <- matrix(lambda_all, nrow = 1L)
  }
  
  if (ncol(lambda_all) != n_all) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: posterior parameter columns must match nrow(data).", call. = FALSE)
  }
  
  if (any(!is.finite(lambda_all))) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: `lambda` contains non-finite values.", call. = FALSE)
  }
  
  if (any(lambda_all <= 0)) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: `lambda` must be strictly positive.", call. = FALSE)
  }
  
  ndraws <- nrow(lambda_all)
  y_sub <- y_all[count_id]
  n_sub <- length(y_sub)
  
  lambda <- if (n_sub > 0L) {
    lambda_all[, count_id, drop = FALSE]
  } else {
    lambda_all[, 0, drop = FALSE]
  }
  
  log_like <- matrix(NA_real_, nrow = ndraws, ncol = n_sub)
  log_surv <- matrix(NA_real_, nrow = ndraws, ncol = n_sub)
  
  if (n_sub > 0L) {
    for (j in seq_len(n_sub)) {
      log_like[, j] <- dtruncpois(
        y = y_sub[j],
        lambda = lambda[, j],
        log.p = TRUE
      )
      
      log_surv[, j] <- ptruncpois(
        y = y_sub[j],
        lambda = lambda[, j],
        lower.tail = FALSE,
        log.p = TRUE
      )
    }
  }
  
  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: `log_like` contains NA.", call. = FALSE)
  }
  
  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_hurdle_poisson_count: `log_surv` contains NA.", call. = FALSE)
  }
  
  is_discrete <- matrix(1L, nrow = 1L, ncol = n_sub)
  
  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}