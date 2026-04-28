log_pointpred_brmsfit_poisson <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_poisson: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  response <- fit$formula$resp
  
  if (!(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_poisson: response variable not found in `data`.", call. = FALSE)
  }
  
  sim.y <- as.numeric(model.data[[response]])
  n <- length(sim.y)
  
  if (n < 1L) {
    stop("log_pointpred_brmsfit_poisson: empty response vector.", call. = FALSE)
  }
  
  if (any(!is.finite(sim.y))) {
    stop("log_pointpred_brmsfit_poisson: response contains non-finite values.", call. = FALSE)
  }
  
  if (any(sim.y < 0 | abs(sim.y - round(sim.y)) > 1e-8)) {
    stop("log_pointpred_brmsfit_poisson: response must be non-negative integers.", call. = FALSE)
  }
  
  lambda <- posterior.pred(
    fit,
    dpar = "mu",
    data = data,
    count.only = FALSE,
    ...
  )
  
  if (!is.matrix(lambda)) {
    lambda <- matrix(lambda, nrow = 1L)
  }
  
  if (ncol(lambda) != n) {
    stop("log_pointpred_brmsfit_poisson: posterior parameter columns must match nrow(data).", call. = FALSE)
  }
  
  if (any(!is.finite(lambda))) {
    stop("log_pointpred_brmsfit_poisson: `lambda` contains non-finite values.", call. = FALSE)
  }
  
  if (any(lambda <= 0)) {
    stop("log_pointpred_brmsfit_poisson: `lambda` must be strictly positive.", call. = FALSE)
  }
  
  ndraws <- nrow(lambda)
  
  y_mat <- matrix(
    rep(sim.y, each = ndraws),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )
  
  log_like <- dpois(
    x = y_mat,
    lambda = lambda,
    log = TRUE
  )
  
  log_surv <- ppois(
    q = y_mat,
    lambda = lambda,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (!identical(dim(log_like), dim(lambda))) {
    stop("log_pointpred_brmsfit_poisson: `log_like` dimension mismatch.", call. = FALSE)
  }
  
  if (!identical(dim(log_surv), dim(lambda))) {
    stop("log_pointpred_brmsfit_poisson: `log_surv` dimension mismatch.", call. = FALSE)
  }
  
  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_poisson: `log_like` contains NA.", call. = FALSE)
  }
  
  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_poisson: `log_surv` contains NA.", call. = FALSE)
  }
  
  is_discrete <- matrix(1L, nrow = 1L, ncol = n)
  
  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}