log_pointpred_brmsfit_hurdle_negbinomial <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  response <- fit$formula$resp
  
  if (!(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: response variable not found in `data`.", call. = FALSE)
  }
  
  sim.y <- as.numeric(model.data[[response]])
  n <- length(sim.y)
  
  if (n < 1L) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: empty response vector.", call. = FALSE)
  }
  
  if (any(!is.finite(sim.y))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: response contains non-finite values.", call. = FALSE)
  }
  
  if (any(sim.y < 0 | abs(sim.y - round(sim.y)) > 1e-8)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: response must be non-negative integers.", call. = FALSE)
  }
  
  hu <- posterior.pred(
    fit,
    dpar = "zero",
    data = data,
    count.only = FALSE,
    ...
  )
  
  mu <- posterior.pred(
    fit,
    dpar = "mu",
    data = data,
    count.only = FALSE,
    ...
  )
  
  shape <- posterior.pred(
    fit,
    dpar = "shape",
    data = data,
    count.only = FALSE,
    ...
  )
  
  if (!is.matrix(hu)) {
    hu <- matrix(hu, nrow = 1L)
  }
  if (!is.matrix(mu)) {
    mu <- matrix(mu, nrow = 1L)
  }
  if (!is.matrix(shape)) {
    shape <- matrix(shape, nrow = 1L)
  }
  
  if (!identical(dim(hu), dim(mu)) || !identical(dim(mu), dim(shape))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `hu`, `mu`, and `shape` must have identical dimensions.", call. = FALSE)
  }
  
  if (ncol(mu) != n) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: posterior parameter columns must match nrow(data).", call. = FALSE)
  }
  
  if (any(!is.finite(hu))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `hu` contains non-finite values.", call. = FALSE)
  }
  
  if (any(hu < 0 | hu > 1)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `hu` must lie in [0, 1].", call. = FALSE)
  }
  
  if (any(!is.finite(mu))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `mu` contains non-finite values.", call. = FALSE)
  }
  
  if (any(mu <= 0)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `mu` must be strictly positive.", call. = FALSE)
  }
  
  if (any(!is.finite(shape))) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `shape` contains non-finite values.", call. = FALSE)
  }
  
  if (any(shape <= 0)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `shape` must be strictly positive.", call. = FALSE)
  }
  
  ndraws <- nrow(mu)
  
  log_like <- matrix(NA_real_, nrow = ndraws, ncol = n)
  log_surv <- matrix(NA_real_, nrow = ndraws, ncol = n)
  
  for (i in seq_len(n)) {
    log_like[, i] <- dhurdlenb(
      y = sim.y[i],
      mu = mu[, i],
      size = shape[, i],
      pi = hu[, i],
      log = TRUE
    )
    
    log_surv[, i] <- phurdlenb(
      y = sim.y[i],
      mu = mu[, i],
      size = shape[, i],
      pi = hu[, i],
      lower.tail = FALSE,
      log.p = TRUE
    )
  }
  
  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `log_like` contains NA.", call. = FALSE)
  }
  
  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_hurdle_negbinomial: `log_surv` contains NA.", call. = FALSE)
  }
  
  is_discrete <- matrix(1L, nrow = 1L, ncol = n)
  
  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}