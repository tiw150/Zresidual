log_pointpred_brmsfit_hurdle_zero <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_hurdle_zero: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  response <- fit$formula$resp
  
  if (!(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_hurdle_zero: response variable not found in `data`.", call. = FALSE)
  }
  
  y <- as.numeric(model.data[[response]])
  n <- length(y)
  
  if (n < 1L) {
    stop("log_pointpred_brmsfit_hurdle_zero: empty response vector.", call. = FALSE)
  }
  
  if (any(!is.finite(y))) {
    stop("log_pointpred_brmsfit_hurdle_zero: response contains non-finite values.", call. = FALSE)
  }
  
  if (any(y < 0 | abs(y - round(y)) > 1e-8)) {
    stop("log_pointpred_brmsfit_hurdle_zero: underlying response must be non-negative integers.", call. = FALSE)
  }
  
  O_i <- as.integer(y == 0)
  
  hu <- posterior.pred(
    fit,
    dpar = "zero",
    data = data,
    count.only = FALSE,
    ...
  )
  
  if (!is.matrix(hu)) {
    hu <- matrix(hu, nrow = 1L)
  }
  
  if (ncol(hu) != n) {
    stop("log_pointpred_brmsfit_hurdle_zero: posterior parameter columns must match nrow(data).", call. = FALSE)
  }
  
  if (any(!is.finite(hu))) {
    stop("log_pointpred_brmsfit_hurdle_zero: `hu` contains non-finite values.", call. = FALSE)
  }
  
  if (any(hu < 0 | hu > 1)) {
    stop("log_pointpred_brmsfit_hurdle_zero: `hu` must lie in [0, 1].", call. = FALSE)
  }
  
  ndraws <- nrow(hu)
  
  o_mat <- matrix(
    rep(O_i, each = ndraws),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )
  
  log_like <- dbinom(
    x = o_mat,
    size = 1L,
    prob = hu,
    log = TRUE
  )
  
  log_surv <- pbinom(
    q = o_mat,
    size = 1L,
    prob = hu,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (!identical(dim(log_like), dim(hu))) {
    stop("log_pointpred_brmsfit_hurdle_zero: `log_like` dimension mismatch.", call. = FALSE)
  }
  
  if (!identical(dim(log_surv), dim(hu))) {
    stop("log_pointpred_brmsfit_hurdle_zero: `log_surv` dimension mismatch.", call. = FALSE)
  }
  
  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_hurdle_zero: `log_like` contains NA.", call. = FALSE)
  }
  
  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_hurdle_zero: `log_surv` contains NA.", call. = FALSE)
  }
  
  is_discrete <- matrix(1L, nrow = 1L, ncol = n)
  
  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}

log_pointpred_brmsfit_hurdle_negbinomial_zero <- log_pointpred_brmsfit_hurdle_zero
log_pointpred_brmsfit_hurdle_poisson_zero <- log_pointpred_brmsfit_hurdle_zero
