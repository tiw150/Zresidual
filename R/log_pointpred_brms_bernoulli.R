log_pointpred_brmsfit_bernoulli <- function(fit, data, ...) {
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_bernoulli: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  response <- fit$formula$resp
  
  if (!(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_bernoulli: response variable not found in `data`.", call. = FALSE)
  }
  
  y_raw <- model.data[[response]]
  
  if (is.logical(y_raw)) {
    y <- as.integer(y_raw)
  } else if (is.factor(y_raw)) {
    if (nlevels(y_raw) != 2L) {
      stop("log_pointpred_brmsfit_bernoulli: factor response must have exactly 2 levels.", call. = FALSE)
    }
    lev <- levels(y_raw)
    if (all(lev %in% c("0", "1"))) {
      y <- as.integer(as.character(y_raw))
    } else {
      y <- as.integer(y_raw) - 1L
    }
  } else {
    y <- suppressWarnings(as.numeric(y_raw))
  }
  
  n <- length(y)
  
  if (n < 1L) {
    stop("log_pointpred_brmsfit_bernoulli: empty response vector.", call. = FALSE)
  }
  
  if (any(!is.finite(y))) {
    stop("log_pointpred_brmsfit_bernoulli: response contains non-finite values.", call. = FALSE)
  }
  
  if (any(!(y %in% c(0, 1)))) {
    stop("log_pointpred_brmsfit_bernoulli: response must be binary 0/1.", call. = FALSE)
  }
  
  y <- as.integer(y)
  
  p_i <- posterior.pred(
    fit,
    dpar = "mu",
    data = data,
    count.only = FALSE,
    ...
  )
  
  if (!is.matrix(p_i)) {
    p_i <- matrix(p_i, nrow = 1L)
  }
  
  if (ncol(p_i) != n) {
    stop("log_pointpred_brmsfit_bernoulli: posterior parameter columns must match nrow(data).", call. = FALSE)
  }
  
  if (any(!is.finite(p_i))) {
    stop("log_pointpred_brmsfit_bernoulli: probability matrix contains non-finite values.", call. = FALSE)
  }
  
  if (any(p_i < 0 | p_i > 1)) {
    stop("log_pointpred_brmsfit_bernoulli: probabilities must lie in [0, 1].", call. = FALSE)
  }
  
  ndraws <- nrow(p_i)
  
  y_mat <- matrix(
    rep(y, each = ndraws),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )
  
  log_like <- dbinom(
    x = y_mat,
    size = 1L,
    prob = p_i,
    log = TRUE
  )
  
  log_surv <- pbinom(
    q = y_mat,
    size = 1L,
    prob = p_i,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (!identical(dim(log_like), dim(p_i))) {
    stop("log_pointpred_brmsfit_bernoulli: `log_like` dimension mismatch.", call. = FALSE)
  }
  
  if (!identical(dim(log_surv), dim(p_i))) {
    stop("log_pointpred_brmsfit_bernoulli: `log_surv` dimension mismatch.", call. = FALSE)
  }
  
  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_bernoulli: `log_like` contains NA.", call. = FALSE)
  }
  
  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_bernoulli: `log_surv` contains NA.", call. = FALSE)
  }
  
  is_discrete <- matrix(1L, nrow = 1L, ncol = n)
  
  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}