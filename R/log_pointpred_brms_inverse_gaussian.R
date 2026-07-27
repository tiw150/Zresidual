.zres_extract_positive_response <- function(fit, data, caller = "inverse.gaussian backend") {
  if (missing(data) || is.null(data)) {
    stop(caller, ": `data` must be provided.", call. = FALSE)
  }

  formula_obj <- tryCatch(fit$formula$formula, error = function(e) NULL)
  if (is.null(formula_obj)) {
    formula_obj <- fit$formula
  }

  mf <- stats::model.frame(formula_obj, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)

  if (inherits(y, "Surv") || is.matrix(y)) {
    stop(caller, ": response must be a univariate positive continuous variable.", call. = FALSE)
  }

  y <- suppressWarnings(as.numeric(y))

  if (length(y) < 1L) {
    stop(caller, ": empty response vector.", call. = FALSE)
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop(caller, ": response contains NA or non-finite values.", call. = FALSE)
  }

  if (any(y <= 0)) {
    stop(caller, ": inverse Gaussian response must be strictly positive.", call. = FALSE)
  }

  y
}

log_pointpred_brmsfit_inverse_gaussian <- function(fit, data, ...) {
  caller <- "log_pointpred_brmsfit_inverse_gaussian"

  y <- .zres_extract_positive_response(fit, data, caller = caller)
  n <- length(y)

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

  mu <- as.matrix(mu)
  shape <- as.matrix(shape)

  if (!identical(dim(mu), dim(shape))) {
    stop(caller, ": `mu` and `shape` must have identical dimensions.", call. = FALSE)
  }

  if (ncol(mu) != n) {
    stop(caller, ": posterior parameter columns must match nrow(data).", call. = FALSE)
  }

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop(caller, ": `mu` must contain only positive finite values.", call. = FALSE)
  }

  if (any(!is.finite(shape)) || any(shape <= 0)) {
    stop(caller, ": `shape` must contain only positive finite values.", call. = FALSE)
  }

  ndraws <- nrow(mu)

  y_mat <- matrix(
    rep(y, each = ndraws),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )

  log_like <- dinvgauss(
    x = y_mat,
    mu = mu,
    shape = shape,
    log = TRUE
  )

  log_surv <- pinvgauss(
    q = y_mat,
    mu = mu,
    shape = shape,
    lower.tail = FALSE,
    log.p = TRUE
  )

  # Important:
  # brms::dinv_gaussian() and brms::pinv_gaussian() may return vectors.
  # Force them back to ndraws x n matrices.
  log_like <- matrix(
    as.numeric(log_like),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )

  log_surv <- matrix(
    as.numeric(log_surv),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )

  if (!identical(dim(log_like), dim(mu))) {
    stop(caller, ": `log_like` dimension mismatch.", call. = FALSE)
  }

  if (!identical(dim(log_surv), dim(mu))) {
    stop(caller, ": `log_surv` dimension mismatch.", call. = FALSE)
  }

  if (anyNA(log_like)) {
    stop(caller, ": `log_like` contains NA.", call. = FALSE)
  }

  if (anyNA(log_surv)) {
    stop(caller, ": `log_surv` contains NA.", call. = FALSE)
  }

  is_discrete <- matrix(0L, nrow = 1L, ncol = n)

  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}


log_pointpred_brms_inverse_gaussian <- log_pointpred_brmsfit_inverse_gaussian
log_pointpred_brmsfit_inverse.gaussian <- log_pointpred_brmsfit_inverse_gaussian
log_pointpred_brms_inverse.gaussian <- log_pointpred_brmsfit_inverse_gaussian