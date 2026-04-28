log_pointpred_survival_survreg <- function(fit, data, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_survival_survreg requires package 'survival'.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_survival_survreg: `data` must be provided.", call. = FALSE)
  }
  
  supported_main <- c(
    "weibull", "exponential", "logistic", "lognormal",
    "loglogistic", "gaussian", "loggaussian", "rayleigh"
  )
  supported_t <- "t"
  
  distr <- fit$dist
  if (!(distr %in% c(supported_main, supported_t))) {
    stop(sprintf(
      "survreg distribution '%s' not supported. Supported: %s",
      distr, paste(c(supported_main, supported_t), collapse = ", ")
    ), call. = FALSE)
  }
  
  mf <- stats::model.frame(fit$terms, data, xlev = fit$xlevels)
  X <- stats::model.matrix(fit$terms, mf, contrasts.arg = fit$contrasts)
  
  y <- mf[[1]]
  if (!inherits(y, "Surv")) {
    stop("log_pointpred_survival_survreg: response is not a Surv object.", call. = FALSE)
  }
  
  if (!identical(attr(y, "type"), "right")) {
    stop("log_pointpred_survival_survreg: only right-censoring is supported.", call. = FALSE)
  }
  
  y_mat <- as.matrix(y)
  if (ncol(y_mat) != 2L) {
    stop("log_pointpred_survival_survreg: Surv must have exactly (time, status).", call. = FALSE)
  }
  
  time <- y_mat[, 1]
  status <- as.integer(y_mat[, 2])
  n <- length(time)
  
  lp <- as.vector(X %*% fit$coefficients)
  off <- stats::model.offset(mf)
  if (!is.null(off)) {
    lp <- lp + off
  }
  
  scale <- fit$scale
  parms <- as.numeric(fit[["parms"]])
  
  psurvreg_safe <- function(q, mean, scale, distribution, parms = NULL) {
    fmls <- names(formals(survival::psurvreg))
    args <- list(q = q, mean = mean, scale = scale, distribution = distribution)
    if (!is.null(parms) && length(parms) > 0 && "parms" %in% fmls) {
      args$parms <- parms
    }
    do.call(survival::psurvreg, args)
  }
  
  dsurvreg_safe <- function(x, mean, scale, distribution, parms = NULL) {
    fmls <- names(formals(survival::dsurvreg))
    args <- list(x = x, mean = mean, scale = scale, distribution = distribution)
    if (!is.null(parms) && length(parms) > 0 && "parms" %in% fmls) {
      args$parms <- parms
    }
    do.call(survival::dsurvreg, args)
  }
  
  if (distr %in% supported_main) {
    Ft <- psurvreg_safe(time, mean = lp, scale = scale, distribution = distr)
    ft <- dsurvreg_safe(time, mean = lp, scale = scale, distribution = distr)
  } else {
    Ft <- psurvreg_safe(time, mean = lp, scale = scale, distribution = distr, parms = parms)
    ft <- dsurvreg_safe(time, mean = lp, scale = scale, distribution = distr, parms = parms)
  }
  
  Ft <- pmin(pmax(Ft, 0), 1)
  ft <- pmax(ft, 0)
  
  logS <- log1p(-pmin(Ft, 1 - .Machine$double.eps))
  logf <- ifelse(ft > 0, log(ft), -Inf)
  
  is_censored <- (status == 0L)
  
  log_surv <- logS
  log_surv[is_censored] <- -Inf
  
  log_like <- logf
  log_like[is_censored] <- logS[is_censored]
  
  is_discrete <- as.integer(is_censored)
  
  list(
    log_surv = matrix(log_surv, nrow = 1L, ncol = n),
    log_like = matrix(log_like, nrow = 1L, ncol = n),
    is_discrete = matrix(is_discrete, nrow = 1L, ncol = n)
  )
}
