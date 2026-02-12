log_pointpred_survival_survreg <- function(fit, data = NULL, ...) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_survreg_survival requires package 'survival'.")
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
    ))
  }
  
  # ---- build model frame / matrix ----
  if (is.null(data)) {
    mf      <- model.frame.survreg(fit)
    fix_var <- model.matrix.survreg(fit)
  } else {
    mf      <- model.frame(fit$terms, data)
    fix_var <- model.matrix(fit$terms, data)
  }
  
  y <- mf[[1]]
  if (!inherits(y, "Surv")) stop("survreg: response is not a Surv object.")
  
  y_mat <- as.matrix(y)
  if (ncol(y_mat) < 2) stop("survreg: Surv object must have (time, status).")
  
  time   <- y_mat[, 1]
  status <- y_mat[, 2]  # 1=event, 0=censored (right-censoring)
  n      <- length(time)
  
  lp    <- as.vector(fix_var %*% fit$coefficients)
  scale <- fit$scale
  parms <- as.numeric(fit[["parms"]])  # only used for t (may be numeric(0))
  
  psurvreg_safe <- function(q, mean, scale, distribution, parms = NULL) {
    fmls <- names(formals(survival::psurvreg))
    args <- list(q, mean = mean, scale = scale, distribution = distribution)
    if (!is.null(parms) && length(parms) > 0 && "parms" %in% fmls) args$parms <- parms
    do.call(survival::psurvreg, args)
  }
  
  dsurvreg_log_safe <- function(x, mean, scale, distribution, parms = NULL) {
    fmls <- names(formals(survival::dsurvreg))
    args <- list(x, mean = mean, scale = scale, distribution = distribution)
    if (!is.null(parms) && length(parms) > 0 && "parms" %in% fmls) args$parms <- parms
    
    if ("log" %in% fmls) {
      args$log <- TRUE
      return(do.call(survival::dsurvreg, args))
    } else {
      dens <- do.call(survival::dsurvreg, args)
      return(log(pmax(dens, .Machine$double.xmin)))
    }
  }
  
  # ---- compute F(t), log f(t) ----
  if (distr %in% supported_main) {
    Ft   <- psurvreg_safe(time, mean = lp, scale = scale, distribution = distr)
    logf <- dsurvreg_log_safe(time, mean = lp, scale = scale, distribution = distr)
  } else { # t
    Ft   <- psurvreg_safe(time, mean = lp, scale = scale, distribution = distr, parms = parms)
    logf <- dsurvreg_log_safe(time, mean = lp, scale = scale, distribution = distr, parms = parms)
  }
  
  Ft   <- pmin(pmax(Ft, 0), 1)
  logS <- log1p(-pmin(Ft, 1 - .Machine$double.eps))
  
  # y_type: 1=event, 0=censored
  y_type <- ifelse(status == 1, 1L, 0L)
  
  list(
    lpmf_hat = matrix(rep(NA_real_, n), nrow = 1, ncol = n),
    lsf_hat  = matrix(logS, nrow = 1, ncol = n),
    y_type   = y_type
  )
}
