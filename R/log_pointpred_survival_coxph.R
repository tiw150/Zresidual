log_pointpred_survival_coxph <- function(fit, data = NULL, ...) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_coxph_survival requires package 'survival'.")
  }
  
  # build model frame / matrix
  if (is.null(data)) {
    mf_new <- model.frame.coxph(fit)
    mm_new <- model.matrix.coxph(fit)
  } else {
    mf_new <- model.frame(fit$terms, data)
    mm_new <- model.matrix.coxph(fit, data)
  }
  
  Y_new <- mf_new[[1]]
  if (!inherits(Y_new, "Surv")) {
    stop("log_pointpred_coxph_survival: response is not a Surv object.")
  }
  
  y_mat <- as.matrix(Y_new)
  if (ncol(y_mat) == 2) {
    time   <- y_mat[, 1]
    status <- y_mat[, 2]  # 1=event, 0=censored
  } else if (ncol(y_mat) == 3) {
    time   <- y_mat[, 2]  # stop time
    status <- y_mat[, 3]  # 1=event, 0=censored
  } else {
    stop("log_pointpred_coxph_survival: unsupported Surv format.")
  }
  n <- length(time)
  
  # baseline cumulative hazard from training fit
  bh <- survival::basehaz(fit, centered = FALSE)
  if (nrow(bh) < 1) stop("log_pointpred_coxph_survival: basehaz returned empty result.")
  
  H0_step <- stats::stepfun(bh$time, c(0, bh$hazard))
  H0_at_t <- H0_step(time)
  
  # linear predictor and survival
  lp  <- as.vector(mm_new %*% fit$coefficients)
  eta <- exp(lp)
  
  SP   <- exp(-eta * H0_at_t)
  logS <- log(pmax(SP, .Machine$double.xmin))
  
  # outputs
  y_type <- ifelse(status == 1, 1L, 0L)
  
  list(
    lpmf_hat = matrix(rep(NA_real_, n), nrow = 1, ncol = n),
    lsf_hat  = matrix(logS, nrow = 1, ncol = n),
    y_type   = y_type
  )
}
