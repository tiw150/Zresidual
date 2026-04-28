log_pointpred_survival_coxph <- function(fit, data, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_survival_coxph requires package 'survival'.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_survival_coxph: `data` must be provided.", call. = FALSE)
  }
  
  mf_new <- stats::model.frame(fit, data = data)
  mm_new <- stats::model.matrix(fit, data = mf_new)
  
  Y_new <- mf_new[[1]]
  if (!inherits(Y_new, "Surv")) {
    stop("log_pointpred_survival_coxph: response is not a Surv object.", call. = FALSE)
  }
  
  y_mat <- as.matrix(Y_new)
  if (ncol(y_mat) == 2L) {
    time <- y_mat[, 1]
    status <- as.integer(y_mat[, 2])
  } else if (ncol(y_mat) == 3L) {
    time <- y_mat[, 2]
    status <- as.integer(y_mat[, 3])
  } else {
    stop("log_pointpred_survival_coxph: unsupported Surv format.", call. = FALSE)
  }
  
  n <- length(time)
  
  bh <- survival::basehaz(fit, centered = FALSE)
  if (nrow(bh) < 1L) {
    stop("log_pointpred_survival_coxph: basehaz returned empty result.", call. = FALSE)
  }
  
  if ("strata" %in% names(bh)) {
    stop("log_pointpred_survival_coxph: stratified coxph is not supported.", call. = FALSE)
  }
  
  lp <- as.vector(mm_new %*% fit$coefficients)
  off <- stats::model.offset(mf_new)
  if (!is.null(off)) {
    lp <- lp + off
  }
  
  eta <- exp(lp)
  
  # --- Robust Interval Lookup ---
  idx_after <- findInterval(time, bh$time)
  exact_idx <- rep(0L, n)
  tol <- 1e-7
  
  # Check for matching hazard jumps using a tolerance to avoid floating-point drift
  for (i in seq_len(n)) {
    if (idx_after[i] > 0L && abs(bh$time[idx_after[i]] - time[i]) <= tol) {
      exact_idx[i] <- idx_after[i]
    } else if (idx_after[i] < nrow(bh) && abs(bh$time[idx_after[i] + 1L] - time[i]) <= tol) {
      idx_after[i] <- idx_after[i] + 1L
      exact_idx[i] <- idx_after[i]
    }
  }
  
  H_after <- ifelse(idx_after > 0L, bh$hazard[idx_after], 0)
  
  # Calculate before-jump hazard (only drops down if it was an exact match)
  idx_before <- idx_after
  exact_match <- exact_idx > 0L
  idx_before[exact_match] <- exact_idx[exact_match] - 1L
  H_before <- ifelse(idx_before > 0L, bh$hazard[idx_before], 0)
  
  logS_after <- -eta * H_after
  logS_before <- -eta * H_before
  
  is_event <- (status == 1L)
  
  # --- Output Construction ---
  
  # 1. Predictive: log_surv
  log_surv <- logS_after
  log_surv[!is_event] <- -Inf
  
  # 2. Generative: log_like 
  log_like <- numeric(n)
  
  # For censored: The observation is discrete (Y > C).
  # The log-likelihood is the probability mass, which is exactly the log-survival tail.
  log_like[!is_event] <- logS_after[!is_event]
  
  # Internal helper for numerical stability: log(exp(a) - exp(b))
  local_log_diff_exp <- function(a, b) {
    diff_ab <- b - a
    res <- a + log1p(-exp(diff_ab))
    res[diff_ab == 0] <- -Inf # If S_before == S_after, density is 0, log-density is -Inf
    return(res)
  }
  
  # For events: The observation is continuous (Y = t).
  # The log-likelihood is the continuous density f(t).
  if (any(is_event)) {
    log_like[is_event] <- local_log_diff_exp(logS_before[is_event], logS_after[is_event])
  }
  
  is_discrete <- as.integer(!is_event)
  
  list(
    log_surv = matrix(log_surv, nrow = 1L, ncol = n),
    log_like = matrix(log_like, nrow = 1L, ncol = n),
    is_discrete = matrix(is_discrete, nrow = 1L, ncol = n)
  )
}


# Tingxuan's old code. I made minor change to avoid an exact event checking.
# log_pointpred_survival_coxph <- function(fit, data, ...) {
#   if (!requireNamespace("survival", quietly = TRUE)) {
#     stop("log_pointpred_survival_coxph requires package 'survival'.", call. = FALSE)
#   }
  
#   if (missing(data) || is.null(data)) {
#     stop("log_pointpred_survival_coxph: `data` must be provided.", call. = FALSE)
#   }
  
#   mf_new <- model.frame.coxph(fit, data = data)
#   mm_new <- model.matrix.coxph(fit, data = mf_new)
  
#   Y_new <- mf_new[[1]]
#   if (!inherits(Y_new, "Surv")) {
#     stop("log_pointpred_survival_coxph: response is not a Surv object.", call. = FALSE)
#   }
  
#   y_mat <- as.matrix(Y_new)
#   if (ncol(y_mat) == 2L) {
#     time <- y_mat[, 1]
#     status <- as.integer(y_mat[, 2])
#   } else if (ncol(y_mat) == 3L) {
#     time <- y_mat[, 2]
#     status <- as.integer(y_mat[, 3])
#   } else {
#     stop("log_pointpred_survival_coxph: unsupported Surv format.", call. = FALSE)
#   }
  
#   n <- length(time)
  
#   bh <- survival::basehaz(fit, centered = FALSE)
#   if (nrow(bh) < 1L) {
#     stop("log_pointpred_survival_coxph: basehaz returned empty result.", call. = FALSE)
#   }
  
#   if ("strata" %in% names(bh)) {
#     stop("log_pointpred_survival_coxph: stratified coxph is not supported.", call. = FALSE)
#   }
  
#   lp <- as.vector(mm_new %*% fit$coefficients)
#   off <- stats::model.offset(mf_new)
#   if (!is.null(off)) {
#     lp <- lp + off
#   }
  
#   eta <- exp(lp)
  
#   idx_after <- findInterval(time, bh$time)
#   H_after <- ifelse(idx_after > 0L, bh$hazard[idx_after], 0)
  
#   is_event <- (status == 1L)
#   exact_event <- is_event & idx_after > 0L & bh$time[idx_after] == time
#   if (any(is_event & !exact_event)) {
#     stop("log_pointpred_survival_coxph: some event times were not found in the baseline hazard grid.", call. = FALSE)
#   }
  
#   idx_before <- idx_after
#   idx_before[exact_event] <- idx_before[exact_event] - 1L
#   H_before <- ifelse(idx_before > 0L, bh$hazard[idx_before], 0)
  
#   logS_after <- -eta * H_after
#   logS_before <- -eta * H_before
  
#   log_surv <- logS_after
#   log_surv[!is_event] <- -Inf
  
#   log_like <- logS_after
#   if (any(is_event)) {
#     log_like[is_event] <- log_diff_exp(logS_before[is_event], logS_after[is_event])
#   }
  
#   is_discrete <- as.integer(!is_event)
  
#   list(
#     log_surv = matrix(log_surv, nrow = 1L, ncol = n),
#     log_like = matrix(log_like, nrow = 1L, ncol = n),
#     is_discrete = matrix(is_discrete, nrow = 1L, ncol = n)
#   )
# }
