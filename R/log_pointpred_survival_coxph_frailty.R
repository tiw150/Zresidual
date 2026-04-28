#' Predictive quantities for survival::coxph with frailty
#' 
#' @param fit A fitted coxph object (usually class coxph.penal).
#' @param data The new data to evaluate.
#' @param traindata The original data used to fit the model (required for baseline hazard).
#' @param ... Additional arguments passed from Zresidual.
#' @export


log_pointpred_survival_coxph.penal  <- function(fit, data, traindata, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_survival_coxph_frailty requires package 'survival'.", call. = FALSE)
  }
  
  # Check if both are missing or NULL
  if ((missing(data) || is.null(data)) && (missing(traindata) || is.null(traindata))) {
    stop("log_pointpred_survival_coxph_frailty: At least one of `data` or `traindata` must be provided.", call. = FALSE)
  }

  # If only one is missing, set it equal to the other
  if (missing(data) || is.null(data)) {
    data <- traindata
  }

  if (missing(traindata) || is.null(traindata)) {
    traindata <- data
  }
  is_surv_mf <- function(x) {
    is.data.frame(x) && ncol(x) >= 1L && inherits(x[[1]], "Surv") && !is.null(attr(x, "terms"))
  }
  
  make_mf <- function(dat) {
    if (is_surv_mf(dat)) return(dat)
    stats::model.frame(fit$terms, dat, xlev = fit$xlevels)
  }
  
  frailty_group_var <- function(fit) {
    tl <- attr(stats::terms(fit), "term.labels")
    lab <- tl[grepl("^frailty\\(", tl)]
    if (length(lab) != 1L) {
      stop("coxph_frailty: expect exactly one frailty() term.", call. = FALSE)
    }
    m <- regexec("^frailty\\(([^,\\)]+)", lab)
    g <- regmatches(lab, m)[[1]]
    if (length(g) < 2L) {
      stop("coxph_frailty: cannot parse frailty() group variable.", call. = FALSE)
    }
    trimws(g[2])
  }
  
  group_var_name <- frailty_group_var(fit)
  
  get_group <- function(dat, name, where) {
    if (!is.data.frame(dat)) {
      stop(sprintf("coxph_frailty: %s must be a data.frame.", where), call. = FALSE)
    }
    if (!name %in% names(dat)) {
      stop(sprintf("coxph_frailty: %s missing group variable '%s'.", where, name), call. = FALSE)
    }
    dat[[name]]
  }
  
  group_train_raw <- get_group(traindata, group_var_name, "traindata")
  group_new_raw <- get_group(data, group_var_name, "data")
  
  frail_levels <- NULL
  if (!is.null(fit$frail) && length(fit$frail) && !is.null(names(fit$frail)) && all(nzchar(names(fit$frail)))) {
    frail_levels <- names(fit$frail)
  }
  
  if (is.null(frail_levels)) {
    group_train <- factor(group_train_raw)
  } else {
    group_train <- factor(as.character(group_train_raw), levels = frail_levels)
  }
  
  if (anyNA(group_train)) {
    stop("coxph_frailty: `traindata` contains group levels not present in `fit$frail`.", call. = FALSE)
  }
  
  group_new <- factor(as.character(group_new_raw), levels = levels(group_train))
  
  traindata[[group_var_name]] <- group_train
  data[[group_var_name]] <- group_new
  
  mf_train <- make_mf(traindata)
  mm_train <- model.matrix.coxph(fit, data = mf_train)
  
  mf_new <- make_mf(data)
  mm_new <- model.matrix.coxph(fit, data = mf_new)
  
  as_counting_surv <- function(Y) {
    if (!inherits(Y, "Surv")) {
      stop("coxph_frailty: response is not a Surv object.", call. = FALSE)
    }
    Ymat <- as.matrix(Y)
    if (ncol(Ymat) != 3L) {
      return(as.matrix(survival::Surv(rep(0, nrow(Ymat)), Ymat[, 1], Ymat[, 2])))
    }
    Ymat
  }
  
  getchz <- function(Y, explp) {
    death <- (Y[, ncol(Y)] == 1)
    dtime <- Y[, ncol(Y) - 1]
    
    time <- sort(unique(dtime))
    if (length(time) < 1L) {
      return(data.frame(time = numeric(0), cumhaz = numeric(0)))
    }
    
    nevent <- as.vector(rowsum(1 * death, dtime))
    nrisk <- rev(cumsum(rev(rowsum(explp, dtime))))
    
    delta <- if (length(time) >= 2L) min(diff(time)) / 2 else 0.5
    etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)
    
    indx <- approx(etime, 1:length(etime), time, method = "constant", rule = 2, f = 1)$y
    
    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
    nrisk <- nrisk - c(esum, 0)[indx]
    nrisk <- pmax(nrisk, .Machine$double.xmin)
    
    cumhaz <- cumsum(nevent / nrisk)
    data.frame(time = time, cumhaz = cumhaz)
  }
  
  mm_nc <- ncol(mm_train)
  if (mm_nc != (length(fit$coefficients) + 1L)) {
    stop("coxph_frailty: cannot identify the frailty group-code column in the design matrix.", call. = FALSE)
  }
  
  fix_var_train <- mm_train[, -mm_nc, drop = FALSE]
  beta <- fit$coefficients
  
  if (!is.null(names(beta)) && !is.null(colnames(fix_var_train))) {
    idx <- match(colnames(fix_var_train), names(beta))
    if (anyNA(idx)) {
      stop("coxph_frailty: coefficient names do not match training design matrix columns.", call. = FALSE)
    }
    beta_use <- beta[idx]
  } else {
    beta_use <- beta
  }
  
  lp_fixed_train <- as.vector(fix_var_train %*% beta_use)
  frail_train <- fit$frail[as.integer(group_train)]
  frail_train[is.na(frail_train)] <- 0
  explp_train <- exp(lp_fixed_train + frail_train)
  
  Y_train <- as_counting_surv(mf_train[[1]])
  df_train <- getchz(Y_train, explp_train)
  
  if (nrow(df_train) < 1L) {
    stop("coxph_frailty: failed to construct baseline cumulative hazard.", call. = FALSE)
  }
  
  Y_new <- as_counting_surv(mf_new[[1]])
  stop_time_new <- Y_new[, 2]
  status_new <- as.integer(Y_new[, 3])
  n_new <- nrow(Y_new)
  
  mm_nc_new <- ncol(mm_new)
  if (mm_nc_new != (ncol(fix_var_train) + 1L)) {
    stop("coxph_frailty: newdata design matrix shape mismatch.", call. = FALSE)
  }
  
  fix_var_new <- mm_new[, -mm_nc_new, drop = FALSE]
  if (!is.null(colnames(fix_var_train)) && !is.null(colnames(fix_var_new))) {
    fix_var_new <- fix_var_new[, colnames(fix_var_train), drop = FALSE]
  }
  
  lp_fixed_new <- as.vector(fix_var_new %*% beta_use)
  frail_new <- fit$frail[as.integer(group_new)]
  frail_new[is.na(frail_new)] <- 0
  
  eta_new <- exp(lp_fixed_new + frail_new)
  
  idx_after <- findInterval(stop_time_new, df_train$time)
  H_after <- ifelse(idx_after > 0L, df_train$cumhaz[idx_after], 0)
  
  is_event <- (status_new == 1L)
  exact_event <- is_event & idx_after > 0L & df_train$time[idx_after] == stop_time_new
  if (any(is_event & !exact_event)) {
    stop("coxph_frailty: some event times were not found in the training baseline hazard grid.", call. = FALSE)
  }
  
  idx_before <- idx_after
  idx_before[exact_event] <- idx_before[exact_event] - 1L
  H_before <- ifelse(idx_before > 0L, df_train$cumhaz[idx_before], 0)
  
  logS_after <- -eta_new * H_after
  logS_before <- -eta_new * H_before
  
  log_surv <- logS_after
  log_surv[!is_event] <- -Inf
  
  log_like <- logS_after
  if (any(is_event)) {
    log_like[is_event] <- log_diff_exp(logS_before[is_event], logS_after[is_event])
  }
  
  is_discrete <- as.integer(!is_event)
  
  list(
    log_surv = matrix(log_surv, nrow = 1L, ncol = n_new),
    log_like = matrix(log_like, nrow = 1L, ncol = n_new),
    is_discrete = matrix(is_discrete, nrow = 1L, ncol = n_new)
  )
}

log_pointpred_survival_coxph_penal <- function(fit, data, traindata, ...) {
  log_pointpred_survival_coxph_penal_frailty (
    fit = fit,
    data = data,
    traindata = traindata,
    ...
  )
}
