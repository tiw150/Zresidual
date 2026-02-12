log_pointpred_survival_coxph_frailty <- function(fit,
                                                 traindata = NULL,
                                                 newdata = NULL,
                                                 data = NULL,
                                                 ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("log_pointpred_survival_coxph_frailty requires package 'survival'.", call. = FALSE)
  }
  
  # Zresidual passes data=..., treat it as newdata
  if (!is.null(data)) newdata <- data
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  is_surv_mf <- function(x) {
    is.data.frame(x) && ncol(x) >= 1L && inherits(x[[1]], "Surv") && !is.null(attr(x, "terms"))
  }
  
  make_mf <- function(dat) {
    if (is_surv_mf(dat)) return(dat)
    stats::model.frame(fit$terms, dat)
  }
  
  # Parse frailty(group, ...) -> "group"
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
  
  # --- default data routing ---
  # If user provided newdata but not traindata, use newdata as training proxy (must contain group var)
  if (!is.null(newdata) && is.null(traindata)) {
    traindata <- newdata
  }
  
  # If still missing, try to recover from fit call, otherwise fallback to model.frame(fit)
  if (is.null(traindata)) {
    td <- NULL
    if (!is.null(fit$call$data)) {
      td <- tryCatch(eval(fit$call$data, envir = environment(stats::formula(fit))), error = function(e) NULL)
    }
    traindata <- td %||% stats::model.frame(fit)
  }
  
  if (is.null(newdata)) newdata <- traindata
  
  # --- extract group from raw data (NOT from mf_*) ---
  get_group <- function(dat, name, where) {
    if (!is.data.frame(dat)) {
      stop(sprintf("coxph_frailty: %s must be a data.frame.", where), call. = FALSE)
    }
    if (!name %in% names(dat)) {
      stop(sprintf("coxph_frailty: %s missing group variable '%s'.", where, name), call. = FALSE)
    }
    dat[[name]]
  }
  
  group_train_raw <- get_group(traindata, group_var_name, "training data")
  group_new_raw   <- get_group(newdata,   group_var_name, "newdata")
  
  # Align factor levels to fit$frail if possible (most stable)
  frail_levels <- NULL
  if (!is.null(fit$frail) && length(fit$frail) && !is.null(names(fit$frail)) && all(nzchar(names(fit$frail)))) {
    frail_levels <- names(fit$frail)
  }
  
  if (is.null(frail_levels)) {
    group_train <- factor(group_train_raw)
  } else {
    group_train <- factor(as.character(group_train_raw), levels = frail_levels)
  }
  
  if (!is.factor(group_train)) group_train <- factor(group_train)
  if (anyNA(group_train)) {
    stop("coxph_frailty: training data contains group levels not present in fit$frail.", call. = FALSE)
  }
  
  group_new <- factor(as.character(group_new_raw), levels = levels(group_train))
  
  # write back factorized group so model.frame(fit$terms, ...) can evaluate frailty() consistently
  traindata[[group_var_name]] <- group_train
  newdata[[group_var_name]]   <- group_new
  
  # --- build model frames + matrices ---
  mf_train <- make_mf(traindata)
  mm_train <- model.matrix.coxph(fit, mf_train)
  
  mf_new <- make_mf(newdata)
  mm_new <- model.matrix.coxph(fit, mf_new)
  
  # --- Surv to counting-process form (start, stop, status) ---
  as_counting_surv <- function(Y) {
    if (!inherits(Y, "Surv")) stop("coxph_frailty: response is not a Surv object.", call. = FALSE)
    Ymat <- as.matrix(Y)
    if (ncol(Ymat) != 3) {
      return(as.matrix(survival::Surv(rep(0, nrow(Ymat)), Ymat[, 1], Ymat[, 2])))
    }
    Ymat
  }
  
  # --- baseline cumulative hazard ---
  getchz <- function(Y, explp) {
    death <- (Y[, ncol(Y)] == 1)
    dtime <- Y[, ncol(Y) - 1]
    
    time <- sort(unique(dtime))
    if (length(time) < 1) return(data.frame(time = numeric(0), cumhaz = numeric(0)))
    
    nevent <- as.vector(rowsum(1 * death, dtime))
    nrisk  <- rev(cumsum(rev(rowsum(explp, dtime))))
    
    delta <- if (length(time) >= 2) min(diff(time)) / 2 else 0.5
    etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)
    
    indx <- approx(etime, 1:length(etime), time,
                   method = "constant", rule = 2, f = 1)$y
    
    esum  <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))
    nrisk <- nrisk - c(esum, 0)[indx]
    nrisk <- pmax(nrisk, .Machine$double.xmin)
    
    cumhaz <- cumsum(nevent / nrisk)
    data.frame(time = time, cumhaz = cumhaz)
  }
  
  # --- training explp (fixed + frailty) ---
  gpnumber <- nlevels(group_train)
  
  mm_nc <- ncol(mm_train)
  
  # Most frailty fits have one extra column in X (group code). If not, assume no extra.
  has_group_code_col <- (mm_nc == (length(fit$coefficients) + 1L))
  
  if (!has_group_code_col) {
    stop("coxph_frailty: cannot identify group-code column in model matrix; expected ncol(X)=length(coef)+1.", call. = FALSE)
  }
  
  fix_var_train <- mm_train[, -mm_nc, drop = FALSE]
  
  beta <- fit$coefficients
  if (!is.null(names(beta)) && !is.null(colnames(fix_var_train))) {
    idx <- match(colnames(fix_var_train), names(beta))
    if (anyNA(idx)) stop("coxph_frailty: coefficient names do not match design matrix columns.", call. = FALSE)
    beta_use <- beta[idx]
  } else {
    beta_use <- beta
  }
  
  lp_fixed_train <- as.vector(fix_var_train %*% beta_use)
  
  gidx_train <- as.integer(group_train)
  frail_train <- fit$frail[gidx_train]
  frail_train[is.na(frail_train)] <- 0
  
  explp_train <- exp(lp_fixed_train + frail_train)
  
  Y_train  <- as_counting_surv(mf_train[[1]])
  df_train <- getchz(Y_train, explp_train)
  
  if (nrow(df_train) >= 2) {
    f_step <- stats::stepfun(df_train$time[-1], df_train$cumhaz)
  } else if (nrow(df_train) == 1) {
    f_step <- function(t) rep(df_train$cumhaz[1], length(t))
  } else {
    f_step <- function(t) rep(0, length(t))
  }
  
  # --- prediction ---
  Y_new <- as_counting_surv(mf_new[[1]])
  stop_time_new <- Y_new[, 2]
  status_new    <- Y_new[, 3]
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
  
  gidx_new <- as.integer(group_new)
  frail_new <- fit$frail[gidx_new]
  frail_new[is.na(frail_new)] <- 0
  
  z_hat_new <- exp(frail_new)
  explp_fixed_new <- exp(lp_fixed_new)
  
  H0_new <- f_step(stop_time_new)
  
  SP <- exp(-z_hat_new * explp_fixed_new * H0_new)
  logS <- log(pmax(SP, .Machine$double.xmin))
  
  y_type <- ifelse(status_new == 1, 1L, 0L)
  
  list(
    lpmf_hat = matrix(rep(NA_real_, n_new), nrow = 1L, ncol = n_new),
    lsf_hat  = matrix(logS, nrow = 1L, ncol = n_new),
    y_type   = y_type
  )
}

log_pointpred_survival_coxph_penal <- function(fit, data = NULL, ...) {
  log_pointpred_survival_coxph_frailty(fit, data = data, ...)
}
