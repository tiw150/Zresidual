
#' Cross-validated Z-residuals
#'
#' Generic function for computing cross-validated Z-residuals.
#'
#' @param object A fitted model object.
#' @param nfolds Integer number of folds. If \code{NULL} and \code{foldlist} is
#'   also \code{NULL}, a default value is chosen.
#' @param foldlist Optional custom list of test-set indices, one element per
#'   fold.
#' @param data Data used for cross-validation refitting. Must be supplied.
#' @param nrep Integer number of Z-residual replicates.
#' @param log_pointpred Optional function or function name passed to
#'   \code{\link{Zresidual}}.
#' @param mcmc_summarize Posterior summarization method for Bayesian fits.
#' @param type Optional component selector used by \code{Zresidual()} and
#'   \code{Zcov()}.
#' @param randomized Logical; whether to generate randomized Z-residuals.
#' @param seed Optional integer seed.
#' @param ... Further arguments passed to \code{\link{Zresidual}}.
#'
#' @return
#' A numeric matrix of class \code{"cvzresid"} and \code{"zresid"}, with one row
#' per observation and one column per replicate.
#'
#' @export
CV.Zresidual <- function(object,
                         nfolds = NULL,
                         foldlist = NULL,
                         data = NULL,
                         nrep = 1,
                         log_pointpred = NULL,
                         mcmc_summarize = c("post", "iscv"),
                         type = NULL,
                         randomized = TRUE,
                         seed = NULL,
                         ...) {
  UseMethod("CV.Zresidual")
}

.cv_surv_status <- function(y) {
  ymat <- as.matrix(y)
  if (ncol(ymat) < 2L) {
    stop("CV.Zresidual: survival response must contain a status column.", call. = FALSE)
  }
  as.integer(ymat[, ncol(ymat)])
}

.cv_surv_time <- function(y) {
  ymat <- as.matrix(y)
  if (ncol(ymat) == 2L) return(ymat[, 1L])
  if (ncol(ymat) >= 3L) return(ymat[, 2L])
  rep(NA_real_, nrow(ymat))
}

.cv_surv_mf_to_data <- function(fit, mf) {
  y <- mf[[1L]]
  ymat <- as.matrix(y)
  nm <- tryCatch(.Zcov_surv_arg_names(fit), error = function(e) list())
  
  if (ncol(ymat) == 2L) {
    lhs_names <- c(
      if (!is.null(nm$time_name)) nm$time_name else "time",
      if (!is.null(nm$status_name)) nm$status_name else "status"
    )
  } else if (ncol(ymat) == 3L) {
    lhs_names <- c(
      if (!is.null(nm$start_name)) nm$start_name else "start",
      if (!is.null(nm$stop_name)) nm$stop_name else "stop",
      if (!is.null(nm$status_name)) nm$status_name else "status"
    )
  } else {
    stop("CV.Zresidual: unsupported Surv response format.", call. = FALSE)
  }
  
  lhs <- as.data.frame(ymat, check.names = FALSE)
  names(lhs) <- lhs_names
  rhs <- mf[, -1L, drop = FALSE]
  out <- data.frame(lhs, rhs, check.names = FALSE)
  rownames(out) <- rownames(mf)
  out
}

.cv_prepare_surv_data <- function(object, data = NULL, where = "CV.Zresidual") {
  if (!is.null(data)) {
    mf_all <- stats::model.frame(object$terms, data, na.action = stats::na.pass)
    keep <- stats::complete.cases(mf_all)
    .Zcov_warn_drop_na(where, nrow(mf_all), keep, mf_all)
    return(data[keep, , drop = FALSE])
  }
  
  mf <- if (inherits(object, "coxph")) {
    model.frame.coxph(object)
  } else if (inherits(object, "survreg")) {
    model.frame.survreg(object)
  } else {
    stop("CV.Zresidual: unsupported survival model class.", call. = FALSE)
  }
  
  .cv_surv_mf_to_data(object, mf)
}

.cv_make_or_check_foldlist <- function(y, covariates, censor, nfolds, foldlist) {
  n <- NROW(covariates)
  
  if (!is.null(foldlist)) {
    foldlist <- lapply(foldlist, function(idx) sort(unique(as.integer(idx))))
    flat <- unlist(foldlist, use.names = FALSE)
    
    if (length(flat) != n || anyNA(flat) || any(flat < 1L) || any(flat > n)) {
      stop("CV.Zresidual: `foldlist` must contain a valid partition of 1:n.", call. = FALSE)
    }
    if (!identical(sort(flat), seq_len(n))) {
      stop("CV.Zresidual: `foldlist` must partition 1:n without overlap or omission.", call. = FALSE)
    }
    return(foldlist)
  }
  
  if (is.null(nfolds)) {
    nfolds <- min(10L, n)
  }
  nfolds <- as.integer(nfolds)
  
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds < 2L) {
    stop("CV.Zresidual: `nfolds` must be an integer >= 2.", call. = FALSE)
  }
  if (nfolds > n) {
    stop("CV.Zresidual: `nfolds` cannot exceed the number of observations.", call. = FALSE)
  }
  
  make_fold(
    fix_var = covariates,
    y = y,
    k = nfolds,
    censor = censor
  )
}

.cv_refit_object <- function(object, train_data) {
  tryCatch(
    suppressWarnings(stats::update(object, data = train_data)),
    error = function(e) NULL
  )
}

.cv_predict_lp_safe <- function(fit, data, type = NULL) {
  zcov <- tryCatch(
    Zcov(fit, data = data, type = type),
    error = function(e) NULL
  )
  if (is.null(zcov) || is.null(zcov$linear_pred)) {
    return(rep(NA_real_, nrow(data)))
  }
  as.numeric(zcov$linear_pred)
}

.cv_default_log_pointpred <- function(object, log_pointpred = NULL, type = NULL) {
  if (!is.null(log_pointpred)) return(log_pointpred)
  
  frailty_fun <- get0(
    "log_pointpred_survival_coxph_frailty",
    mode = "function",
    inherits = TRUE
  )
  
  if (inherits(object, "coxph")) {
    frailty_terms <- attr(object$terms, "specials")$frailty
    if (!is.null(frailty_terms) && length(frailty_terms) > 0L) {
      return(frailty_fun)
    }
  }
  
  if (!is.null(type) && identical(tolower(type), "frailty")) {
    return(frailty_fun)
  }
  
  NULL
}

.cv_run_survival_cv <- function(object,
                                nfolds = NULL,
                                foldlist = NULL,
                                data = NULL,
                                nrep = 1,
                                log_pointpred = NULL,
                                mcmc_summarize = c("post", "iscv"),
                                type = NULL,
                                randomized = TRUE,
                                seed = NULL,
                                ...) {
  if (!is.null(seed)) set.seed(seed)
  
  mcmc_summarize <- match.arg(tolower(as.character(mcmc_summarize)), c("post", "iscv"))
  nrep <- as.integer(nrep)
  if (length(nrep) != 1L || is.na(nrep) || nrep < 1L) {
    stop("CV.Zresidual: `nrep` must be an integer >= 1.", call. = FALSE)
  }
  nrep_eff <- if (isTRUE(randomized)) nrep else 1L
  
  data_used <- .cv_prepare_surv_data(object, data = data)
  mf_full <- stats::model.frame(object$terms, data_used)
  
  y_full <- mf_full[[1L]]
  cov_full <- mf_full[, -1L, drop = FALSE]
  censor_full <- .cv_surv_status(y_full)
  time_full <- .cv_surv_time(y_full)
  
  foldlist <- .cv_make_or_check_foldlist(
    y = y_full,
    covariates = cov_full,
    censor = censor_full,
    nfolds = nfolds,
    foldlist = foldlist
  )
  
  lp_fun <- .cv_default_log_pointpred(object, log_pointpred = log_pointpred, type = type)
  
  n <- nrow(data_used)
  z_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  rsp_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  lp_full <- rep(NA_real_, n)
  surv_prob_full <- rep(NA_real_, n)
  
  for (fid in seq_along(foldlist)) {
    idx <- foldlist[[fid]]
    
    train_data <- data_used[-idx, , drop = FALSE]
    test_data  <- data_used[ idx, , drop = FALSE]
    
    fit_train <- .cv_refit_object(object, train_data)
    if (is.null(fit_train)) next
    
    extra_args <- list(...)
    frailty_fun <- get0(
      "log_pointpred_survival_coxph_frailty",
      mode = "function",
      inherits = TRUE
    )
    
    if (identical(lp_fun, frailty_fun) || identical(tolower(type), "frailty")) {
      extra_args$traindata <- train_data
    }
    
    z_call <- c(list(
      fit = fit_train,
      data = test_data,
      log_pointpred = lp_fun,
      mcmc_summarize = mcmc_summarize,
      type = type,
      randomized = randomized,
      nrep = nrep_eff
    ), extra_args)
    z_fold <- do.call(Zresidual, z_call)
    
    pre_call <- c(list(
      fit = fit_train,
      data = test_data,
      log_pointpred = lp_fun,
      mcmc_summarize = mcmc_summarize,
      type = type
    ), extra_args)
    pre_fold <- do.call(log_summary_pred, pre_call)
    
    z_full[idx, ] <- as.matrix(z_fold)
    
    rsp_fold <- attr(z_fold, "rsp")
    if (!is.null(rsp_fold)) {
      rsp_full[idx, ] <- as.matrix(rsp_fold)
    }
    
    lp_full[idx] <- .cv_predict_lp_safe(fit_train, test_data, type = type)
    surv_prob_full[idx] <- exp(pre_fold$log_surv_hat)
  }
  
  colnames(z_full) <- paste0("rep", seq_len(nrep_eff))
  class(z_full) <- c("cvzresid", "zresid", "matrix", "array")
  
  attr(z_full, "rsp") <- rsp_full
  attr(z_full, "type") <- "survival"
  attr(z_full, "Survival.Prob") <- surv_prob_full
  attr(z_full, "linear.pred") <- lp_full
  attr(z_full, "linear_pred") <- lp_full
  attr(z_full, "censored.status") <- censor_full
  attr(z_full, "covariates") <- cov_full
  attr(z_full, "object.model.frame") <- mf_full
  
  attr(z_full, "zcov") <- list(
    type = NULL,
    family = class(object)[1L],
    response_name = names(mf_full)[1L],
    response = y_full,
    covariates = cov_full,
    linear_pred = lp_full,
    obs_id = seq_len(n),
    y_type = ifelse(censor_full == 1L, 1L, 0L),
    y_type_kind = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra = list(
      time = time_full,
      status = censor_full,
      Survival.Prob = surv_prob_full,
      cv = TRUE
    )
  )
  
  z_full
}
