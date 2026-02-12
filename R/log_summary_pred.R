log_summary_pred <- function(fit,
                             data = NULL,
                             config = list(),
                             log_pointpred_fn = NULL,
                             ...) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  lp_norm_key <- function(x) {
    x <- tolower(as.character(x %||% ""))
    x <- gsub("[^a-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }
  
  as_mat <- function(z) {
    if (is.null(z)) return(NULL)
    if (is.vector(z)) return(matrix(as.vector(z), nrow = 1L))
    if (is.matrix(z)) return(z)
    stop("log_summary_pred: lsf_hat/lpmf_hat must be vector or matrix.", call. = FALSE)
  }
  
  # ---- embedded posterior summary: E[S], E[p] ----
  post_hat <- function(log_surv_m, log_pmf_m = NULL) {
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      stop("log_summary_pred: Bayesian summaries require package 'matrixStats'.", call. = FALSE)
    }
    Tdraw <- nrow(log_surv_m)
    log_surv_hat <- matrixStats::colLogSumExps(log_surv_m) - log(Tdraw)
    
    if (is.null(log_pmf_m) || all(is.na(log_pmf_m))) {
      log_pmf_hat <- rep(NA_real_, ncol(log_surv_m))
    } else {
      if (nrow(log_pmf_m) != Tdraw || ncol(log_pmf_m) != ncol(log_surv_m)) {
        stop("post: log_surv/log_pmf dimension mismatch.", call. = FALSE)
      }
      log_pmf_hat <- matrixStats::colLogSumExps(log_pmf_m) - log(Tdraw)
    }
    list(log_surv_hat = log_surv_hat, log_pmf_hat = log_pmf_hat)
  }
  
  # ---- embedded ISCV summary: weights w=1/p(y|theta) ----
  iscv_hat <- function(log_surv_m, log_pmf_m) {
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      stop("log_summary_pred: Bayesian summaries require package 'matrixStats'.", call. = FALSE)
    }
    if (is.null(log_pmf_m) || all(is.na(log_pmf_m))) {
      stop("iscv: requires log_pmf (importance weights 1/p).", call. = FALSE)
    }
    if (anyNA(log_pmf_m)) {
      stop("iscv: log_pmf contains NA; use -Inf for zero probabilities.", call. = FALSE)
    }
    Tdraw <- nrow(log_pmf_m)
    
    # log Σ S/p  and log Σ 1/p
    log_sum_S_over_p <- matrixStats::colLogSumExps(log_surv_m - log_pmf_m)
    log_sum_inv_p    <- matrixStats::colLogSumExps(-log_pmf_m)
    
    log_surv_hat <- log_sum_S_over_p - log_sum_inv_p
    log_pmf_hat  <- log(Tdraw)       - log_sum_inv_p
    
    list(log_surv_hat = log_surv_hat, log_pmf_hat = log_pmf_hat)
  }
  
  # ---- main ----
  config <- config %||% list()
  method <- lp_norm_key(config$method %||% "post")   # default post
  
  lp_fn <- log_pointpred_fn %||% log_pointpred
  pp <- lp_fn(fit, data = data, config = config, ...)
  
  if (!is.list(pp) || is.null(pp$lsf_hat)) {
    stop("log_summary_pred: log_pointpred must return a list with lsf_hat.", call. = FALSE)
  }
  
  log_surv0 <- as_mat(pp$lsf_hat)
  log_pmf0  <- as_mat(pp$lpmf_hat)
  y_type    <- pp$y_type %||% NULL
  
  pmf_available <- !is.null(log_pmf0) && !all(is.na(log_pmf0))
  is_bayes <- nrow(log_surv0) > 1L || (pmf_available && nrow(log_pmf0) > 1L)
  
  # Frequentist: 1×n -> vector
  if (!is_bayes) {
    log_surv_hat <- as.numeric(log_surv0[1, ])
    log_pmf_hat  <- if (pmf_available) as.numeric(log_pmf0[1, ]) else rep(NA_real_, length(log_surv_hat))
    return(list(log_surv_hat = log_surv_hat, log_pmf_hat = log_pmf_hat, y_type = y_type))
  }
  
  # Bayesian: draws×n -> vector
  if (method == "post") {
    hat <- post_hat(log_surv0, log_pmf0)
  } else if (method == "iscv") {
    hat <- iscv_hat(log_surv0, log_pmf0)
  } else {
    stop("log_summary_pred: config$method must be 'post' or 'iscv'.", call. = FALSE)
  }
  
  list(log_surv_hat = as.numeric(hat$log_surv_hat),
       log_pmf_hat  = as.numeric(hat$log_pmf_hat),
       y_type       = y_type)
}
