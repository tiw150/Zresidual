#' Zcov: extract covariates / lp / response / y_type for diagnostics
#'
#' @param fit model object (brmsfit, coxph, survreg, ...)
#' @param data optional newdata
#' @param mode optional; default "auto"
#' @param detail optional; "basic" or "full"
#' @export
Zcov <- function(fit,
                 data = NULL,
                 mode = c("auto", "hurdle_zero"),
                 detail = c("basic", "full"),
                 ...) {
  UseMethod("Zcov")
}

# ---- helpers ----

.Zcov_pick_data <- function(fit, data) {
  if (!is.null(data)) return(data)
  if (!is.null(fit[["data"]])) return(fit[["data"]])
  NULL
}

.Zcov_warn_drop_na <- function(where, data_n, keep, mf) {
  if (all(keep)) return(invisible(NULL))
  
  drop_id <- which(!keep)
  drop_n  <- length(drop_id)
  
  # columns with missing values (in the model frame)
  na_cols <- names(mf)[vapply(mf, function(x) any(is.na(x)), logical(1))]
  
  # compact row index list
  show_k <- 20L
  show_rows <- if (drop_n <= show_k) {
    paste(drop_id, collapse = ",")
  } else {
    paste0(paste(drop_id[seq_len(show_k)], collapse = ","), ", ... (+", drop_n - show_k, " more)")
  }
  
  # per-row NA columns (first few rows only)
  show_r <- min(8L, drop_n)
  if (show_r > 0L) {
    det_rows <- drop_id[seq_len(show_r)]
    det_txt <- vapply(det_rows, function(i) {
      cols_i <- names(mf)[is.na(mf[i, , drop = FALSE])[1, ]]
      paste0(i, ":{", paste(cols_i, collapse = ","), "}")
    }, character(1))
    row_detail <- paste(det_txt, collapse = "; ")
  } else {
    row_detail <- ""
  }
  
  warning(
    sprintf(
      "%s: dropped %d/%d rows with NA in model variables. Dropped rows: %s. NA details (first %d rows): %s. Columns with any NA: %s",
      where, drop_n, data_n, show_rows, show_r,
      if (nzchar(row_detail)) row_detail else "<none>",
      if (length(na_cols)) paste(na_cols, collapse = ",") else "<unknown>"
    ),
    call. = FALSE
  )
  
  invisible(NULL)
}





# ---- brmsfit method ----
#' @rawNamespace S3method(Zcov, brmsfit)
Zcov.brmsfit <- function(fit,
                         data = NULL,
                         mode = c("auto", "hurdle_zero"),
                         detail = c("basic", "full"),
                         ...) {
  
  mode   <- match.arg(mode)
  detail <- match.arg(detail)
  
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Zcov.brmsfit requires package 'brms'.", call. = FALSE)
  }
  
  data_in <- .Zcov_pick_data(fit, data)
  if (is.null(data_in)) {
    stop("Zcov.brmsfit: no data available (provide data or fit$data).", call. = FALSE)
  }
  
  # Keep NA rows first; then explicitly drop so obs_id aligns to input rows.
  mf_all <- model.frame(fit$formula, data = data_in, na.action = stats::na.pass)
  response_name <- .Zcov_resp_name_from_brms(fit, mf_all)
  
  keep <- stats::complete.cases(mf_all)
  .Zcov_warn_drop_na("Zcov.brmsfit", nrow(mf_all), keep, mf_all)
  
  obs_id <- which(keep)  # indices in data_in kept
  mf <- mf_all[keep, , drop = FALSE]
  data_used <- data_in[keep, , drop = FALSE]
  
  fam <- fit$family$family
  is_hurdle <- grepl("^hurdle_", fam)
  is_trunc  <- grepl("^truncated_", fam)
  
  # response (keep original type for auto mode)
  y_raw <- mf[[response_name]]
  y_num <- suppressWarnings(as.numeric(y_raw))
  
  # covariates: mimic old code style (data columns excluding response), aligned to kept rows
  cov_all <- data_in[, setdiff(names(data_in), response_name), drop = FALSE]
  cov_used <- cov_all[obs_id, , drop = FALSE]
  
  # hurdle bookkeeping
  zero_id <- which(!is.na(y_num) & (y_num == 0))         # indices within kept rows
  zero_obs_id <- obs_id[zero_id]                         # indices in original data_in
  
  y_type_kind   <- "plain"
  y_type_levels <- c(obs = 1L)
  y_type <- rep.int(1L, length(obs_id))
  if (is_hurdle) {
    y_type_kind   <- "hurdle"
    y_type_levels <- c(zero = 0L, count = 1L)
    y_type <- ifelse(!is.na(y_num) & (y_num == 0), 0L, 1L)
  }
  
  # helper: fitted mean on response scale (matches your old "linear.pred" semantics)
  .brms_fitted_conditional <- function(fit, newdata) {
    # try fitted(type="conditional") first, then fallback
    out <- tryCatch(
      brms::fitted(fit, newdata = newdata, type = "conditional", re_formula = NULL),
      error = function(e1) {
        tryCatch(
          brms::fitted(fit, newdata = newdata, type = "conditional"),
          error = function(e2) {
            tryCatch(
              brms::fitted(fit, newdata = newdata, re_formula = NULL),
              error = function(e3) brms::fitted(fit, newdata = newdata)
            )
          }
        )
      }
    )
    as.vector(out[, 1])
  }
  
  # default linear_pred for auto mode: computed on *data_used* (not training data)
  linear_pred <- .brms_fitted_conditional(fit, data_used)
  
  extra <- list(
    family = fam,
    zero_id = zero_id,
    zero_obs_id = zero_obs_id
  )
  
  # truncated_*: keep only y>0 (optional, consistent with your earlier intent)
  if (is_trunc) {
    keep2 <- which(!is.na(y_num) & (y_num > 0))
    obs_id <- obs_id[keep2]
    mf <- mf[keep2, , drop = FALSE]
    data_used <- data_used[keep2, , drop = FALSE]
    cov_used <- cov_used[keep2, , drop = FALSE]
    y_raw <- y_raw[keep2]
    y_num <- y_num[keep2]
    linear_pred <- linear_pred[keep2]
    
    y_type_kind   <- "trunc"
    y_type_levels <- c(count = 1L)
    y_type <- rep.int(1L, length(obs_id))
  }
  
  # mode == hurdle_zero: keep only y==0 rows; response must be 0 (per your requirement)
  if (mode == "hurdle_zero") {
    if (!is_hurdle) {
      stop("Zcov.brmsfit(mode='hurdle_zero') requires a hurdle_* brms family.", call. = FALSE)
    }
    if (is_trunc) {
      stop("Zcov.brmsfit(mode='hurdle_zero') is not meaningful for truncated_* families.", call. = FALSE)
    }
    if (length(zero_id) == 0L) {
      stop("Zcov.brmsfit(mode='hurdle_zero'): no zero observations after NA dropping.", call. = FALSE)
    }
    
    idx <- zero_id
    obs_id_z <- obs_id[idx]
    cov_z <- cov_used[idx, , drop = FALSE]
    data_z <- data_used[idx, , drop = FALSE]
    
    # response: all zeros
    y_z <- rep.int(0, length(idx))
    
    # linear_pred: prefer hu linear predictor mean; fallback to conditional fitted mean on zeros
    lp_hu <- tryCatch(
      {
        mat <- brms::posterior_linpred(
          fit, newdata = data_z, dpar = "hu", transform = FALSE, re_formula = NULL
        )
        colMeans(mat)
      },
      error = function(e) NULL
    )
    if (is.null(lp_hu)) {
      lp_hu <- .brms_fitted_conditional(fit, data_z)
    }
    
    if (detail == "full") {
      mu_hu <- tryCatch(
        {
          mat <- brms::posterior_epred(fit, newdata = data_z, dpar = "hu", re_formula = NULL)
          colMeans(mat)
        },
        error = function(e) NULL
      )
      extra$mu_hu <- mu_hu
    }
    extra$mode <- "hurdle_zero"
    
    return(list(
      family        = fam,
      response_name = response_name,
      response      = y_z,
      covariates    = cov_z,
      linear_pred   = lp_hu,
      obs_id        = obs_id_z,
      y_type        = rep.int(0L, length(idx)),
      y_type_kind   = "hurdle",
      y_type_levels = c(zero = 0L, count = 1L),
      extra         = extra
    ))
  }
  
  # mode == auto
  if (detail == "full") {
    # optional: attach mu_hat (mean response) for extra diagnostics
    mu_hat <- tryCatch(
      {
        mat <- brms::posterior_epred(fit, newdata = data_used, re_formula = NULL)
        colMeans(mat)
      },
      error = function(e) NULL
    )
    extra$mu_hat <- mu_hat
  }
  
  list(
    family        = fam,
    response_name = response_name,
    response      = y_raw,
    covariates    = cov_used,
    linear_pred   = linear_pred,
    obs_id        = obs_id,
    y_type        = y_type,
    y_type_kind   = y_type_kind,
    y_type_levels = y_type_levels,
    extra         = extra
  )
}


####################################################################################

.Zcov_surv_response_name <- function(fit, mf) {
  rn <- NULL
  if (!is.null(mf) && ncol(mf) >= 1) rn <- names(mf)[1]
  if (!is.null(rn) && nzchar(rn)) return(rn)
  
  f <- tryCatch(stats::formula(fit), error = function(e) NULL)
  if (!is.null(f) && length(f) >= 2) return(deparse(f[[2]]))
  
  "Surv(...)"
}

.Zcov_surv_arg_names <- function(fit) {
  f <- tryCatch(stats::formula(fit), error = function(e) NULL)
  if (is.null(f) || length(f) < 2) return(list())
  
  lhs <- f[[2]]
  if (!is.call(lhs) || !identical(lhs[[1]], as.name("Surv"))) return(list())
  
  args <- as.list(lhs)[-1]
  out <- list(surv_expr = deparse(lhs), surv_args = lapply(args, deparse))
  
  # heuristic mapping
  if (length(args) == 2) {
    out$time_name   <- deparse(args[[1]])
    out$status_name <- deparse(args[[2]])
  } else if (length(args) >= 3) {
    out$start_name  <- deparse(args[[1]])
    out$stop_name   <- deparse(args[[2]])
    out$status_name <- deparse(args[[3]])
  }
  out
}




# ---- coxph method ----
#' @rawNamespace S3method(Zcov, coxph)
Zcov.coxph <- function(fit,
                       data = NULL,
                       detail = c("basic", "full"),
                       ...) {
  detail <- match.arg(detail)
  
  dots <- list(...)
  if ("mode" %in% names(dots)) {
    m <- dots$mode
    if (!is.character(m) || length(m) != 1L) {
      stop("Zcov.coxph: 'mode' must be a single string; only brmsfit hurdle models support mode='hurdle_zero'.")
    }
    if (m != "auto") {
      stop("Zcov.coxph: 'mode' is only supported for brmsfit hurdle models (use mode='auto').")
    }
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.coxph requires package 'survival'.")
  }
  
  if (!is.null(data)) {
    mf_all <- model.frame(fit$terms, data, na.action = stats::na.pass)
    keep   <- stats::complete.cases(mf_all)
    .Zcov_warn_drop_na("Zcov.coxph", nrow(mf_all), keep, mf_all)
    
    data_used <- data[keep, , drop = FALSE]
    mf_new <- model.frame(fit$terms, data_used)
    Y_new  <- mf_new[[1]]
    
    lp <- tryCatch(
      as.vector(predict(fit, newdata = data_used, type = "lp")),
      error = function(e) {
        mm_new <- model.matrix.coxph(fit, data_used)
        as.vector(mm_new %*% fit$coefficients)
      }
    )
    
    obs_id <- which(keep)
  } else {
    mf_new <- model.frame.coxph(fit)
    Y_new  <- mf_new[[1]]
    
    lp <- tryCatch(
      as.vector(predict(fit, type = "lp")),
      error = function(e) {
        mm_new <- model.matrix.coxph(fit)
        as.vector(mm_new %*% fit$coefficients)
      }
    )
    
    obs_id <- seq_len(nrow(mf_new))
  }
  
  if (!inherits(Y_new, "Surv")) stop("Zcov.coxph: response is not a Surv object.")
  
  y_mat <- as.matrix(Y_new)
  if (ncol(y_mat) == 2) {
    time   <- y_mat[, 1]
    status <- y_mat[, 2]
  } else if (ncol(y_mat) == 3) {
    time   <- y_mat[, 2]
    status <- y_mat[, 3]
  } else {
    stop("Zcov.coxph: unsupported Surv format.")
  }
  
  y_type <- ifelse(status == 1, 1L, 0L) # 1=event, 0=censored
  
  extra <- list(
    time = time,
    status = status,
    event_id = which(status == 1),
    censor_id = which(status == 0)
  )
  
  response_name <- .Zcov_surv_response_name(fit, mf_new)
  extra <- c(extra, .Zcov_surv_arg_names(fit))
  
  
  list(
    family        = "coxph",
    response_name = response_name,
    response      = Y_new,
    covariates    = mf_new[, -1, drop = FALSE],
    linear_pred   = lp,
    obs_id        = obs_id,
    y_type        = y_type,
    y_type_kind   = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra         = extra
  )
}

# ---- survreg method (mode only meaningful for brms hurdle; accept via ...) ----
#' @rawNamespace S3method(Zcov, survreg)
Zcov.survreg <- function(fit,
                         data = NULL,
                         detail = c("basic", "full"),
                         ...) {
  detail <- match.arg(detail)
  
  dots <- list(...)
  if ("mode" %in% names(dots)) {
    m <- dots$mode
    if (!is.character(m) || length(m) != 1L) {
      stop("Zcov.survreg: 'mode' must be a single string; only brmsfit hurdle models support mode='hurdle_zero'.")
    }
    if (m != "auto") {
      stop("Zcov.survreg: 'mode' is only supported for brmsfit hurdle models (use mode='auto').")
    }
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.survreg requires package 'survival'.")
  }
  
  if (!is.null(data)) {
    mf_all <- model.frame(fit$terms, data, na.action = stats::na.pass)
    keep   <- stats::complete.cases(mf_all)
    .Zcov_warn_drop_na("Zcov.survreg", nrow(mf_all), keep, mf_all)
    
    data_used <- data[keep, , drop = FALSE]
    mf <- model.frame(fit$terms, data_used)
    y  <- mf[[1]]
    
    lp <- tryCatch(
      as.vector(predict(fit, newdata = data_used, type = "lp")),
      error = function(e) {
        mm <- model.matrix(fit$terms, data_used)
        as.vector(mm %*% fit$coefficients)
      }
    )
    
    obs_id <- which(keep)
  } else {
    mf <- model.frame.survreg(fit)
    y  <- mf[[1]]
    
    lp <- tryCatch(
      as.vector(predict(fit, type = "lp")),
      error = function(e) {
        mm <- model.matrix.survreg(fit)
        as.vector(mm %*% fit$coefficients)
      }
    )
    
    obs_id <- seq_len(nrow(mf))
  }
  
  if (!inherits(y, "Surv")) stop("Zcov.survreg: response is not a Surv object.")
  
  y_mat <- as.matrix(y)
  if (ncol(y_mat) < 2) stop("Zcov.survreg: Surv object must have (time, status).")
  
  time   <- y_mat[, 1]
  status <- y_mat[, 2]
  
  y_type <- ifelse(status == 1, 1L, 0L)
  
  extra <- list(
    time = time,
    status = status,
    event_id = which(status == 1),
    censor_id = which(status == 0),
    dist = fit$dist,
    scale = fit$scale
  )
  
  response_name <- .Zcov_surv_response_name(fit, mf)
  extra <- c(extra, .Zcov_surv_arg_names(fit))
  
  
  list(
    family        = "survreg",
    response_name = response_name,
    response      = y,
    covariates    = mf[, -1, drop = FALSE],
    linear_pred   = lp,
    obs_id        = obs_id,
    y_type        = y_type,
    y_type_kind   = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra         = extra
  )
}

#' Default method for unsupported model classes
#' @noRd
Zcov.default <- function(fit, data = NULL, ...) {
  stop(sprintf(
    "Zcov: unsupported model class: %s",
    paste(class(fit), collapse = "/")
  ), call. = FALSE)
}
