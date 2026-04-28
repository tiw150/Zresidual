#' Extract aligned model metadata for Z-residual diagnostics
#'
#' @description
#' Extracts model-aligned response information, response type, covariates, and
#' linear predictors for use in Z-residual diagnostics and plotting functions.
#' The exact components returned depend on the fitted model class and the value
#' of \code{detail}.
#'
#' @param fit A fitted model object, such as \code{brmsfit},
#'   \code{survival::coxph}, or \code{survival::survreg}.
#' @param data Data used to extract aligned metadata. Must be provided.
#' @param detail Character string; either \code{"basic"} or \code{"full"}.
#' @param type Optional component selector. For example, in \code{brms} hurdle
#'   models, use \code{type = "zero"} to isolate the binary hurdle component,
#'   \code{type = "count"} for the positive-count component, and
#'   \code{type = "hurdle"} for the full hurdle model.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return
#' A named list. Depending on the model and \code{detail}, this typically
#' includes aligned response information, response-type metadata, covariates,
#' linear predictors, and optional model-frame details.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 20
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.2 * x))
#'   t_cens  <- rexp(n, rate = 0.5)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   info <- Zcov(fit, data = dat)
#'   names(info)
#' }
#'
#' @export
Zcov <- function(fit,
                 data,
                 detail = c("basic", "full"),
                 type = NULL,
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
  
  na_cols <- names(mf)[vapply(mf, function(x) any(is.na(x)), logical(1))]
  
  show_k <- 20L
  show_rows <- if (drop_n <= show_k) {
    paste(drop_id, collapse = ",")
  } else {
    paste0(paste(drop_id[seq_len(show_k)], collapse = ","), ", ... (+", drop_n - show_k, " more)")
  }
  
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

# Normalize config in a future-proof way:
# - keep unknown fields
# - only interpret a small stable subset now
.Zcov_config_normalize <- function(config, family = NULL, is_hurdle = FALSE) {
  if (is.null(config)) config <- list()
  if (!is.list(config)) stop("Zcov: config must be a list.", call. = FALSE)
  
  cfg_raw <- config
  
  # stable fields (future-proof contract)
  type_req <- NULL
  if (!is.null(config$type)) {
    type_req <- config$type
    if (!is.character(type_req) || length(type_req) != 1L) {
      stop("Zcov: config$type must be a single string.", call. = FALSE)
    }
    type_req <- tolower(type_req)
    if (type_req %in% c("overall", "whole", "all", "model")) type_req <- "hurdle"
    if (type_req %in% c("logistic", "logit")) type_req <- "zero"
  } else {
    # default behavior: if hurdle family -> whole hurdle metadata; else no component selection
    type_req <- if (is_hurdle) "hurdle" else NULL
  }
  
  # lp_scale: controls whether linear_pred is mean/prob vs eta
  lp_scale <- NULL
  if (!is.null(config$lp_scale)) {
    lp_scale <- tolower(as.character(config$lp_scale)[1])
    if (!(lp_scale %in% c("mean", "prob", "eta"))) {
      stop("Zcov: config$lp_scale must be one of 'mean'/'prob'/'eta'.", call. = FALSE)
    }
  }
  
  
  # brms sub-config (future use)
  brms_cfg <- list()
  if (!is.null(config$brms)) {
    if (!is.list(config$brms)) stop("Zcov: config$brms must be a list.", call. = FALSE)
    brms_cfg <- config$brms
  }
  
  # Defaults based on type (only for hurdle brms)
  if (is_hurdle) {
    if (is.null(type_req) || !(type_req %in% c("hurdle", "count", "zero"))) {
      stop("Zcov: for hurdle_* families, config$type must be 'hurdle'/'count'/'zero'.", call. = FALSE)
    }
    
    if (is.null(lp_scale)) {
      lp_scale <- if (type_req == "zero") "prob" else "mean"
    }
    
    if (is.null(brms_cfg$dpar)) {
      brms_cfg$dpar <- if (type_req == "zero") "hu" else if (type_req == "count") "mu" else NULL
    }
    if (is.null(brms_cfg$re_formula)) brms_cfg$re_formula <- NULL
  }
  
  list(
    type     = type_req,
    lp_scale = lp_scale,
    brms     = brms_cfg,
    family   = family,
    raw      = cfg_raw
  )
}

.Zcov_resp_name_from_brms <- function(fit, mf_all) {
  if (is.null(mf_all) || !is.data.frame(mf_all) || ncol(mf_all) < 1L) {
    stop("Zcov.brmsfit: unable to determine response from model frame.", call. = FALSE)
  }
  
  names(mf_all)[1L]
}

# ---- brmsfit method ----
#' @describeIn Zcov Method for \code{brmsfit} objects. 
#'   Currently, the top-level \code{type} argument is 
#'   implemented \strong{only} for hurdle or zero-inflated models. 
#'   For these models, setting \code{type = "zero"} isolates the 
#'   binary zero-inflation/hurdle process, \code{type = "count"} isolates the 
#'   truncated count process, and \code{type = "hurdle"} evaluates the entire model. 
#'   For all other standard families (e.g., standard logistic, Poisson, Gaussian), 
#'   the \code{type} argument is unnecessary and will be ignored with a warning. 
#'   For these standard models, the function simply returns the standard linear 
#'   predictors, covariates, and response.
#' @export
Zcov.brmsfit <- function(fit,
                         data,
                         detail = c("basic", "full"),
                         type = NULL,
                         ...) {
  
  detail <- match.arg(detail)
  
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Zcov.brmsfit requires package 'brms'.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zcov.brmsfit: `data` must be provided.", call. = FALSE)
  }
  
  type_req <- NULL
  if (!is.null(type)) {
    type_req <- tolower(as.character(type)[1])
    if (type_req %in% c("overall", "whole", "all", "model")) type_req <- "hurdle"
    if (type_req %in% c("logistic", "logit")) type_req <- "zero"
  }
  
  data_in <- data
  
  fam <- fit$family$family
  is_hurdle <- grepl("^hurdle(_|$)", fam)
  is_trunc  <- grepl("^truncated(_|$)", fam)
  
  # -------- HURDLE CHECK & WARNING --------
  # If not a hurdle model, warn the user if they tried to pass a 'type' argument
  if (!is_hurdle) {
    if (!is.null(type_req)) {
      warning(
        sprintf("Zcov.brmsfit: 'type' is implemented for hurdle/zero-inflated models only. Ignoring type='%s' for family='%s' and returning standard linear predictors.", 
                type_req, fam), 
        call. = FALSE
      )
    }
    type_req <- NULL
  } else if (is.null(type_req) || !type_req %in% c("hurdle", "count", "zero")) {
    type_req <- "hurdle"
  }
  
  # -------- build model frame and drop NA rows (align obs_id) --------
  mf_all <- stats::model.frame(fit$formula, data = data_in, na.action = stats::na.pass)
  
  # response name
  response_name <- .Zcov_resp_name_from_brms(fit, mf_all)
  
  keep <- stats::complete.cases(mf_all)
  obs_id <- which(keep)
  mf <- mf_all[keep, , drop = FALSE]
  data_used <- data_in[keep, , drop = FALSE]
  
  # response
  y_raw <- mf[[response_name]]
  y_num <- suppressWarnings(as.numeric(y_raw))
  
  # covariates (drop response column)
  cov_all  <- data_in[, setdiff(names(data_in), response_name), drop = FALSE]
  cov_used <- cov_all[obs_id, , drop = FALSE]
  
  # hurdle bookkeeping (indices are within kept rows)
  zero_id <- which(!is.na(y_num) & (y_num == 0))
  zero_obs_id <- obs_id[zero_id]
  
  # y_type: 0=zero, 1=count (for hurdle); else plain=1
  if (is_hurdle) {
    y_type_kind   <- "hurdle"
    y_type_levels <- c(zero = 0L, count = 1L)
    y_type <- ifelse(!is.na(y_num) & (y_num == 0), 0L, 1L)
  } else {
    y_type_kind   <- "plain"
    y_type_levels <- c(obs = 1L)
    y_type <- rep.int(1L, length(obs_id))
  }
  
  # -------- LP: always use conditional predict like old code --------
  lp_cond <- tryCatch({
    pr <- stats::predict(fit, newdata = data_used, type = "conditional")
    if (is.matrix(pr) || is.data.frame(pr)) as.numeric(pr[, 1]) else as.numeric(pr)
  }, error = function(e1) {
    # fallback: fitted() generic (brmsfit has method)
    out <- tryCatch(stats::fitted(fit, newdata = data_used, type = "conditional", re_formula = NULL),
                    error = function(e2) NULL)
    if (is.null(out)) NULL
    else if (is.matrix(out) || is.data.frame(out)) as.numeric(out[, 1]) else as.numeric(out)
  })
  
  if (is.null(lp_cond)) {
    stop("Zcov.brmsfit: failed to compute linear_pred via predict/fitted(type='conditional').", call. = FALSE)
  }
  
  extra <- list(
    family = fam,
    zero_id = zero_id,
    zero_obs_id = zero_obs_id,
    type_raw = type
  )
  
  # -------- truncated_* handling (keep y>0) --------
  if (is_trunc) {
    keep2 <- which(!is.na(y_num) & (y_num > 0))
    obs_id <- obs_id[keep2]
    cov_used <- cov_used[keep2, , drop = FALSE]
    y_raw <- y_raw[keep2]
    y_num <- y_num[keep2]
    lp_cond <- lp_cond[keep2]
    
    y_type_kind   <- "trunc"
    y_type_levels <- c(count = 1L)
    y_type <- rep.int(1L, length(obs_id))
    
    return(list(
      type          = "count",
      family        = fam,
      response_name = response_name,
      response      = y_raw,
      covariates    = cov_used,
      linear_pred   = lp_cond,
      obs_id        = obs_id,
      y_type        = y_type,
      y_type_kind   = y_type_kind,
      y_type_levels = y_type_levels,
      extra         = extra
    ))
  }
  
  # -------- hurdle component selection --------
  if (is_hurdle && identical(type_req, "count")) {
    keep_cnt <- which(!is.na(y_num) & (y_num > 0))
    obs_id_cnt <- obs_id[keep_cnt]
    cov_cnt    <- cov_used[keep_cnt, , drop = FALSE]
    y_cnt      <- y_raw[keep_cnt]
    lp_cnt     <- lp_cond[keep_cnt]
    
    return(list(
      type          = "count",
      family        = fam,
      response_name = response_name,
      response      = y_cnt,
      covariates    = cov_cnt,
      linear_pred   = lp_cnt,          
      obs_id        = obs_id_cnt,
      y_type        = rep.int(1L, length(keep_cnt)),
      y_type_kind   = "trunc",
      y_type_levels = c(count = 1L),
      extra         = extra
    ))
  }
  
  if (is_hurdle && identical(type_req, "zero")) {
    # logistic target: O_i = 1[y==0]
    o_i <- as.integer(!is.na(y_num) & (y_num == 0))
    
    # keep hu_hat optionally for users, but DO NOT use as linear_pred (to match old behavior)
    if (detail == "full") {
      extra$hu_hat <- tryCatch({
        mat <- brms::posterior_epred(fit, newdata = data_used, dpar = "hu", re_formula = NULL)
        as.numeric(colMeans(mat))
      }, error = function(e) NULL)
    }
    
    return(list(
      type          = "zero",
      family        = fam,
      response_name = response_name,
      response      = o_i,
      covariates    = cov_used,
      linear_pred   = lp_cond,         
      obs_id        = obs_id,
      y_type        = y_type,          # 0=zero, 1=count
      y_type_kind   = "hurdle",
      y_type_levels = c(zero = 0L, count = 1L),
      extra         = extra
    ))
  }
  
  # default: hurdle whole model (or non-hurdle standard model)
  list(
    type          = if (is_hurdle) "hurdle" else NULL,
    family        = fam,
    response_name = response_name,
    response      = y_raw,
    covariates    = cov_used,
    linear_pred   = lp_cond,          
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
#' @describeIn Zcov Method for \code{coxph} objects from the \pkg{survival} package.
#'   Extracts the \code{Surv} response, linear predictors, and censoring metadata.
#' @export
Zcov.coxph <- function(fit,
                       data,
                       detail = c("basic", "full"),
                       type = NULL,
                       ...) {
  detail <- match.arg(detail)
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.coxph requires package 'survival'.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zcov.coxph: `data` must be provided.", call. = FALSE)
  }
  
  extra_cfg <- type
  
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
  
  if (!inherits(Y_new, "Surv")) stop("Zcov.coxph: response is not a Surv object.", call. = FALSE)
  
  y_mat <- as.matrix(Y_new)
  if (ncol(y_mat) == 2) {
    time   <- y_mat[, 1]
    status <- y_mat[, 2]
  } else if (ncol(y_mat) == 3) {
    time   <- y_mat[, 2]
    status <- y_mat[, 3]
  } else {
    stop("Zcov.coxph: unsupported Surv format.", call. = FALSE)
  }
  
  y_type <- ifelse(status == 1, 1L, 0L)
  
  extra <- list(
    time = time,
    status = status,
    event_id = which(status == 1),
    censor_id = which(status == 0),
    type_raw = extra_cfg
  )
  response_name <- .Zcov_surv_response_name(fit, mf_new)
  extra <- c(extra, .Zcov_surv_arg_names(fit))
  
  list(
    type          = NULL,
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

# ---- survreg method ----
#' @describeIn Zcov Method for \code{survreg} objects from the \pkg{survival} package.
#'   Extracts the \code{Surv} response, linear predictors, distribution scale, and
#'   censoring metadata.
#' @export
Zcov.survreg <- function(fit,
                         data,
                         detail = c("basic", "full"),
                         type = NULL,
                         ...) {
  detail <- match.arg(detail)
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.survreg requires package 'survival'.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("Zcov.survreg: `data` must be provided.", call. = FALSE)
  }
  
  extra_cfg <- type
  
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
  
  if (!inherits(y, "Surv")) stop("Zcov.survreg: response is not a Surv object.", call. = FALSE)
  
  y_mat <- as.matrix(y)
  if (ncol(y_mat) < 2) stop("Zcov.survreg: Surv object must have (time, status).", call. = FALSE)
  
  time   <- y_mat[, 1]
  status <- y_mat[, 2]
  
  y_type <- ifelse(status == 1, 1L, 0L)
  
  extra <- list(
    time = time,
    status = status,
    event_id = which(status == 1),
    censor_id = which(status == 0),
    dist = fit$dist,
    scale = fit$scale,
    type_raw = extra_cfg
  )
  
  response_name <- .Zcov_surv_response_name(fit, mf)
  extra <- c(extra, .Zcov_surv_arg_names(fit))
  
  list(
    type          = NULL,
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


# ---- glm method ----
#' @describeIn Zcov Method for \code{glm} objects.
#'   Standard generalized linear models are single-process models. 
#'   Therefore, the \code{type} argument is unnecessary and will be 
#'   ignored with a warning if provided. The function returns the 
#'   standard linear predictors on the link scale, covariates, and response.
#' @export
Zcov.glm <- function(fit,
                     data,
                     detail = c("basic", "full"),
                     type = NULL,
                     ...) {
  
  detail <- match.arg(detail)
  
  if (missing(data) || is.null(data)) {
    stop("Zcov.glm: `data` must be provided.", call. = FALSE)
  }
  
  fam <- fit$family$family
  
  # -------- TYPE CHECK & WARNING --------
  # glm models do not natively support hurdle/zero-inflation sub-components
  if (!is.null(type)) {
    warning(
      sprintf("Zcov.glm: 'type' argument is implemented for hurdle/zero-inflated models. Ignoring type='%s' for standard glm family='%s'.", 
              type, fam), 
      call. = FALSE
    )
  }
  
  # -------- build model frame and drop NA rows (align obs_id) --------
  mf_all <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  
  # response name (typically the first column in the model frame)
  response_name <- names(mf_all)[1L]
  
  keep <- stats::complete.cases(mf_all)
  obs_id <- which(keep)
  mf <- mf_all[keep, , drop = FALSE]
  data_used <- data[keep, , drop = FALSE]
  
  # response
  y_raw <- stats::model.response(mf)
  
  # covariates (drop response column from the original data)
  cov_all  <- data[, setdiff(names(data), response_name), drop = FALSE]
  cov_used <- cov_all[obs_id, , drop = FALSE]
  
  # y_type: glm models are plain single-component models
  y_type_kind   <- "plain"
  y_type_levels <- c(obs = 1L)
  y_type <- rep.int(1L, length(obs_id))
  
  # -------- LP: compute linear predictor on the link scale --------
  lp_cond <- tryCatch({
    pr <- stats::predict(fit, newdata = data_used, type = "link")
    if (is.matrix(pr) || is.data.frame(pr)) as.numeric(pr[, 1]) else as.numeric(pr)
  }, error = function(e1) {
    stop("Zcov.glm: failed to compute linear_pred via predict(type='link').", call. = FALSE)
  })
  
  extra <- list(
    family = fam,
    type_raw = type
  )
  
  # return aligned metadata
  list(
    type          = NULL,
    family        = fam,
    response_name = response_name,
    response      = y_raw,
    covariates    = cov_used,
    linear_pred   = lp_cond,          
    obs_id        = obs_id,
    y_type        = y_type,
    y_type_kind   = y_type_kind,
    y_type_levels = y_type_levels,
    extra         = extra
  )
}

#' @method Zcov default
#' @export
#' @noRd
Zcov.default <- function(fit,
                         data = NULL,
                         detail = c("basic", "full"),
                         type = NULL,
                         ...) {
  stop(sprintf(
    "Zcov: unsupported model class: %s",
    paste(class(fit), collapse = "/")
  ), call. = FALSE)
}


