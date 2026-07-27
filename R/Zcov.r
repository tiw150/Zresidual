#' Extract aligned covariates and optional point-level predictive quantities
#'
#' @description
#' `Zcov()` extracts model-aligned response information, covariates, response
#' type metadata, and linear predictors for Z-residual diagnostics. When
#' `point_details` is supplied, it also returns point-level quantities such as
#' mean, variance, and lp matrices for predictive checks.
#'
#' @param fit A fitted model object.
#' @param data Data used to extract aligned metadata. Must be provided for most
#'   model classes.
#' @param type Optional component selector. For hurdle models, common values are
#'   `"hurdle"`, `"count"`, and `"zero"`.
#' @param point_details Optional character vector specifying point-level
#'   quantities to return. Supported values are `"mean"`, `"var"`,
#'   `"variance"`, `"lp"`, and `"covariate"`. The point-level matrices follow
#'   the same orientation as `log_pointpred`: rows are draws/replications and
#'   columns are observations. For frequentist models, a single-row matrix is
#'   returned.
#' @param ndraws Optional number of posterior draws for Bayesian models.
#' @param draw_ids Optional posterior draw indices for Bayesian models.
#' @param ... Additional arguments passed to class-specific methods.
#'
#' @return A named list containing the original Z-residual metadata and,
#'   optionally, point-level predictive quantities.
#'
#' @export
Zcov <- function(fit,
                 data,
                 type = NULL,
                 point_details = NULL,
                 ndraws = NULL,
                 draw_ids = NULL,
                 ...) {
  UseMethod("Zcov")
}

#' @method Zcov default
#' @export
Zcov.default <- function(fit,
                         data = NULL,
                         type = NULL,
                         point_details = NULL,
                         ndraws = NULL,
                         draw_ids = NULL,
                         ...) {
  stop(sprintf(
    "Zcov: unsupported model class: %s",
    paste(class(fit), collapse = "/")
  ), call. = FALSE)
}

# ---- shared helpers ----

.Zcov_normalize_point_details <- function(point_details) {
  if (is.null(point_details)) return(NULL)

  point_details <- tolower(as.character(point_details))
  point_details[point_details == "var"] <- "variance"
  point_details[point_details == "lp_point"] <- "lp"
  point_details <- unique(point_details)

  allowed <- c("mean", "variance", "lp", "covariate")
  bad <- setdiff(point_details, allowed)
  if (length(bad) > 0L) {
    stop(
      sprintf(
        "Zcov: unsupported point_details value(s): %s. Supported values are: %s.",
        paste(bad, collapse = ", "),
        paste(allowed, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  point_details
}

.Zcov_matrix_1xn <- function(x) {
  matrix(as.numeric(x), nrow = 1L)
}

.Zcov_warn_drop_na <- function(where, data_n, keep, mf) {
  if (all(keep)) return(invisible(NULL))

  drop_id <- which(!keep)
  drop_n <- length(drop_id)

  na_cols <- names(mf)[vapply(mf, function(x) any(is.na(x)), logical(1))]

  show_k <- 20L
  show_rows <- if (drop_n <= show_k) {
    paste(drop_id, collapse = ",")
  } else {
    paste0(
      paste(drop_id[seq_len(show_k)], collapse = ","),
      ", ... (+", drop_n - show_k, " more)"
    )
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
      where,
      drop_n,
      data_n,
      show_rows,
      show_r,
      if (nzchar(row_detail)) row_detail else "<none>",
      if (length(na_cols)) paste(na_cols, collapse = ",") else "<unknown>"
    ),
    call. = FALSE
  )

  invisible(NULL)
}

.Zcov_safe_model_matrix <- function(terms_or_formula, data, where = "Zcov") {
  out <- tryCatch(
    stats::model.matrix(terms_or_formula, data = data),
    error = function(e) NULL
  )

  if (is.null(out)) {
    warning(
      sprintf("%s: failed to construct model matrix; returning numeric(0) covariate matrix.", where),
      call. = FALSE
    )
    out <- matrix(numeric(0), nrow = nrow(data), ncol = 0L)
  }

  out
}

.Zcov_safe_covariate_matrix <- function(covariates, where = "Zcov") {
  if (is.null(covariates)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  
  covariates <- as.data.frame(covariates)
  
  if (ncol(covariates) < 1L) {
    return(matrix(numeric(0), nrow = nrow(covariates), ncol = 0L))
  }
  
  mm <- tryCatch(
    stats::model.matrix(~ . - 1, data = covariates),
    error = function(e) NULL
  )
  
  if (!is.null(mm)) {
    return(mm)
  }
  
  numeric_cols <- vapply(
    covariates,
    function(x) is.numeric(x) || is.integer(x) || is.logical(x),
    logical(1L)
  )
  
  if (!any(numeric_cols)) {
    warning(
      sprintf(
        "%s: failed to construct covariate matrix and no numeric covariates are available; returning numeric(0) covariate matrix.",
        where
      ),
      call. = FALSE
    )
    
    return(matrix(numeric(0), nrow = nrow(covariates), ncol = 0L))
  }
  
  warning(
    sprintf(
      "%s: failed to construct full covariate matrix; returning numeric covariates only.",
      where
    ),
    call. = FALSE
  )

  as.matrix(covariates[, numeric_cols, drop = FALSE])
}

.Zcov_type_normalize <- function(type, is_hurdle = FALSE) {
  if (is.null(type)) return(if (is_hurdle) "hurdle" else NULL)

  type_req <- tolower(as.character(type)[1L])
  if (type_req %in% c("overall", "whole", "all", "model")) type_req <- "hurdle"
  if (type_req %in% c("logistic", "logit")) type_req <- "zero"

  type_req
}

.Zcov_surv_response_name <- function(fit, mf) {
  rn <- NULL
  if (!is.null(mf) && ncol(mf) >= 1L) rn <- names(mf)[1L]
  if (!is.null(rn) && nzchar(rn)) return(rn)

  f <- tryCatch(stats::formula(fit), error = function(e) NULL)
  if (!is.null(f) && length(f) >= 2L) return(deparse(f[[2L]]))

  "Surv(...)"
}

.Zcov_surv_arg_names <- function(fit) {
  f <- tryCatch(stats::formula(fit), error = function(e) NULL)
  if (is.null(f) || length(f) < 2L) return(list())

  lhs <- f[[2L]]
  if (!is.call(lhs) || !identical(lhs[[1L]], as.name("Surv"))) return(list())

  args <- as.list(lhs)[-1L]
  out <- list(surv_expr = deparse(lhs), surv_args = lapply(args, deparse))

  if (length(args) == 2L) {
    out$time_name <- deparse(args[[1L]])
    out$status_name <- deparse(args[[2L]])
  } else if (length(args) >= 3L) {
    out$start_name <- deparse(args[[1L]])
    out$stop_name <- deparse(args[[2L]])
    out$status_name <- deparse(args[[3L]])
  }

  out
}
