# =========================================================
# Helpers for predictive checks
# =========================================================

.infer_predcheck_backend <- function(fit) {
  if (inherits(fit, "brmsfit")) {
    fam <- fit$family$family
    return(paste0("predcheck_pointpred_brms_", fam))
  }

  stop(
    "Cannot infer predictive-check backend from `fit`. ",
    "Please supply `predcheck_pointpred` explicitly."
  )
}

.validate_predcheck_object <- function(obj) {
  needed <- c("support", "family", "y", "n", "ndraws", "pmf", "tail", "rng", "moments")
  miss <- setdiff(needed, names(obj))

  if (length(miss) > 0) {
    stop("Invalid predcheck object. Missing fields: ", paste(miss, collapse = ", "))
  }

  if (!is.function(obj$pmf)) stop("`pmf` must be a function.")
  if (!is.function(obj$tail)) stop("`tail` must be a function.")
  if (!is.function(obj$rng)) stop("`rng` must be a function.")
  if (!is.function(obj$moments)) stop("`moments` must be a function.")

  invisible(obj)
}

.predcheck_clip_prob <- function(p, eps = .Machine$double.eps) {
  pmin(pmax(p, eps), 1 - eps)
}

.predcheck_get_draw_ids <- function(fit, ndraws = NULL) {
  if (!inherits(fit, "brmsfit")) {
    if (is.null(ndraws)) {
      return(NULL)
    }

    ndraws <- suppressWarnings(as.integer(ndraws))
    if (!is.finite(ndraws) || ndraws < 1L) {
      stop("`ndraws` must be a positive integer.")
    }

    return(seq_len(ndraws))
  }

  nd_all <- brms::ndraws(fit)

  if (is.null(ndraws)) {
    return(seq_len(nd_all))
  }

  ndraws <- suppressWarnings(as.integer(ndraws))
  if (!is.finite(ndraws) || ndraws < 1L) {
    stop("`ndraws` must be a positive integer.")
  }

  nd_use <- min(ndraws, nd_all)
  sort(sample.int(nd_all, size = nd_use, replace = FALSE))
}

.predcheck_extract_response <- function(fit, data, eps = sqrt(.Machine$double.eps)) {
  if (is.null(data)) {
    stop("`data` must be provided.")
  }

  formula_obj <- fit$formula$formula
  mf <- stats::model.frame(formula_obj, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)

  if (inherits(y, "Surv") || is.matrix(y)) {
    stop("Predictive checks currently require a univariate count response.")
  }

  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop("The response contains NA or non-finite values.")
  }

  if (any(y < 0)) {
    stop("The response must be non-negative.")
  }

  if (any(abs(y - round(y)) > eps)) {
    stop("The response must contain integer counts.")
  }

  as.integer(round(y))
}

.predcheck_extract_x <- function(data, x, n_expected) {
  if (is.null(x)) {
    return(NULL)
  }

  xv <- if (is.character(x) && length(x) == 1L) {
    if (is.null(data) || !(x %in% names(data))) {
      stop("`x` is not a column in `data`.")
    }
    data[[x]]
  } else {
    x
  }

  if (length(xv) != n_expected) {
    stop("`x` must have length equal to nrow(data).")
  }

  if (is.factor(xv)) {
    return(droplevels(xv))
  }

  if (is.character(xv)) {
    xv_num <- suppressWarnings(as.numeric(xv))
    xv <- if (all(!is.na(xv_num))) xv_num else factor(xv)
  }

  if (is.logical(xv)) {
    xv <- as.integer(xv)
  }

  if (is.numeric(xv)) {
    if (anyNA(xv) || any(!is.finite(xv))) {
      stop("`x` contains NA or non-finite numeric values.")
    }
    return(xv)
  }

  if (anyNA(xv)) {
    stop("`x` contains NA values.")
  }

  droplevels(as.factor(xv))
}

.predcheck_extract_shape_matrix <- function(fit,
                                            data,
                                            draw_ids,
                                            n_obs,
                                            eps = .Machine$double.eps) {
  shape_raw <- tryCatch(
    brms::posterior_epred(
      fit,
      newdata = data,
      dpar = "shape",
      draw_ids = draw_ids
    ),
    error = function(e) NULL
  )

  if (!is.null(shape_raw)) {
    if (is.null(dim(shape_raw))) {
      shape_mat <- matrix(shape_raw, nrow = length(draw_ids), ncol = n_obs)
    } else {
      shape_mat <- shape_raw
    }
  } else {
    dd <- posterior::as_draws_df(fit)
    shape_col <- intersect(c("shape", "shape_Intercept", "b_shape_Intercept"), names(dd))

    if (length(shape_col) == 0) {
      stop("Cannot find shape draws in `fit`.")
    }

    shape_vec <- as.numeric(dd[[shape_col[1]]])
    if (length(shape_vec) < length(draw_ids)) {
      shape_vec <- sample(shape_vec, size = length(draw_ids), replace = TRUE)
    } else {
      shape_vec <- shape_vec[draw_ids]
    }

    shape_mat <- matrix(shape_vec, nrow = length(draw_ids), ncol = n_obs)
  }

  pmax(shape_mat, eps)
}
