# =========================================================
# Internal cross-validation fold helpers
# =========================================================

.cv_validate_k <- function(k, n, caller = "cross-validation") {
  k <- suppressWarnings(as.integer(k))

  if (length(k) != 1L || is.na(k) || k < 2L) {
    stop(sprintf("%s: `k` must be an integer >= 2.", caller), call. = FALSE)
  }

  if (n < 2L) {
    stop(sprintf("%s: at least two observations are required.", caller), call. = FALSE)
  }

  if (k > n) {
    stop(sprintf("%s: `k` cannot exceed the number of observations.", caller), call. = FALSE)
  }

  k
}

.cv_quantile_bin <- function(x, max_bins = 5L) {
  x <- as.numeric(x)
  n <- length(x)

  if (n < 2L || length(unique(x[is.finite(x)])) < 2L) {
    return(factor(rep("all", n)))
  }

  n_bins <- min(
    as.integer(max_bins),
    max(2L, floor(sqrt(n))),
    length(unique(x[is.finite(x)]))
  )

  probs <- seq(0, 1, length.out = n_bins + 1L)
  breaks <- unique(stats::quantile(x, probs = probs, na.rm = TRUE, type = 8))

  if (length(breaks) < 3L) {
    return(factor(rep("all", n)))
  }

  cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)
}

.cv_strata <- function(y, k, censor = NULL) {
  n <- NROW(y)

  if (inherits(y, "Surv")) {
    ymat <- as.matrix(y)
    time <- ymat[, ncol(ymat) - 1L]
    status <- if (is.null(censor)) ymat[, ncol(ymat)] else censor
    time_bin <- .cv_quantile_bin(time, max_bins = min(5L, k))
    return(interaction(factor(status), time_bin, drop = TRUE, lex.order = TRUE))
  }

  if (!is.null(censor)) {
    if (length(censor) != n) {
      stop("kfold_fn: `censor` must have the same length as `y`.", call. = FALSE)
    }

    if (is.numeric(y)) {
      y_bin <- .cv_quantile_bin(y, max_bins = min(5L, k))
    } else {
      y_bin <- factor(y)
    }

    return(interaction(factor(censor), y_bin, drop = TRUE, lex.order = TRUE))
  }

  if (is.numeric(y)) {
    return(.cv_quantile_bin(y, max_bins = min(5L, k)))
  }

  factor(y)
}

.cv_assign_stratified_folds <- function(strata, k) {
  n <- length(strata)
  fold_id <- integer(n)
  fold_size <- integer(k)

  strata_index <- split(seq_len(n), strata, drop = TRUE)
  strata_index <- strata_index[order(lengths(strata_index), decreasing = TRUE)]

  for (idx in strata_index) {
    idx <- idx[
      sample.int(length(idx), size = length(idx), replace = FALSE)
    ]

    for (obs in idx) {
      smallest <- which(fold_size == min(fold_size))
      fid <- if (length(smallest) == 1L) smallest else sample(smallest, 1L)
      fold_id[obs] <- fid
      fold_size[fid] <- fold_size[fid] + 1L
    }
  }

  if (any(fold_size == 0L)) {
    stop("kfold_fn: unable to create non-empty folds.", call. = FALSE)
  }

  fold_id
}

#' Create stratified k-fold indices
#'
#' Generate non-empty cross-validation folds. For survival outcomes,
#' stratification uses censoring/event status together with quantile groups of
#' observed follow-up time. Continuous outcomes are stratified using quantile
#' groups rather than treating every distinct value as a factor level.
#'
#' @param y Outcome used for stratification. It may be a \code{Surv} object,
#'   factor, character vector, logical vector, or numeric vector.
#' @param k Number of folds.
#' @param list Logical; return a list of indices when \code{TRUE}, otherwise an
#'   integer vector of fold labels.
#' @param returnTrain Logical; when \code{list = TRUE}, return training indices
#'   instead of test indices.
#' @param censor Optional censoring/status vector used with non-\code{Surv}
#'   outcomes.
#'
#' @return A list of index vectors or an integer vector of fold labels.
#'
#' @keywords internal
#' @noRd
kfold_fn <- function(y,
                     k,
                     list = TRUE,
                     returnTrain = FALSE,
                     censor = NULL) {
  n <- NROW(y)
  k <- .cv_validate_k(k, n, caller = "kfold_fn")

  if (!is.null(censor) && length(censor) != n) {
    stop("kfold_fn: `censor` must have the same length as `y`.", call. = FALSE)
  }

  strata <- .cv_strata(y = y, k = k, censor = censor)
  fold_id <- .cv_assign_stratified_folds(strata = strata, k = k)

  if (!isTRUE(list)) return(fold_id)

  out <- lapply(seq_len(k), function(fid) which(fold_id == fid))
  names(out) <- sprintf("Fold%02d", seq_len(k))

  if (isTRUE(returnTrain)) {
    all_index <- seq_len(n)
    out <- lapply(out, function(test_index) setdiff(all_index, test_index))
  }

  out
}

.cv_categorical_columns <- function(fix_var) {
  if (is.null(fix_var) || NCOL(fix_var) == 0L) return(character(0))

  fix_var <- as.data.frame(fix_var)
  names(fix_var)[vapply(
    fix_var,
    function(x) is.factor(x) || is.character(x) || is.logical(x),
    logical(1L)
  )]
}

.cv_validate_categorical_feasibility <- function(fix_var, categorical_names) {
  if (length(categorical_names) == 0L) return(invisible(TRUE))

  for (name in categorical_names) {
    x <- as.character(fix_var[[name]])
    counts <- table(x, useNA = "no")
    rare <- names(counts)[counts < 2L]

    if (length(rare) > 0L) {
      stop(
        sprintf(
          paste0(
            "make_fold: categorical variable '%s' contains level(s) with fewer ",
            "than two observations: %s. Every test-fold level must also be ",
            "present in its corresponding training fold."
          ),
          name,
          paste(rare, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.cv_fold_has_known_levels <- function(fold_list, fix_var, categorical_names) {
  if (length(categorical_names) == 0L) return(TRUE)

  n <- nrow(fix_var)
  all_index <- seq_len(n)

  for (test_index in fold_list) {
    train_index <- setdiff(all_index, test_index)

    for (name in categorical_names) {
      test_levels <- unique(as.character(fix_var[[name]][test_index]))
      train_levels <- unique(as.character(fix_var[[name]][train_index]))
      test_levels <- test_levels[!is.na(test_levels)]
      train_levels <- train_levels[!is.na(train_levels)]

      if (!all(test_levels %in% train_levels)) return(FALSE)
    }
  }

  TRUE
}

#' Build cross-validation folds with categorical-level checks
#'
#' Generate stratified test folds and ensure that every categorical level in a
#' test fold is represented in the corresponding training fold. Frailty group
#' restrictions are checked separately by the Cox frailty CV backend.
#'
#' @param fix_var Matrix or data frame of covariates.
#' @param y Outcome used for stratification.
#' @param k Number of folds.
#' @param censor Censoring/status vector with the same number of observations as
#'   \code{y}.
#' @param max_attempts Maximum number of candidate fold assignments.
#'
#' @return A list of test-fold indices.
#'
#' @keywords internal
#' @noRd
make_fold <- function(fix_var,
                      y,
                      k,
                      censor,
                      max_attempts = 200L) {
  fix_var <- as.data.frame(fix_var)
  n <- NROW(y)
  k <- .cv_validate_k(k, n, caller = "make_fold")

  if (nrow(fix_var) != n) {
    stop("make_fold: `fix_var` and `y` must contain the same number of observations.", call. = FALSE)
  }

  if (length(censor) != n) {
    stop("make_fold: `censor` must have the same number of observations as `y`.", call. = FALSE)
  }

  max_attempts <- suppressWarnings(as.integer(max_attempts))
  if (length(max_attempts) != 1L || is.na(max_attempts) || max_attempts < 1L) {
    stop("make_fold: `max_attempts` must be a positive integer.", call. = FALSE)
  }

  categorical_names <- .cv_categorical_columns(fix_var)
  .cv_validate_categorical_feasibility(fix_var, categorical_names)

  for (attempt in seq_len(max_attempts)) {
    fold_list <- kfold_fn(
      y = y,
      k = k,
      list = TRUE,
      returnTrain = FALSE,
      censor = censor
    )

    if (.cv_fold_has_known_levels(fold_list, fix_var, categorical_names)) {
      return(fold_list)
    }
  }

  stop(
    sprintf(
      paste0(
        "make_fold: unable to construct %d folds after %d attempts while ",
        "keeping every test-fold categorical level in its training fold. ",
        "Reduce `k`, combine rare levels, or supply a valid `foldlist`."
      ),
      k,
      max_attempts
    ),
    call. = FALSE
  )
}
