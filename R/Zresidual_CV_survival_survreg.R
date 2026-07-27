#' Cross-validated Z-residuals for parametric survival models
#'
#' Internal \code{survreg} backend used by \code{Zresidual_CV()}. Each fold is
#' refitted on its training observations and evaluated on its held-out
#' observations. The current predictive backend supports right-censored
#' \code{Surv} responses.
#'
#' @param object A fitted \code{survival::survreg} object.
#' @param data Data containing all variables used in \code{object}.
#' @param nfolds Number of folds when \code{foldlist} is not supplied.
#' @param foldlist Optional list of held-out indices relative to the internally
#'   retained complete-case data.
#' @param nrep Number of randomized residual replicates.
#' @param randomized Logical; generate randomized residuals for censored
#'   observations when \code{TRUE}.
#' @param type Optional model subtype passed to \code{Zresidual()}.
#' @param log_pointpred Optional predictive backend function or function name.
#' @param seed Optional integer seed.
#' @param ... Additional arguments passed to \code{Zresidual()}.
#'
#' @return A \code{cvzresid} matrix.
#'
#' @keywords internal
Zresidual_CV_survival_survreg <- function(object,
                                          data,
                                          nfolds = NULL,
                                          foldlist = NULL,
                                          nrep = 1L,
                                          randomized = TRUE,
                                          type = NULL,
                                          log_pointpred = NULL,
                                          seed = NULL,
                                          ...) {
  caller <- "Zresidual_CV_survival_survreg"

  if (!is.null(seed)) set.seed(seed)

  if (missing(object) || is.null(object) || !inherits(object, "survreg")) {
    stop(paste0(caller, ": `object` must inherit from 'survreg'."), call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop(paste0(caller, ": `data` must be provided."), call. = FALSE)
  }

  data <- as.data.frame(data)

  nrep_value <- suppressWarnings(as.numeric(nrep))
  if (length(nrep_value) != 1L || !is.finite(nrep_value) ||
      nrep_value < 1 || nrep_value != floor(nrep_value)) {
    stop(paste0(caller, ": `nrep` must be a positive integer."), call. = FALSE)
  }
  nrep <- as.integer(nrep_value)

  if (length(randomized) != 1L || is.na(randomized) || !is.logical(randomized)) {
    stop(paste0(caller, ": `randomized` must be TRUE or FALSE."), call. = FALSE)
  }
  nrep_eff <- if (isTRUE(randomized)) nrep else 1L

  mf_all <- stats::model.frame(
    object$terms,
    data = data,
    na.action = stats::na.pass
  )
  keep <- stats::complete.cases(mf_all)

  if (!all(keep)) {
    warning(
      sprintf("%s: removed %d rows with missing model variables before cross-validation.", caller, sum(!keep)),
      call. = FALSE
    )
  }

  data_used <- data[keep, , drop = FALSE]
  original_row_index <- which(keep)
  n <- nrow(data_used)

  if (n < 2L) {
    stop(paste0(caller, ": at least two complete observations are required."), call. = FALSE)
  }

  mf <- stats::model.frame(
    object$terms,
    data = data_used,
    na.action = stats::na.pass
  )
  y <- stats::model.response(mf)

  if (!inherits(y, "Surv")) {
    stop(paste0(caller, ": the response must be a Surv object."), call. = FALSE)
  }

  surv_type <- attr(y, "type")
  if (is.null(surv_type) || !identical(surv_type, "right")) {
    stop(
      sprintf(
        "%s: only right-censored Surv responses are currently supported; received type '%s'.",
        caller,
        as.character(surv_type %||% "unknown")
      ),
      call. = FALSE
    )
  }

  ymat <- as.matrix(y)
  if (ncol(ymat) != 2L) {
    stop(paste0(caller, ": unexpected right-censored Surv response structure."), call. = FALSE)
  }

  time <- ymat[, 1L]
  status <- as.integer(ymat[, 2L])

  if (anyNA(status) || any(!status %in% c(0L, 1L))) {
    stop(paste0(caller, ": event status must be coded as 0/1."), call. = FALSE)
  }

  covariates <- mf[, -1L, drop = FALSE]

  normalize_foldlist <- function(x, n) {
    if (!is.list(x) || length(x) < 2L) {
      stop(paste0(caller, ": `foldlist` must contain at least two folds."), call. = FALSE)
    }

    out <- lapply(x, function(index) {
      value <- suppressWarnings(as.numeric(index))
      if (length(value) == 0L || any(!is.finite(value)) ||
          any(value < 1) || any(value != floor(value))) {
        stop(paste0(caller, ": every fold must contain positive integer indices."), call. = FALSE)
      }
      sort(unique(as.integer(value)))
    })

    flat <- unlist(out, use.names = FALSE)
    if (anyDuplicated(flat) || !identical(sort(flat), seq_len(n))) {
      stop(
        paste0(caller, ": `foldlist` must be a non-overlapping partition of 1:nrow(data_used)."),
        call. = FALSE
      )
    }

    names(out) <- sprintf("Fold%02d", seq_along(out))
    out
  }

  if (is.null(foldlist)) {
    if (is.null(nfolds)) nfolds <- min(10L, n)
    nfolds <- .cv_validate_k(nfolds, n, caller = caller)

    foldlist <- make_fold(
      fix_var = covariates,
      y = y,
      k = nfolds,
      censor = status
    )
  } else {
    foldlist <- normalize_foldlist(foldlist, n)
  }

  z_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  rsp_full <- matrix(NA_real_, nrow = n, ncol = nrep_eff)
  lp_full <- rep(NA_real_, n)
  fold_id <- integer(n)
  extra_args <- list(...)

  for (fid in seq_along(foldlist)) {
    test_index <- foldlist[[fid]]
    train_data <- data_used[-test_index, , drop = FALSE]
    test_data <- data_used[test_index, , drop = FALSE]

    fit_train <- tryCatch(
      stats::update(object, data = train_data),
      error = function(e) {
        stop(
          sprintf("%s: fold %d model refitting failed: %s", caller, fid, conditionMessage(e)),
          call. = FALSE
        )
      }
    )

    z_call <- c(
      list(
        fit = fit_train,
        data = test_data,
        log_pointpred = log_pointpred,
        type = type,
        randomized = randomized,
        nrep = nrep_eff
      ),
      extra_args
    )

    z_fold <- tryCatch(
      do.call(Zresidual, z_call),
      error = function(e) {
        stop(
          sprintf("%s: fold %d Zresidual calculation failed: %s", caller, fid, conditionMessage(e)),
          call. = FALSE
        )
      }
    )

    z_matrix <- as.matrix(z_fold)
    expected_dim <- c(length(test_index), nrep_eff)

    if (!identical(dim(z_matrix), expected_dim)) {
      stop(
        sprintf(
          "%s: fold %d returned Z-residual dimension %s; expected %s.",
          caller,
          fid,
          paste(dim(z_matrix), collapse = " x "),
          paste(expected_dim, collapse = " x ")
        ),
        call. = FALSE
      )
    }

    if (anyNA(z_matrix) || any(!is.finite(z_matrix))) {
      stop(sprintf("%s: fold %d returned non-finite Z-residuals.", caller, fid), call. = FALSE)
    }

    rsp_fold <- attr(z_fold, "rsp")
    if (is.null(rsp_fold)) {
      stop(sprintf("%s: fold %d Z-residual output is missing the `rsp` attribute.", caller, fid), call. = FALSE)
    }

    rsp_matrix <- as.matrix(rsp_fold)
    if (!identical(dim(rsp_matrix), expected_dim) || anyNA(rsp_matrix) ||
        any(!is.finite(rsp_matrix))) {
      stop(sprintf("%s: fold %d returned invalid probability-scale residuals.", caller, fid), call. = FALSE)
    }

    lp_fold <- attr(z_fold, "linear_pred") %||% attr(z_fold, "linear.pred")
    if (is.null(lp_fold)) {
      lp_fold <- tryCatch(
        as.vector(stats::predict(fit_train, newdata = test_data, type = "lp")),
        error = function(e) NULL
      )
    }

    if (is.null(lp_fold) || length(lp_fold) != length(test_index) ||
        anyNA(lp_fold) || any(!is.finite(lp_fold))) {
      stop(sprintf("%s: fold %d returned an invalid linear predictor.", caller, fid), call. = FALSE)
    }

    z_full[test_index, ] <- z_matrix
    rsp_full[test_index, ] <- rsp_matrix
    lp_full[test_index] <- as.numeric(lp_fold)
    fold_id[test_index] <- fid
  }

  if (anyNA(z_full) || any(!is.finite(z_full))) {
    stop(paste0(caller, ": combined CV Z-residual matrix is incomplete."), call. = FALSE)
  }
  if (anyNA(rsp_full) || any(!is.finite(rsp_full))) {
    stop(paste0(caller, ": combined CV probability-scale residual matrix is incomplete."), call. = FALSE)
  }
  if (anyNA(lp_full) || any(!is.finite(lp_full))) {
    stop(paste0(caller, ": combined CV linear predictor is incomplete."), call. = FALSE)
  }

  rownames(z_full) <- rownames(data_used)
  colnames(z_full) <- paste0("rep", seq_len(nrep_eff))
  rownames(rsp_full) <- rownames(data_used)
  colnames(rsp_full) <- colnames(z_full)
  class(z_full) <- c("cvzresid", "zresid", "matrix", "array")

  attr(z_full, "rsp") <- rsp_full
  attr(z_full, "type") <- "survival"
  attr(z_full, "linear.pred") <- lp_full
  attr(z_full, "linear_pred") <- lp_full
  attr(z_full, "censored.status") <- as.integer(status == 0L)
  attr(z_full, "covariates") <- covariates
  attr(z_full, "object.model.frame") <- mf
  attr(z_full, "foldlist") <- foldlist
  attr(z_full, "fold_id") <- fold_id
  attr(z_full, "original_row_index") <- original_row_index
  attr(z_full, "surv_type") <- surv_type
  attr(z_full, "cv_scope") <- "observation_level"

  attr(z_full, "zcov") <- list(
    type = type,
    family = "survreg",
    response_name = names(mf)[1L],
    response = y,
    covariates = covariates,
    linear_pred = lp_full,
    obs_id = seq_len(n),
    y_type = ifelse(status == 1L, 1L, 0L),
    y_type_kind = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra = list(
      time = time,
      status = status,
      surv_type = surv_type,
      dist = object$dist,
      cv = TRUE,
      cv_scope = "observation_level",
      foldlist = foldlist,
      fold_id = fold_id,
      original_row_index = original_row_index
    )
  )

  z_full
}
