# =========================================================
# Top functions
# =========================================================

#' Predictive checks
#'
#' @description
#' Performs predictive checks on evaluation data. The same computational
#' algorithm is used for in-sample and holdout evaluation.
#'
#' @param fit A fitted model object.
#' @param data Data used to evaluate the fitted model.
#' @param data_role Role of \code{data}: \code{"training"} if the data were used
#'   to fit the model, or \code{"holdout"} otherwise.
#' @param log_pointpred Optional pointwise predictive bridge.
#' @param point_Zcov Optional point-level covariate/moment bridge.
#' @param postpred_simulate Optional posterior predictive simulation bridge.
#' @param test Optional list of Z-residual test specifications.
#' @param x Optional covariate for the default ANOVA-style diagnostic.
#' @param ndraws Optional number of posterior draws to use.
#' @param seed Optional random seed.
#' @param type Optional component selector.
#' @param randomized Logical; use randomized residuals for discrete outcomes.
#' @param k_anova Maximum number of bins for numeric ANOVA covariates.
#' @param k_bl Maximum number of bins for numeric Bartlett covariates.
#' @param ... Additional arguments passed to bridge functions.
#' @return An object of class \code{"predcheck"}.
#' @examples
#' \dontrun{
#' training_result <- predcheck(
#'   fit_full, full_data, data_role = "training"
#' )
#' holdout_result <- predcheck(
#'   fit_training, test_data, data_role = "holdout"
#' )
#' }
#' @export
predcheck <- function(fit, data,
                      data_role = c("training", "holdout"),
                      log_pointpred = NULL,
                      point_Zcov = NULL,
                      postpred_simulate = NULL,
                      test = NULL, x = NULL, ndraws = NULL, seed = NULL,
                      type = NULL, randomized = TRUE,
                      k_anova = 10, k_bl = 10, ...) {
  data_role <- match.arg(data_role)
  .predictive_check_core(
    fit = fit, data = data, data_role = data_role,
    log_pointpred = log_pointpred, point_Zcov = point_Zcov,
    postpred_simulate = postpred_simulate, test = test, x = x,
    ndraws = ndraws, seed = seed, type = type,
    randomized = randomized, k_anova = k_anova, k_bl = k_bl, ...
  )
}

#' Print a predictive-check summary
#'
#' @param x An object of class \code{"predcheck"}.
#' @param ... Reserved for future use.
#'
#' @return The input object, invisibly.
#'
#' @export
print.predcheck <- function(x, ...) {
  cat("Predictive check\n")
  evaluation_label <- if (identical(x$data_role, "training")) {
    "in-sample (PPC)"
  } else {
    "holdout (HPC)"
  }
  cat("  evaluation: ", evaluation_label, "\n", sep = "")
  cat("  family:     ", x$family, "\n", sep = "")
  if (!is.null(x$type) && !is.na(x$type)) cat("  type:   ", x$type, "\n", sep = "")
  cat("  n:      ", x$n, "\n", sep = "")
  cat("  ndraws: ", x$ndraws, "\n\n", sep = "")

  mode_label <- if (identical(x$data_role, "training")) "PPC" else "HPC"

  get_results <- function(obj, name) {
    if (is.null(obj) || is.null(obj$results)) return(list())
    obj$results[vapply(
      obj$results,
      function(r) identical(tolower(as.character(r$name)[1L]), name),
      logical(1L)
    )]
  }

  fmt_label <- function(res, default) {
    label <- res$label %||% default
    as.character(label)[1L]
  }

  cat("  chi-square ", mode_label, " p-value:                 ",
      format(x$chisq, digits = 4), "\n", sep = "")
  cat("  Z-SW randomized mean observed p-value: ",
      format(x$mean_sw_pv, digits = 4), "\n", sep = "")
  cat("  Z-SW midpoint ", mode_label, " p-value:        ",
      format(x$sw_mpp, digits = 4), "\n", sep = "")

  aov_results <- get_results(x$zresid, "aov")
  if (length(aov_results) > 0L) {
    for (res in aov_results) {
      cat("  Z-", fmt_label(res, "AOV"), " randomized ", mode_label, " p-value:      ",
          format(res$p_value, digits = 4), "\n", sep = "")
    }
  } else if (!is.null(x$aov_rpp) && any(!is.na(x$aov_rpp))) {
    cat("  Z-AOV randomized ", mode_label, " p-value:      ",
        format(x$aov_rpp, digits = 4), "\n", sep = "")
  }

  bl_results <- get_results(x$zresid, "bartlett")
  if (length(bl_results) > 0L) {
    for (res in bl_results) {
      cat("  Z-", fmt_label(res, "BL"), " randomized ", mode_label, " p-value:       ",
          format(res$p_value, digits = 4), "\n", sep = "")
    }
  } else if (!is.null(x$bartlett_rpp) && any(!is.na(x$bartlett_rpp))) {
    cat("  Z-BL randomized ", mode_label, " p-value:       ",
        format(x$bartlett_rpp, digits = 4), "\n", sep = "")
  }

  invisible(x)
}

#' Print a Z-residual predictive-check summary
#'
#' @param x An object of class \code{"z_predcheck"}.
#' @param ... Reserved for future use.
#'
#' @return The input object, invisibly.
#'
#' @export
print.z_predcheck <- function(x, ...) {
  cat("Z-residual predictive check\n")
  evaluation_label <- if (identical(x$data_role, "training")) {
    "in-sample (PPC)"
  } else {
    "holdout (HPC)"
  }
  cat("  evaluation: ", evaluation_label, "\n", sep = "")
  cat("  family:     ", x$family, "\n", sep = "")
  if (!is.null(x$type) && !is.na(x$type)) cat("  type:       ", x$type, "\n", sep = "")
  cat("  n:          ", x$n, "\n", sep = "")
  cat("  ndraws:     ", x$ndraws, "\n", sep = "")
  cat("  randomized: ", x$randomized, "\n\n", sep = "")

  mode_label <- if (identical(x$data_role, "training")) "PPC" else "HPC"

  for (res in x$results) {
    label <- res$label %||% res$name
    if (isTRUE(x$randomized)) {
      cat("  ", label,
          ": mean observed p = ", format(res$mean_obs_p, digits = 4),
          "\n", sep = "")
    } else {
      cat("  ", label,
          ": ", mode_label, " p = ", format(res$p_value, digits = 4),
          ", mean observed p = ", format(res$mean_obs_p, digits = 4),
          "\n", sep = "")
    }
  }

  invisible(x)
}

#' Print a moment predictive-check summary
#'
#' @param x An object of class \code{"moment_predcheck"}.
#' @param ... Reserved for future use.
#'
#' @return The input object, invisibly.
#' @export
print.moment_predcheck <- function(x, ...) {
  cat("Moment predictive check\n")
  evaluation_label <- if (identical(x$data_role, "training")) {
    "in-sample (PPC)"
  } else {
    "holdout (HPC)"
  }
  cat("  evaluation: ", evaluation_label, "\n", sep = "")
  cat("  family:     ", x$family, "\n", sep = "")
  if (!is.null(x$type) && !is.na(x$type)) cat("  type:   ", x$type, "\n", sep = "")
  cat("  n:      ", x$n, "\n", sep = "")
  cat("  ndraws: ", x$ndraws, "\n\n", sep = "")
  cat("  p-value: ", format(x$p_value, digits = 4), "\n", sep = "")
  invisible(x)
}
