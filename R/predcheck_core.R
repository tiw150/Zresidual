#' Z-residual predictive checks
#'
#' @description
#' Computes probability-space predictive diagnostics from pointwise predictive
#' log-survival and log-density/probability matrices returned by
#' \code{log_pointpred}. The final PPC/HPC p-value is computed by comparing
#' the test statistic/p-value from posterior predictive replicated data against
#' the same quantity from the observed data, using the same posterior draw.
#'
#' @param fit A fitted model object.
#' @param data Data used to evaluate the fitted model.
#' @param data_role Role of \code{data}: \code{"training"} if the data were used
#'   to fit the model, or \code{"holdout"} otherwise.
#' @param test A list of test specifications. Supported aliases include
#'   \code{"sw"}, \code{"aov"}, and \code{"bartlett"}.
#' @param log_pointpred Optional pointwise predictive bridge.
#' @param point_Zcov Optional point-level covariate bridge.
#' @param postpred_simulate Optional posterior predictive simulation bridge.
#' @param type Optional component selector.
#' @param x Optional covariate for the default ANOVA-style test.
#' @param ndraws Optional number of posterior draws.
#' @param draw_ids Optional posterior draw indices.
#' @param seed Optional random seed.
#' @param randomized Logical; if \code{TRUE}, use randomized residuals for
#'   discrete outcomes. If \code{FALSE}, use midpoint residuals with \code{U=0.5}.
#' @param k_anova Maximum number of bins for numeric ANOVA covariates.
#' @param k_bl Maximum number of bins for numeric Bartlett covariates.
#' @param eps Probability clipping constant.
#' @param ... Additional arguments passed to bridges.
#'
#' @return An object of class \code{"z_predcheck"}.
#'
#' @export
predcheck_zresid <- function(fit,
                             data,
                             data_role = c("training", "holdout"),
                             test = NULL,
                             log_pointpred = NULL,
                             point_Zcov = NULL,
                             postpred_simulate = NULL,
                             type = NULL,
                             x = NULL,
                             ndraws = NULL,
                             draw_ids = NULL,
                             seed = NULL,
                             randomized = TRUE,
                             k_anova = 10,
                             k_bl = 10,
                             eps = 1e-12,
                             ...) {
  if (!is.null(seed)) set.seed(seed)

  data_role <- match.arg(data_role)
  eval_data <- as.data.frame(data)

  if (is.null(draw_ids)) {
    draw_ids <- .predcheck_get_draw_ids(fit, ndraws = ndraws)
  }

  lp_fn <- .predcheck_resolve_log_pointpred(
    log_pointpred = log_pointpred,
    fit = fit,
    pred_method = "analytic",
    type = type
  )

  zcov_fn <- .predcheck_resolve_point_Zcov(point_Zcov)
  sim_fn <- .predcheck_resolve_postpred_simulate(postpred_simulate)

  log_pred_obs <- .predcheck_call_log_pointpred(
    fn = lp_fn,
    fit = fit,
    data = eval_data,
    type = type,
    draw_ids = draw_ids,
    ...
  )

  zcov_out <- zcov_fn(
    fit = fit,
    data = eval_data,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )

  z_obs_matrix <- .predcheck_compute_z_matrix(
    log_pred = log_pred_obs,
    randomized = randomized,
    eps = eps
  )

  y_rep <- sim_fn(
    fit = fit,
    data = eval_data,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )
  y_rep <- .predcheck_as_matrix(y_rep, "y_rep")

  N <- nrow(z_obs_matrix)
  M_z <- ncol(z_obs_matrix)
  M_y <- nrow(y_rep)

  if (ncol(y_rep) != N) {
    stop("predcheck_zresid: `y_rep` columns do not match the Z-residual length.", call. = FALSE)
  }

  if (M_z != M_y) {
    if (M_z == 1L) {
      z_obs_matrix <- matrix(rep(z_obs_matrix[, 1L], M_y), nrow = N, ncol = M_y)
      M_z <- M_y
    } else if (M_y == 1L) {
      y_rep <- matrix(rep(y_rep[1L, ], each = M_z), nrow = M_z, ncol = N, byrow = FALSE)
      M_y <- M_z
    } else {
      stop("predcheck_zresid: observed Z-residual draws and posterior predictive draws have incompatible dimensions.", call. = FALSE)
    }
  }

  M <- ncol(z_obs_matrix)
  tests <- .predcheck_normalize_tests(
    test = test,
    x = x,
    k_anova = k_anova,
    k_bl = k_bl
  )

  obs_p_mat <- matrix(NA_real_, nrow = M, ncol = length(tests))
  rep_p_mat <- matrix(NA_real_, nrow = M, ncol = length(tests))

  z_rep_matrix <- .predcheck_fast_rep_z_matrix(
    fit = fit,
    data = eval_data,
    zcov_out = zcov_out,
    y_rep = y_rep,
    type = type,
    draw_ids = draw_ids,
    randomized = randomized,
    eps = eps,
    ...
  )

  if (is.null(z_rep_matrix)) {
    stop(
      paste0(
        "predcheck_zresid: no matrix replicated Z-residual implementation is available ",
        "for this model/family/type. PPC/HPC no longer falls back to the old ",
        "per-draw log_pointpred loop. Add support in `.predcheck_fast_rep_z_matrix()`."
      ),
      call. = FALSE
    )
  }

  if (!identical(dim(z_rep_matrix), c(N, M))) {
    stop(
      sprintf(
        "predcheck_zresid: matrix replicated Z-residual path returned dimension %s but expected %s.",
        paste(dim(z_rep_matrix), collapse = " x "),
        paste(c(N, M), collapse = " x ")
      ),
      call. = FALSE
    )
  }

  .eval_test_matrix <- function(zmat) {
    out <- lapply(seq_along(tests), function(k) {
      spec <- tests[[k]]
      vapply(seq_len(M), function(j) {
        tryCatch(
          .predcheck_eval_one_test(
            z = zmat[, j],
            zcov_out = zcov_out,
            spec = spec,
            draw = j,
            n_draws = M
          ),
          error = function(e) NA_real_
        )
      }, numeric(1L))
    })

    out <- do.call(cbind, out)
    if (is.null(dim(out))) out <- matrix(out, ncol = length(tests))
    out
  }

  obs_p_mat <- .eval_test_matrix(z_obs_matrix)
  rep_p_mat <- .eval_test_matrix(z_rep_matrix)

  results <- vector("list", length(tests))

  for (k in seq_along(tests)) {
    spec <- tests[[k]]
    obs_p_values <- obs_p_mat[, k]
    rep_p_values <- rep_p_mat[, k]
    comparisons <- rep_p_values <= obs_p_values

    results[[k]] <- list(
      name = as.character(spec$name)[1L],
      label = .predcheck_test_label(spec),
      spec = spec,
      p_value = mean(comparisons, na.rm = TRUE),
      comparisons = comparisons,
      obs_p_values = obs_p_values,
      rep_p_values = rep_p_values,
      mean_obs_p = mean(obs_p_values, na.rm = TRUE),
      mean_rep_p = mean(rep_p_values, na.rm = TRUE),
      median_obs_p = stats::median(obs_p_values, na.rm = TRUE),
      min_obs_p = suppressWarnings(min(obs_p_values, na.rm = TRUE)),
      max_obs_p = suppressWarnings(max(obs_p_values, na.rm = TRUE)),
      mean_p = mean(obs_p_values, na.rm = TRUE),
      median_p = stats::median(obs_p_values, na.rm = TRUE),
      min_p = suppressWarnings(min(obs_p_values, na.rm = TRUE)),
      max_p = suppressWarnings(max(obs_p_values, na.rm = TRUE))
    )
  }

  out <- list(
    data_role = data_role,
    evaluation = if (identical(data_role, "training")) "in_sample" else "holdout",
    family = zcov_out$family %||% tryCatch(family(fit)$family, error = function(e) NA_character_),
    type = zcov_out$type %||% type,
    ndraws = M,
    n = N,
    randomized = randomized,
    results = results,
    z_matrix = z_obs_matrix,
    zcov = zcov_out,
    y_rep = y_rep,
    draw_ids = draw_ids
  )

  class(out) <- "z_predcheck"
  out
}

#' Pearson chi-square predictive check
#'
#' @description
#' Computes a moment-based predictive check by comparing the observed Pearson
#' discrepancy with the same discrepancy computed from posterior predictive
#' replicated data.
#'
#' @param fit A fitted model object.
#' @param data Data used to evaluate the fitted model.
#' @param data_role Role of \code{data}: \code{"training"} if the data were used
#'   to fit the model, or \code{"holdout"} otherwise.
#' @param point_Zcov Optional point-level covariate/moment bridge.
#' @param postpred_simulate Optional posterior predictive simulation bridge.
#' @param type Optional component selector.
#' @param ndraws Optional number of posterior draws.
#' @param draw_ids Optional posterior draw indices.
#' @param seed Optional random seed.
#' @param variance_floor Minimum positive variance used for numerical stability.
#' @param ... Additional arguments passed to bridges.
#'
#' @return An object of class \code{"moment_predcheck"}.
#'
#' @export
predcheck_chisquare <- function(fit,
                                data,
                                data_role = c("training", "holdout"),
                                point_Zcov = NULL,
                                postpred_simulate = NULL,
                                type = NULL,
                                ndraws = NULL,
                                draw_ids = NULL,
                                seed = NULL,
                                variance_floor = 1e-6,
                                ...) {
  if (!is.null(seed)) set.seed(seed)

  data_role <- match.arg(data_role)
  eval_data <- as.data.frame(data)

  if (is.null(draw_ids)) {
    draw_ids <- .predcheck_get_draw_ids(fit, ndraws = ndraws)
  }

  zcov_fn <- .predcheck_resolve_point_Zcov(point_Zcov)
  sim_fn <- .predcheck_resolve_postpred_simulate(postpred_simulate)

  zcov_out <- zcov_fn(
    fit = fit,
    data = eval_data,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )

  y_obs <- as.numeric(zcov_out$response)
  if (length(y_obs) < 1L || anyNA(y_obs) || any(!is.finite(y_obs))) {
    stop("predcheck_chisquare: invalid observed response in point_Zcov output.", call. = FALSE)
  }

  y_rep <- sim_fn(
    fit = fit,
    data = eval_data,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )
  y_rep <- .predcheck_as_matrix(y_rep, "y_rep")

  mean_mat <- .predcheck_as_matrix(zcov_out$mean, "mean")
  var_mat <- .predcheck_as_matrix(zcov_out$variance, "variance")

  M <- nrow(y_rep)
  N <- length(y_obs)

  if (ncol(y_rep) != N) {
    stop("predcheck_chisquare: `y_rep` columns do not match observed response length.", call. = FALSE)
  }

  mean_mat <- .predcheck_recycle_matrix_to_draws(mean_mat, M, "mean")
  var_mat <- .predcheck_recycle_matrix_to_draws(var_mat, M, "variance")

  if (ncol(mean_mat) != N || ncol(var_mat) != N) {
    stop("predcheck_chisquare: moment matrix columns do not match observed response length.", call. = FALSE)
  }

  var_mat[!is.finite(var_mat) | var_mat <= 0] <- variance_floor

  y_obs_mat <- matrix(
    rep(y_obs, each = M),
    nrow = M,
    ncol = N,
    byrow = FALSE
  )

  obs_discrepancy <- rowSums((y_obs_mat - mean_mat)^2 / var_mat)
  rep_discrepancy <- rowSums((y_rep - mean_mat)^2 / var_mat)

  out <- list(
    data_role = data_role,
    evaluation = if (identical(data_role, "training")) "in_sample" else "holdout",
    family = zcov_out$family %||% tryCatch(family(fit)$family, error = function(e) NA_character_),
    type = zcov_out$type %||% type,
    ndraws = M,
    n = N,
    p_value = mean(rep_discrepancy >= obs_discrepancy, na.rm = TRUE),
    obs_discrepancy = obs_discrepancy,
    rep_discrepancy = rep_discrepancy,
    zcov = zcov_out,
    draw_ids = draw_ids
  )

  class(out) <- "moment_predcheck"
  out
}

.predictive_check_core <- function(fit,
                                   data,
                                   data_role = c("training", "holdout"),
                                   log_pointpred = NULL,
                                   point_Zcov = NULL,
                                   postpred_simulate = NULL,
                                   test = NULL,
                                   x = NULL,
                                   ndraws = NULL,
                                   seed = NULL,
                                   type = NULL,
                                   randomized = TRUE,
                                   k_anova = 10,
                                   k_bl = 10,
                                   ...) {
  data_role <- match.arg(data_role)
  eval_data <- as.data.frame(data)

  if (!is.null(seed)) set.seed(seed)
  draw_ids <- .predcheck_get_draw_ids(fit, ndraws = ndraws)

  z_rpp <- predcheck_zresid(
    fit = fit,
    data = eval_data,
    data_role = data_role,
    test = test,
    log_pointpred = log_pointpred,
    point_Zcov = point_Zcov,
    postpred_simulate = postpred_simulate,
    type = type,
    x = x,
    ndraws = ndraws,
    draw_ids = draw_ids,
    randomized = randomized,
    k_anova = k_anova,
    k_bl = k_bl,
    ...
  )

  z_mpp <- predcheck_zresid(
    fit = fit,
    data = eval_data,
    data_role = data_role,
    test = list(list(name = "sw")),
    log_pointpred = log_pointpred,
    point_Zcov = point_Zcov,
    postpred_simulate = postpred_simulate,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    randomized = FALSE,
    k_anova = k_anova,
    k_bl = k_bl,
    ...
  )

  chisq <- predcheck_chisquare(
    fit = fit,
    data = eval_data,
    data_role = data_role,
    point_Zcov = point_Zcov,
    postpred_simulate = postpred_simulate,
    type = type,
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )

  get_result_field <- function(obj, name, field) {
    idx <- which(vapply(obj$results, function(r) identical(tolower(r$name), name), logical(1L)))[1L]
    if (is.na(idx)) return(NA_real_)
    obj$results[[idx]][[field]] %||% NA_real_
  }

  get_result_fields <- function(obj, name, field) {
    idx <- which(vapply(obj$results, function(r) identical(tolower(r$name), name), logical(1L)))
    if (length(idx) == 0L) return(NA_real_)

    vals <- vapply(idx, function(i) obj$results[[i]][[field]] %||% NA_real_, numeric(1L))
    labs <- vapply(idx, function(i) obj$results[[i]]$label %||% obj$results[[i]]$name, character(1L))
    names(vals) <- labs
    vals
  }

  get_pred_p <- function(obj, name) get_result_field(obj, name, "p_value")
  get_obs_mean <- function(obj, name) get_result_field(obj, name, "mean_obs_p")

  out <- list(
    data_role = data_role,
    evaluation = if (identical(data_role, "training")) "in_sample" else "holdout",
    family = z_rpp$family,
    type = z_rpp$type,
    ndraws = z_rpp$ndraws,
    n = z_rpp$n,
    chisq = chisq$p_value,
    # User-facing randomized Z-SW summary follows the original PPC function:
    # mean observed randomized Z-residual SW p-value across posterior draws.
    # The posterior predictive comparison p-value is still retained below and
    # inside zresid$results, but is not used as the main printed output.
    sw_rpp = get_obs_mean(z_rpp, "sw"),
    sw_rpp_comparison = get_pred_p(z_rpp, "sw"),
    aov_rpp = get_pred_p(z_rpp, "aov"),
    aov_rpp_all = get_result_fields(z_rpp, "aov", "p_value"),
    bartlett_rpp = get_pred_p(z_rpp, "bartlett"),
    bartlett_rpp_all = get_result_fields(z_rpp, "bartlett", "p_value"),
    sw_mpp = get_pred_p(z_mpp, "sw"),
    mean_sw_pv = get_obs_mean(z_rpp, "sw"),
    mean_aov_pv = get_obs_mean(z_rpp, "aov"),
    mean_aov_pv_all = get_result_fields(z_rpp, "aov", "mean_obs_p"),
    zresid = z_rpp,
    zresid_midpoint = z_mpp,
    chisquare = chisq,
    draw_ids = draw_ids
  )

  class(out) <- "predcheck"
  out
}
