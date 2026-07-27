#' @describeIn Zcov Method for `brmsfit` objects.
#' @export
Zcov.brmsfit <- function(fit,
                         data,
                         type = NULL,
                         point_details = NULL,
                         ndraws = NULL,
                         draw_ids = NULL,
                         ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Zcov.brmsfit requires package 'brms'.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("Zcov.brmsfit: `data` must be provided.", call. = FALSE)
  }

  point_details <- .Zcov_normalize_point_details(point_details)

  data_in <- data
  fam <- fit$family$family
  is_hurdle <- grepl("^hurdle(_|$)", fam)
  is_zi <- grepl("^zero_inflated(_|$)", fam)
  is_trunc <- grepl("^truncated(_|$)", fam)
  has_component <- is_hurdle || is_zi

  type_req <- .Zcov_type_normalize(type, is_hurdle = has_component)

  if (!has_component) {
    if (!is.null(type_req)) {
      warning(
        sprintf(
          "Zcov.brmsfit: ignoring type='%s' for non-hurdle/non-zero-inflated family='%s'.",
          type_req,
          fam
        ),
        call. = FALSE
      )
    }
    type_req <- NULL
  } else if (is.null(type_req) || !(type_req %in% c("hurdle", "count", "zero"))) {
    type_req <- "hurdle"
  }

  mf_all <- stats::model.frame(fit$formula, data = data_in, na.action = stats::na.pass)
  response_name <- .Zcov_resp_name_from_brms(fit, mf_all)

  keep <- stats::complete.cases(mf_all)
  .Zcov_warn_drop_na("Zcov.brmsfit", nrow(mf_all), keep, mf_all)

  obs_id <- which(keep)
  mf <- mf_all[keep, , drop = FALSE]
  data_used <- data_in[keep, , drop = FALSE]

  y_raw <- mf[[response_name]]
  y_num <- suppressWarnings(as.numeric(y_raw))

  cov_all <- data_in[, setdiff(names(data_in), response_name), drop = FALSE]
  cov_used <- cov_all[obs_id, , drop = FALSE]

  zero_id <- which(!is.na(y_num) & y_num == 0)
  zero_obs_id <- obs_id[zero_id]

  if (has_component) {
    y_type_kind <- if (is_hurdle) "hurdle" else "zero_inflated"
    y_type_levels <- c(zero = 0L, count = 1L)
    y_type <- ifelse(!is.na(y_num) & y_num == 0, 0L, 1L)
  } else {
    y_type_kind <- "plain"
    y_type_levels <- c(obs = 1L)
    y_type <- rep.int(1L, length(obs_id))
  }

  lp_cond <- .Zcov_brms_linear_pred_old(
    fit = fit,
    data = data_used,
    type_req = type_req,
    fam = fam
  )

  extra <- list(
    family = fam,
    zero_id = zero_id,
    zero_obs_id = zero_obs_id,
    type_raw = type,
    type_resolved = type_req
  )

  if (is_trunc) {
    keep2 <- which(!is.na(y_num) & y_num > 0)
    out <- list(
      type = "count",
      family = fam,
      response_name = response_name,
      response = y_raw[keep2],
      covariates = cov_used[keep2, , drop = FALSE],
      linear_pred = lp_cond[keep2],
      obs_id = obs_id[keep2],
      y_type = rep.int(1L, length(keep2)),
      y_type_kind = "trunc",
      y_type_levels = c(count = 1L),
      extra = extra
    )

    return(.Zcov_brms_add_point_details(
      out = out,
      fit = fit,
      data = data_used[keep2, , drop = FALSE],
      fam = fam,
      type_req = "count",
      point_details = point_details,
      ndraws = ndraws,
      draw_ids = draw_ids
    ))
  }

  if (has_component && identical(type_req, "count")) {
    keep_cnt <- which(!is.na(y_num) & y_num > 0)
    out <- list(
      type = "count",
      family = fam,
      response_name = response_name,
      response = y_raw[keep_cnt],
      covariates = cov_used[keep_cnt, , drop = FALSE],
      linear_pred = lp_cond[keep_cnt],
      obs_id = obs_id[keep_cnt],
      y_type = rep.int(1L, length(keep_cnt)),
      y_type_kind = "trunc",
      y_type_levels = c(count = 1L),
      extra = extra
    )

    return(.Zcov_brms_add_point_details(
      out = out,
      fit = fit,
      data = data_used[keep_cnt, , drop = FALSE],
      fam = fam,
      type_req = "count",
      point_details = point_details,
      ndraws = ndraws,
      draw_ids = draw_ids
    ))
  }

  if (has_component && identical(type_req, "zero")) {
    o_i <- as.integer(!is.na(y_num) & y_num == 0)

    out <- list(
      type = "zero",
      family = fam,
      response_name = response_name,
      response = o_i,
      covariates = cov_used,
      linear_pred = lp_cond,
      obs_id = obs_id,
      y_type = y_type,
      y_type_kind = y_type_kind,
      y_type_levels = y_type_levels,
      extra = extra
    )

    return(.Zcov_brms_add_point_details(
      out = out,
      fit = fit,
      data = data_used,
      fam = fam,
      type_req = "zero",
      point_details = point_details,
      ndraws = ndraws,
      draw_ids = draw_ids
    ))
  }

  out <- list(
    type = if (has_component) "hurdle" else NULL,
    family = fam,
    response_name = response_name,
    response = y_raw,
    covariates = cov_used,
    linear_pred = lp_cond,
    obs_id = obs_id,
    y_type = y_type,
    y_type_kind = y_type_kind,
    y_type_levels = y_type_levels,
    extra = extra
  )

  .Zcov_brms_add_point_details(
    out = out,
    fit = fit,
    data = data_used,
    fam = fam,
    type_req = type_req,
    point_details = point_details,
    ndraws = ndraws,
    draw_ids = draw_ids
  )
}

.Zcov_resp_name_from_brms <- function(fit, mf_all) {
  if (is.null(mf_all) || !is.data.frame(mf_all) || ncol(mf_all) < 1L) {
    stop("Zcov.brmsfit: unable to determine response from model frame.", call. = FALSE)
  }

  names(mf_all)[1L]
}

.Zcov_brms_linear_pred_old <- function(fit, data, type_req = NULL, fam = NULL) {
  dpar <- NULL

  if (identical(type_req, "zero")) {
    dpar <- .Zcov_brms_zero_dpar(fam)
  } else if (identical(type_req, "count")) {
    dpar <- "mu"
  }

  lp_mat <- tryCatch(
    brms::posterior_linpred(
      fit,
      newdata = data,
      dpar = dpar,
      transform = FALSE,
      re_formula = NULL
    ),
    error = function(e) NULL
  )

  if (is.null(lp_mat)) {
    lp_mat <- tryCatch(
      brms::posterior_linpred(
        fit,
        newdata = data,
        transform = FALSE,
        re_formula = NULL
      ),
      error = function(e) NULL
    )
  }

  if (is.null(lp_mat)) {
    stop(
      "Zcov.brmsfit: failed to compute linear_pred via brms::posterior_linpred(transform = FALSE).",
      call. = FALSE
    )
  }

  lp_mat <- as.matrix(lp_mat)

  if (ncol(lp_mat) != nrow(data)) {
    stop(
      "Zcov.brmsfit: linear_pred dimension mismatch. ncol(linear_pred) must equal nrow(data).",
      call. = FALSE
    )
  }

  as.numeric(colMeans(lp_mat))
}

.Zcov_brms_add_point_details <- function(out,
                                         fit,
                                         data,
                                         fam,
                                         type_req = NULL,
                                         point_details = NULL,
                                         ndraws = NULL,
                                         draw_ids = NULL) {
  if (is.null(point_details)) return(out)

  draw_ids <- .Zcov_brms_draw_ids(fit, ndraws = ndraws, draw_ids = draw_ids)

  if ("mean" %in% point_details || "variance" %in% point_details) {
    out$mean <- .Zcov_brms_point_mean(
      fit = fit,
      data = data,
      fam = fam,
      type_req = type_req,
      draw_ids = draw_ids
    )
  }

  if ("lp" %in% point_details) {
    out$lp <- .Zcov_brms_point_lp(
      fit = fit,
      data = data,
      type_req = type_req,
      draw_ids = draw_ids
    )
  }

  if ("variance" %in% point_details) {
    out$variance <- .Zcov_brms_point_variance(
      fit = fit,
      data = data,
      fam = fam,
      type_req = type_req,
      draw_ids = draw_ids,
      mean_mat = out$mean
    )
  }

  if ("covariate" %in% point_details || length(point_details) > 0L) {
    out$covariate <- .Zcov_safe_covariate_matrix(
      covariates = out$covariates,
      where = "Zcov.brmsfit"
    )
  }

  out
}

.Zcov_brms_draw_ids <- function(fit, ndraws = NULL, draw_ids = NULL) {
  draws_total <- nrow(as.matrix(fit))
  if (draws_total < 1L) {
    stop("Zcov.brmsfit: no posterior draws found in brmsfit object.", call. = FALSE)
  }

  if (!is.null(draw_ids)) {
    draw_ids <- as.integer(draw_ids)
    if (anyNA(draw_ids) || any(draw_ids < 1L) || any(draw_ids > draws_total)) {
      stop("Zcov.brmsfit: `draw_ids` contains invalid posterior draw indices.", call. = FALSE)
    }
    return(draw_ids)
  }

  if (is.null(ndraws)) return(NULL)

  ndraws <- as.integer(ndraws[1L])
  if (is.na(ndraws) || ndraws < 1L) {
    stop("Zcov.brmsfit: ndraws must be a positive integer.", call. = FALSE)
  }

  seq_len(min(ndraws, draws_total))
}

.Zcov_brms_point_mean <- function(fit, data, fam, type_req = NULL, draw_ids = NULL) {
  if (identical(type_req, "zero")) {
    return(.Zcov_brms_dpar_matrix(fit, data, dpar = .Zcov_brms_zero_dpar(fam), draw_ids = draw_ids))
  }

  if (identical(type_req, "count") && grepl("^hurdle", fam)) {
    mu <- .Zcov_brms_dpar_matrix(fit, data, dpar = "mu", draw_ids = draw_ids)
    return(.Zcov_brms_hurdle_count_mean(fit, data, fam, mu, draw_ids))
  }

  out <- tryCatch(
    brms::posterior_epred(fit, newdata = data, draw_ids = draw_ids, re_formula = NULL),
    error = function(e) NULL
  )

  if (is.null(out)) {
    stop("Zcov.brmsfit: failed to compute point_details mean via posterior_epred().", call. = FALSE)
  }

  as.matrix(out)
}

.Zcov_brms_point_lp <- function(fit, data, type_req = NULL, draw_ids = NULL) {
  dpar <- NULL
  transform <- FALSE

  if (identical(type_req, "zero")) {
    dpar <- "hu"
  } else if (identical(type_req, "count")) {
    dpar <- "mu"
  }

  out <- tryCatch(
    brms::posterior_linpred(
      fit,
      newdata = data,
      draw_ids = draw_ids,
      re_formula = NULL,
      dpar = dpar,
      transform = transform
    ),
    error = function(e) NULL
  )

  if (is.null(out)) {
    out <- tryCatch(
      brms::posterior_linpred(
        fit,
        newdata = data,
        draw_ids = draw_ids,
        re_formula = NULL,
        transform = transform
      ),
      error = function(e) NULL
    )
  }

  if (is.null(out)) {
    stop("Zcov.brmsfit: failed to compute point_details lp via posterior_linpred().", call. = FALSE)
  }

  as.matrix(out)
}

.Zcov_brms_point_variance <- function(fit,
                                      data,
                                      fam,
                                      type_req = NULL,
                                      draw_ids = NULL,
                                      mean_mat = NULL) {
  if (is.null(mean_mat)) {
    mean_mat <- .Zcov_brms_point_mean(fit, data, fam, type_req, draw_ids)
  }

  if (identical(type_req, "zero")) {
    return(pmax(mean_mat * (1 - mean_mat), .Machine$double.eps))
  }

  if (grepl("^hurdle", fam)) {
    return(.Zcov_brms_hurdle_variance(fit, data, fam, type_req, draw_ids, mean_mat))
  }

  if (grepl("^zero_inflated", fam)) {
    warning(
      sprintf("Zcov.brmsfit: variance for family='%s' is not implemented yet; returning NA matrix.", fam),
      call. = FALSE
    )
    return(matrix(NA_real_, nrow = nrow(mean_mat), ncol = ncol(mean_mat)))
  }

  .Zcov_brms_family_variance(fit, data, fam, mean_mat, draw_ids)
}

.Zcov_brms_family_variance <- function(fit, data, fam, mu, draw_ids = NULL) {
  if (identical(fam, "gaussian")) {
    sigma <- .Zcov_brms_dpar_matrix(fit, data, dpar = "sigma", draw_ids = draw_ids)
    return(pmax(sigma^2, .Machine$double.eps))
  }

  if (identical(fam, "poisson")) {
    return(pmax(mu, .Machine$double.eps))
  }

  if (identical(fam, "negbinomial")) {
    shape <- .Zcov_brms_dpar_matrix(fit, data, dpar = "shape", draw_ids = draw_ids)
    return(pmax(mu + mu^2 / shape, .Machine$double.eps))
  }

  if (identical(fam, "bernoulli")) {
    return(pmax(mu * (1 - mu), .Machine$double.eps))
  }

  if (identical(fam, "Gamma")) {
    shape <- .Zcov_brms_dpar_matrix(fit, data, dpar = "shape", draw_ids = draw_ids)
    return(pmax(mu^2 / shape, .Machine$double.eps))
  }

  if (identical(fam, "inverse.gaussian")) {
    shape <- .Zcov_brms_dpar_matrix(fit, data, dpar = "shape", draw_ids = draw_ids)
    return(pmax(mu^3 / shape, .Machine$double.eps))
  }

  if (identical(fam, "beta")) {
    phi <- .Zcov_brms_dpar_matrix(fit, data, dpar = "phi", draw_ids = draw_ids)
    return(pmax(mu * (1 - mu) / (phi + 1), .Machine$double.eps))
  }

  warning(
    sprintf("Zcov.brmsfit: variance is not implemented for brms family='%s'; returning NA matrix.", fam),
    call. = FALSE
  )
  matrix(NA_real_, nrow = nrow(mu), ncol = ncol(mu))
}

.Zcov_brms_dpar_matrix <- function(fit, data, dpar, draw_ids = NULL) {
  out <- tryCatch(
    brms::posterior_linpred(
      fit,
      newdata = data,
      draw_ids = draw_ids,
      re_formula = NULL,
      dpar = dpar,
      transform = TRUE
    ),
    error = function(e) NULL
  )

  if (!is.null(out)) return(as.matrix(out))

  par_draws <- .Zcov_brms_scalar_parameter(fit, dpar, draw_ids)
  if (!is.null(par_draws)) {
    return(matrix(
      rep(par_draws, times = nrow(data)),
      nrow = length(par_draws),
      ncol = nrow(data),
      byrow = FALSE
    ))
  }

  stop(sprintf("Zcov.brmsfit: failed to extract distributional parameter '%s'.", dpar), call. = FALSE)
}

.Zcov_brms_scalar_parameter <- function(fit, variable, draw_ids = NULL) {
  mat <- tryCatch(as.matrix(fit, variable = variable), error = function(e) NULL)
  if (is.null(mat) || ncol(mat) < 1L) return(NULL)

  vals <- as.numeric(mat[, 1L])
  if (!is.null(draw_ids)) vals <- vals[draw_ids]
  vals
}

.Zcov_brms_zero_dpar <- function(fam) {
  if (grepl("^zero_inflated", fam)) return("zi")
  "hu"
}

.Zcov_brms_hurdle_count_mean <- function(fit, data, fam, mu, draw_ids = NULL) {
  f0 <- .Zcov_brms_base_zero_prob(fit, data, fam, mu, draw_ids)
  mu / pmax(1 - f0, .Machine$double.eps)
}

.Zcov_brms_hurdle_variance <- function(fit, data, fam, type_req, draw_ids = NULL, mean_mat = NULL) {
  mu <- .Zcov_brms_dpar_matrix(fit, data, dpar = "mu", draw_ids = draw_ids)
  base_var <- .Zcov_brms_base_count_variance(fit, data, fam, mu, draw_ids)
  f0 <- .Zcov_brms_base_zero_prob(fit, data, fam, mu, draw_ids)

  e2_pos <- (base_var + mu^2) / pmax(1 - f0, .Machine$double.eps)
  mean_pos <- mu / pmax(1 - f0, .Machine$double.eps)

  if (identical(type_req, "count")) {
    return(pmax(e2_pos - mean_pos^2, .Machine$double.eps))
  }

  hu <- .Zcov_brms_dpar_matrix(fit, data, dpar = "hu", draw_ids = draw_ids)
  mean_all <- if (is.null(mean_mat)) (1 - hu) * mean_pos else mean_mat
  e2_all <- (1 - hu) * e2_pos

  pmax(e2_all - mean_all^2, .Machine$double.eps)
}

.Zcov_brms_base_count_variance <- function(fit, data, fam, mu, draw_ids = NULL) {
  if (grepl("poisson", fam)) {
    return(pmax(mu, .Machine$double.eps))
  }

  if (grepl("negbinomial", fam)) {
    shape <- .Zcov_brms_dpar_matrix(fit, data, dpar = "shape", draw_ids = draw_ids)
    return(pmax(mu + mu^2 / shape, .Machine$double.eps))
  }

  warning(
    sprintf("Zcov.brmsfit: hurdle count variance is not implemented for family='%s'; returning NA matrix.", fam),
    call. = FALSE
  )
  matrix(NA_real_, nrow = nrow(mu), ncol = ncol(mu))
}

.Zcov_brms_base_zero_prob <- function(fit, data, fam, mu, draw_ids = NULL) {
  if (grepl("poisson", fam)) {
    return(exp(-mu))
  }

  if (grepl("negbinomial", fam)) {
    shape <- .Zcov_brms_dpar_matrix(fit, data, dpar = "shape", draw_ids = draw_ids)
    return((shape / (shape + mu))^shape)
  }

  matrix(0, nrow = nrow(mu), ncol = ncol(mu))
}
