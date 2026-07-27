#' @describeIn Zcov Method for `coxph` objects from the survival package.
#' @export
Zcov.coxph <- function(fit,
                       data,
                       type = NULL,
                       point_details = NULL,
                       ndraws = NULL,
                       ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.coxph requires package 'survival'.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("Zcov.coxph: `data` must be provided.", call. = FALSE)
  }

  point_details <- .Zcov_normalize_point_details(point_details)

  mf_all <- stats::model.frame(fit$terms, data, na.action = stats::na.pass)
  keep <- stats::complete.cases(mf_all)
  .Zcov_warn_drop_na("Zcov.coxph", nrow(mf_all), keep, mf_all)

  data_used <- data[keep, , drop = FALSE]
  mf_new <- stats::model.frame(fit$terms, data_used)
  y_new <- mf_new[[1L]]

  lp <- tryCatch(
    as.vector(stats::predict(fit, newdata = data_used, type = "lp")),
    error = function(e) {
      mm_new <- model.matrix.coxph(fit, data_used)
      as.vector(mm_new %*% fit$coefficients)
    }
  )

  obs_id <- which(keep)

  if (!inherits(y_new, "Surv")) {
    stop("Zcov.coxph: response is not a Surv object.", call. = FALSE)
  }

  y_mat <- as.matrix(y_new)
  if (ncol(y_mat) == 2L) {
    time <- y_mat[, 1L]
    status <- y_mat[, 2L]
  } else if (ncol(y_mat) == 3L) {
    time <- y_mat[, 2L]
    status <- y_mat[, 3L]
  } else {
    stop("Zcov.coxph: unsupported Surv format.", call. = FALSE)
  }

  out <- list(
    type = NULL,
    family = "coxph",
    response_name = .Zcov_surv_response_name(fit, mf_new),
    response = y_new,
    covariates = mf_new[, -1L, drop = FALSE],
    linear_pred = lp,
    obs_id = obs_id,
    y_type = ifelse(status == 1, 1L, 0L),
    y_type_kind = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra = c(
      list(
        time = time,
        status = status,
        event_id = which(status == 1),
        censor_id = which(status == 0),
        type_raw = type
      ),
      .Zcov_surv_arg_names(fit)
    )
  )

  .Zcov_coxph_add_point_details(
    out = out,
    fit = fit,
    data = data_used,
    point_details = point_details
  )
}

.Zcov_coxph_add_point_details <- function(out,
                                          fit,
                                          data,
                                          point_details = NULL) {
  if (is.null(point_details)) return(out)

  if ("lp" %in% point_details) {
    out$lp <- .Zcov_matrix_1xn(out$linear_pred)
  }

  if ("mean" %in% point_details || "variance" %in% point_details) {
    warning(
      "Zcov.coxph: point_details 'mean' and 'variance' are not defined for semi-parametric Cox models. Only lp and covariate are returned.",
      call. = FALSE
    )
  }

  if ("covariate" %in% point_details || length(point_details) > 0L) {
    out$covariate <- .Zcov_safe_model_matrix(fit$terms, data, "Zcov.coxph")
    if (ncol(out$covariate) > 0L && colnames(out$covariate)[1L] == "(Intercept)") {
      out$covariate <- out$covariate[, -1L, drop = FALSE]
    }
  }

  out
}

#' @describeIn Zcov Method for `survreg` objects from the survival package.
#' @export
Zcov.survreg <- function(fit,
                         data,
                         type = NULL,
                         point_details = NULL,
                         ndraws = NULL,
                         ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Zcov.survreg requires package 'survival'.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("Zcov.survreg: `data` must be provided.", call. = FALSE)
  }

  point_details <- .Zcov_normalize_point_details(point_details)

  mf_all <- stats::model.frame(fit$terms, data, na.action = stats::na.pass)
  keep <- stats::complete.cases(mf_all)
  .Zcov_warn_drop_na("Zcov.survreg", nrow(mf_all), keep, mf_all)

  data_used <- data[keep, , drop = FALSE]
  mf <- stats::model.frame(fit$terms, data_used)
  y <- mf[[1L]]

  lp <- tryCatch(
    as.vector(stats::predict(fit, newdata = data_used, type = "lp")),
    error = function(e) {
      mm <- stats::model.matrix(fit$terms, data_used)
      as.vector(mm %*% fit$coefficients)
    }
  )

  obs_id <- which(keep)

  if (!inherits(y, "Surv")) {
    stop("Zcov.survreg: response is not a Surv object.", call. = FALSE)
  }

  y_mat <- as.matrix(y)
  if (ncol(y_mat) < 2L) {
    stop("Zcov.survreg: Surv object must have at least time and status columns.", call. = FALSE)
  }

  time <- y_mat[, 1L]
  status <- y_mat[, 2L]

  out <- list(
    type = NULL,
    family = "survreg",
    response_name = .Zcov_surv_response_name(fit, mf),
    response = y,
    covariates = mf[, -1L, drop = FALSE],
    linear_pred = lp,
    obs_id = obs_id,
    y_type = ifelse(status == 1, 1L, 0L),
    y_type_kind = "censor",
    y_type_levels = c(censored = 0L, event = 1L),
    extra = c(
      list(
        time = time,
        status = status,
        event_id = which(status == 1),
        censor_id = which(status == 0),
        dist = fit$dist,
        scale = fit$scale,
        type_raw = type
      ),
      .Zcov_surv_arg_names(fit)
    )
  )

  .Zcov_survreg_add_point_details(
    out = out,
    fit = fit,
    data = data_used,
    point_details = point_details
  )
}

.Zcov_survreg_add_point_details <- function(out,
                                            fit,
                                            data,
                                            point_details = NULL) {
  if (is.null(point_details)) return(out)

  eta <- out$linear_pred

  if ("lp" %in% point_details) {
    out$lp <- .Zcov_matrix_1xn(eta)
  }

  if ("mean" %in% point_details || "variance" %in% point_details) {
    moments <- .Zcov_survreg_moments(fit, eta)
    if ("mean" %in% point_details) out$mean <- .Zcov_matrix_1xn(moments$mean)
    if ("variance" %in% point_details) out$variance <- .Zcov_matrix_1xn(moments$variance)
  }

  if ("covariate" %in% point_details || length(point_details) > 0L) {
    out$covariate <- .Zcov_safe_model_matrix(fit$terms, data, "Zcov.survreg")
    if (ncol(out$covariate) > 0L && colnames(out$covariate)[1L] == "(Intercept)") {
      out$covariate <- out$covariate[, -1L, drop = FALSE]
    }
  }

  out
}

.Zcov_survreg_moments <- function(fit, eta) {
  dist <- fit$dist
  scale <- fit$scale

  if (identical(dist, "weibull")) {
    mean <- exp(eta) * gamma(1 + scale)
    variance <- exp(2 * eta) * (gamma(1 + 2 * scale) - gamma(1 + scale)^2)
  } else if (identical(dist, "exponential")) {
    mean <- exp(eta)
    variance <- exp(2 * eta)
  } else if (identical(dist, "lognormal")) {
    mean <- exp(eta + 0.5 * scale^2)
    variance <- (exp(scale^2) - 1) * exp(2 * eta + scale^2)
  } else if (identical(dist, "loglogistic")) {
    mean <- rep(NA_real_, length(eta))
    variance <- rep(NA_real_, length(eta))

    if (scale < 1) {
      mean <- exp(eta) * pi * scale / sin(pi * scale)
    }
    if (scale < 0.5) {
      e2 <- exp(2 * eta) * 2 * pi * scale / sin(2 * pi * scale)
      variance <- e2 - mean^2
    }
  } else if (identical(dist, "gaussian")) {
    mean <- eta
    variance <- rep(scale^2, length(eta))
  } else {
    warning(
      sprintf("Zcov.survreg: moments are not implemented for dist='%s'.", dist),
      call. = FALSE
    )
    mean <- rep(NA_real_, length(eta))
    variance <- rep(NA_real_, length(eta))
  }

  list(mean = mean, variance = pmax(variance, .Machine$double.eps))
}
