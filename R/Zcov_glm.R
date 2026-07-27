#' @describeIn Zcov Method for `glm` objects.
#' @export
Zcov.glm <- function(fit,
                     data,
                     type = NULL,
                     point_details = NULL,
                     ndraws = NULL,
                     draw_ids = NULL,
                     ...) {
  if (missing(data) || is.null(data)) {
    stop("Zcov.glm: `data` must be provided.", call. = FALSE)
  }

  point_details <- .Zcov_normalize_point_details(point_details)
  fam <- fit$family$family

  if (!is.null(type)) {
    warning(
      sprintf(
        "Zcov.glm: ignoring type='%s' for standard glm family='%s'.",
        type,
        fam
      ),
      call. = FALSE
    )
  }

  mf_all <- stats::model.frame(stats::formula(fit), data = data, na.action = stats::na.pass)
  response_name <- names(mf_all)[1L]

  keep <- stats::complete.cases(mf_all)
  .Zcov_warn_drop_na("Zcov.glm", nrow(mf_all), keep, mf_all)

  obs_id <- which(keep)
  mf <- mf_all[keep, , drop = FALSE]
  data_used <- data[keep, , drop = FALSE]

  y_raw <- stats::model.response(mf)

  cov_all <- data[, setdiff(names(data), response_name), drop = FALSE]
  cov_used <- cov_all[obs_id, , drop = FALSE]

  lp_cond <- tryCatch(
    as.numeric(stats::predict(fit, newdata = data_used, type = "link")),
    error = function(e) {
      stop("Zcov.glm: failed to compute linear_pred via predict(type='link').", call. = FALSE)
    }
  )

  out <- list(
    type = NULL,
    family = fam,
    response_name = response_name,
    response = y_raw,
    covariates = cov_used,
    linear_pred = lp_cond,
    obs_id = obs_id,
    y_type = rep.int(1L, length(obs_id)),
    y_type_kind = "plain",
    y_type_levels = c(obs = 1L),
    extra = list(
      family = fam,
      type_raw = type
    )
  )

  .Zcov_glm_add_point_details(
    out = out,
    fit = fit,
    data = data_used,
    point_details = point_details
  )
}

.Zcov_glm_add_point_details <- function(out,
                                        fit,
                                        data,
                                        point_details = NULL) {
  if (is.null(point_details)) return(out)

  fam <- fit$family$family

  if ("lp" %in% point_details) {
    out$lp <- .Zcov_matrix_1xn(stats::predict(fit, newdata = data, type = "link"))
  }

  if ("mean" %in% point_details || "variance" %in% point_details) {
    mu <- as.numeric(stats::predict(fit, newdata = data, type = "response"))
    out$mean <- .Zcov_matrix_1xn(mu)
  }

  if ("variance" %in% point_details) {
    out$variance <- .Zcov_matrix_1xn(.Zcov_glm_variance(fit, mu, data))
  }

  if ("covariate" %in% point_details || length(point_details) > 0L) {
    out$covariate <- .Zcov_safe_model_matrix(stats::delete.response(stats::terms(fit)), data, "Zcov.glm")
  }

  out
}

.Zcov_glm_variance <- function(fit, mu, data) {
  fam <- fit$family$family
  dispersion <- tryCatch(summary(fit)$dispersion, error = function(e) 1)
  if (!is.finite(dispersion) || length(dispersion) != 1L) dispersion <- 1

  switch(
    fam,
    gaussian = rep(dispersion, length(mu)),
    poisson = mu,
    quasipoisson = dispersion * mu,
    binomial = pmax(mu * (1 - mu), .Machine$double.eps),
    quasibinomial = dispersion * pmax(mu * (1 - mu), .Machine$double.eps),
    Gamma = dispersion * mu^2,
    inverse.gaussian = dispersion * mu^3,
    {
      warning(
        sprintf("Zcov.glm: variance is not implemented for family='%s'.", fam),
        call. = FALSE
      )
      rep(NA_real_, length(mu))
    }
  )
}
