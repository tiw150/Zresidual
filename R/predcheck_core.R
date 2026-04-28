# =========================================================
# Core infrastructure for predictive checks
# =========================================================

get_predcheck <- function(fit,
                          data,
                          predcheck_pointpred = NULL,
                          draw_ids = NULL,
                          ...) {
  if (!is.null(predcheck_pointpred)) {
    fun <- match.fun(predcheck_pointpred)
    obj <- fun(fit = fit, data = data, draw_ids = draw_ids, ...)
    .validate_predcheck_object(obj)
    return(obj)
  }

  backend_name <- .infer_predcheck_backend(fit)

  if (!exists(backend_name, mode = "function", inherits = TRUE)) {
    stop(
      "Backend function not found: ", backend_name,
      ". Please define it or pass `predcheck_pointpred` explicitly."
    )
  }

  fun <- get(backend_name, mode = "function", inherits = TRUE)
  obj <- fun(fit = fit, data = data, draw_ids = draw_ids, ...)
  .validate_predcheck_object(obj)
  obj
}

.predcheck_compute_z <- function(y, pred, draw, randomized = TRUE, eps = .Machine$double.eps) {
  pmf_vec <- pred$pmf(y, draw = draw)
  tail_vec <- pred$tail(y, draw = draw)

  u <- if (randomized) stats::runif(pred$n) else rep(0.5, pred$n)
  pval <- .predcheck_clip_prob(tail_vec + u * pmf_vec, eps = eps)

  -stats::qnorm(pval)
}

.predictive_check_core <- function(fit,
                                   data,
                                   mode = c("ppc", "hpc"),
                                   predcheck_pointpred = NULL,
                                   x = NULL,
                                   ndraws = NULL,
                                   seed = NULL,
                                   k_anova = 10,
                                   ...) {
  mode <- match.arg(mode)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  draw_ids <- .predcheck_get_draw_ids(fit, ndraws = ndraws)

  pred <- get_predcheck(
    fit = fit,
    data = data,
    predcheck_pointpred = predcheck_pointpred,
    draw_ids = draw_ids,
    ...
  )

  y <- pred$y
  y_rep <- pred$rng()

  if (!is.matrix(y_rep)) {
    y_rep <- as.matrix(y_rep)
  }

  if (nrow(y_rep) != pred$ndraws || ncol(y_rep) != pred$n) {
    stop("`pred$rng()` returned a matrix with inconsistent dimensions.")
  }

  include_aov <- !is.null(x)
  if (include_aov) {
    x_vec <- .predcheck_extract_x(data = data, x = x, n_expected = pred$n)
  }

  chisq_comp <- logical(pred$ndraws)
  sw_rpp_comp <- logical(pred$ndraws)
  sw_mpp_comp <- logical(pred$ndraws)
  pv_sw_obs <- numeric(pred$ndraws)

  if (include_aov) {
    aov_rpp_comp <- logical(pred$ndraws)
    pv_aov_obs <- numeric(pred$ndraws)
  } else {
    aov_rpp_comp <- rep(NA, pred$ndraws)
    pv_aov_obs <- rep(NA_real_, pred$ndraws)
  }

  for (i in seq_len(pred$ndraws)) {
    moment_i <- pred$moments(draw = i)
    mean_i <- moment_i$mean
    var_i <- moment_i$var
    var_i[!is.finite(var_i) | var_i <= 0] <- 1e-6

    stat_rep <- sum((y_rep[i, ] - mean_i)^2 / var_i)
    stat_obs <- sum((y        - mean_i)^2 / var_i)
    chisq_comp[i] <- (stat_obs <= stat_rep)

    z_rep_rpp <- .predcheck_compute_z(y = y_rep[i, ], pred = pred, draw = i, randomized = TRUE)
    z_obs_rpp <- .predcheck_compute_z(y = y,         pred = pred, draw = i, randomized = TRUE)

    sw_rep <- sw.test.zresid(z_rep_rpp)
    sw_obs <- sw.test.zresid(z_obs_rpp)

    sw_rpp_comp[i] <- (sw_rep <= sw_obs)
    pv_sw_obs[i] <- sw_obs

    if (include_aov) {
      aov_rep <- test.nl.aov(z_rep_rpp, x_vec, k.anova = k_anova)
      aov_obs <- test.nl.aov(z_obs_rpp, x_vec, k.anova = k_anova)

      aov_rpp_comp[i] <- (aov_rep <= aov_obs)
      pv_aov_obs[i] <- aov_obs
    }

    z_rep_mpp <- .predcheck_compute_z(y = y_rep[i, ], pred = pred, draw = i, randomized = FALSE)
    z_obs_mpp <- .predcheck_compute_z(y = y,         pred = pred, draw = i, randomized = FALSE)

    sw_rep_m <- sw.test.zresid(z_rep_mpp)
    sw_obs_m <- sw.test.zresid(z_obs_mpp)

    sw_mpp_comp[i] <- (sw_rep_m <= sw_obs_m)
  }

  out <- list(
    mode = mode,
    family = pred$family,
    ndraws = pred$ndraws,
    n = pred$n,
    chisq = mean(chisq_comp),
    sw_rpp = mean(sw_rpp_comp),
    aov_rpp = if (include_aov) mean(aov_rpp_comp) else NA_real_,
    sw_mpp = mean(sw_mpp_comp),
    mean_sw_pv = mean(pv_sw_obs),
    mean_aov_pv = if (include_aov) mean(pv_aov_obs) else NA_real_
  )

  class(out) <- "predcheck"
  out
}
