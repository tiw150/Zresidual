#' Simulate posterior predictive draws for predictive checks
#'
#' @description
#' Generic bridge used by \code{predcheck_chisquare}. Methods should return an
#' M x N matrix, where rows are posterior draws or replications and columns are
#' observations.
#'
#' @param fit A fitted model object.
#' @param data Evaluation data.
#' @param type Optional component selector.
#' @param ndraws Optional number of draws.
#' @param draw_ids Optional posterior draw indices.
#' @param ... Additional arguments passed to methods.
#'
#' @export
postpred_simulate <- function(fit,
                              data,
                              type = NULL,
                              ndraws = NULL,
                              draw_ids = NULL,
                              ...) {
  UseMethod("postpred_simulate")
}

#' @method postpred_simulate brmsfit
#' @export
postpred_simulate.brmsfit <- function(fit,
                                      data,
                                      type = NULL,
                                      ndraws = NULL,
                                      draw_ids = NULL,
                                      ...) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("postpred_simulate.brmsfit requires package 'brms'.", call. = FALSE)
  }

  zcov_out <- Zcov(fit, data = data, type = type)
  data_used <- as.data.frame(data)[zcov_out$obs_id, , drop = FALSE]
  fam <- fit$family$family
  type_req <- zcov_out$type

  if (is.null(draw_ids)) {
    draw_ids <- .predcheck_get_draw_ids(fit, ndraws = ndraws)
  }

  if (identical(type_req, "zero")) {
    hu <- posterior.pred(
      fit = fit,
      dpar = "zero",
      data = data_used,
      count.only = FALSE,
      draw_ids = draw_ids,
      ...
    )
    return(matrix(stats::rbinom(length(hu), size = 1L, prob = as.vector(hu)),
                  nrow = nrow(hu), ncol = ncol(hu)))
  }

  if (identical(type_req, "count") && grepl("^hurdle", fam)) {
    mu <- posterior.pred(
      fit = fit,
      dpar = "mu",
      data = data_used,
      count.only = FALSE,
      draw_ids = draw_ids,
      ...
    )

    if (grepl("poisson", fam)) {
      return(.postpred_rtruncpois_matrix(mu))
    }

    if (grepl("negbinomial", fam)) {
      shape <- posterior.pred(
        fit = fit,
        dpar = "shape",
        data = data_used,
        count.only = FALSE,
        draw_ids = draw_ids,
        ...
      )
      return(.postpred_rtruncnb_matrix(mu = mu, size = shape))
    }
  }

  if (identical(fam, "inverse.gaussian")) {
    mu <- posterior.pred(
      fit = fit,
      dpar = "mu",
      data = data_used,
      count.only = FALSE,
      draw_ids = draw_ids,
      ...
    )

    shape <- posterior.pred(
      fit = fit,
      dpar = "shape",
      data = data_used,
      count.only = FALSE,
      draw_ids = draw_ids,
      ...
    )

    return(rinvgauss(mu = mu, shape = shape))
  }

  args <- list(fit, newdata = data_used, draw_ids = draw_ids)
  dots <- list(...)
  if (length(dots) > 0L) args <- c(args, dots)

  out <- do.call(brms::posterior_predict, args)
  as.matrix(out)
}

#' @method postpred_simulate glm
#' @export
postpred_simulate.glm <- function(fit,
                                  data,
                                  type = NULL,
                                  ndraws = NULL,
                                  draw_ids = NULL,
                                  ...) {
  if (!is.null(type)) {
    warning("postpred_simulate.glm: `type` is ignored for standard glm objects.", call. = FALSE)
  }

  zcov_out <- Zcov(fit, data = data, type = NULL)
  data_used <- as.data.frame(data)[zcov_out$obs_id, , drop = FALSE]
  mu <- as.numeric(stats::predict(fit, newdata = data_used, type = "response"))
  fam <- fit$family$family

  M <- if (is.null(ndraws)) 1L else as.integer(ndraws[1L])
  if (!is.finite(M) || M < 1L) M <- 1L
  N <- length(mu)

  out <- switch(
    fam,
    gaussian = matrix(stats::rnorm(M * N, mean = rep(mu, each = M), sd = sqrt(summary(fit)$dispersion)), nrow = M, ncol = N),
    poisson = matrix(stats::rpois(M * N, lambda = rep(mu, each = M)), nrow = M, ncol = N),
    binomial = matrix(stats::rbinom(M * N, size = 1L, prob = rep(mu, each = M)), nrow = M, ncol = N),
    quasibinomial = matrix(stats::rbinom(M * N, size = 1L, prob = rep(mu, each = M)), nrow = M, ncol = N),
    quasipoisson = matrix(stats::rpois(M * N, lambda = rep(mu, each = M)), nrow = M, ncol = N),
    stop(sprintf("postpred_simulate.glm: simulation is not implemented for family='%s'.", fam), call. = FALSE)
  )

  out
}

.postpred_rtruncpois_matrix <- function(lambda) {
  lambda <- as.matrix(lambda)
  f0 <- stats::ppois(0, lambda = lambda)
  u <- f0 + stats::runif(length(lambda)) * pmax(1 - f0, .Machine$double.eps)
  matrix(stats::qpois(u, lambda = as.vector(lambda)), nrow = nrow(lambda), ncol = ncol(lambda))
}

.postpred_rtruncnb_matrix <- function(mu, size) {
  mu <- as.matrix(mu)
  size <- as.matrix(size)
  f0 <- stats::pnbinom(0, mu = mu, size = size)
  u <- f0 + stats::runif(length(mu)) * pmax(1 - f0, .Machine$double.eps)
  matrix(stats::qnbinom(u, mu = as.vector(mu), size = as.vector(size)), nrow = nrow(mu), ncol = ncol(mu))
}
