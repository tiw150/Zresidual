predcheck_pointpred_brms_negbinomial <- function(fit,
                                                 data,
                                                 draw_ids = NULL,
                                                 ...) {
  eps <- .Machine$double.eps

  if (is.null(draw_ids)) {
    draw_ids <- .predcheck_get_draw_ids(fit)
  }

  y <- .predcheck_extract_response(fit, data)

  mu <- brms::posterior_epred(
    fit,
    newdata = data,
    dpar = "mu",
    draw_ids = draw_ids
  )
  mu <- pmax(mu, eps)

  shape <- .predcheck_extract_shape_matrix(
    fit = fit,
    data = data,
    draw_ids = draw_ids,
    n_obs = nrow(data),
    eps = eps
  )

  nd <- nrow(mu)
  n <- ncol(mu)

  pmf_fun <- function(y_new, draw = NULL) {
    y_new <- as.integer(y_new)
    if (is.null(draw)) {
      y_mat <- matrix(rep(y_new, each = nd), nrow = nd)
      return(stats::dnbinom(y_mat, mu = mu, size = shape))
    }
    stats::dnbinom(y_new, mu = mu[draw, ], size = shape[draw, ])
  }

  tail_fun <- function(y_new, draw = NULL) {
    y_new <- as.integer(y_new)
    if (is.null(draw)) {
      y_mat <- matrix(rep(y_new, each = nd), nrow = nd)
      return(stats::pnbinom(y_mat, mu = mu, size = shape, lower.tail = FALSE))
    }
    stats::pnbinom(y_new, mu = mu[draw, ], size = shape[draw, ], lower.tail = FALSE)
  }

  rng_fun <- function() {
    brms::posterior_predict(
      fit,
      newdata = data,
      draw_ids = draw_ids
    )
  }

  moment_fun <- function(draw = NULL) {
    if (is.null(draw)) {
      var_mat <- mu + mu^2 / shape
      var_mat[!is.finite(var_mat) | var_mat <= 0] <- 1e-6
      return(list(mean = mu, var = var_mat))
    }

    mean_vec <- mu[draw, ]
    var_vec <- mean_vec + mean_vec^2 / shape[draw, ]
    var_vec[!is.finite(var_vec) | var_vec <= 0] <- 1e-6
    list(mean = mean_vec, var = var_vec)
  }

  list(
    support = "discrete",
    family = "negbinomial",
    y = y,
    n = n,
    ndraws = nd,
    params = list(mu = mu, shape = shape),
    pmf = pmf_fun,
    tail = tail_fun,
    rng = rng_fun,
    moments = moment_fun
  )
}
