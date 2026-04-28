predcheck_pointpred_brms_poisson <- function(fit,
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

  nd <- nrow(mu)
  n <- ncol(mu)

  pmf_fun <- function(y_new, draw = NULL) {
    y_new <- as.integer(y_new)
    if (is.null(draw)) {
      y_mat <- matrix(rep(y_new, each = nd), nrow = nd)
      return(stats::dpois(y_mat, lambda = mu))
    }
    stats::dpois(y_new, lambda = mu[draw, ])
  }

  tail_fun <- function(y_new, draw = NULL) {
    y_new <- as.integer(y_new)
    if (is.null(draw)) {
      y_mat <- matrix(rep(y_new, each = nd), nrow = nd)
      return(stats::ppois(y_mat, lambda = mu, lower.tail = FALSE))
    }
    stats::ppois(y_new, lambda = mu[draw, ], lower.tail = FALSE)
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
      return(list(mean = mu, var = mu))
    }
    list(mean = mu[draw, ], var = mu[draw, ])
  }

  list(
    support = "discrete",
    family = "poisson",
    y = y,
    n = n,
    ndraws = nd,
    params = list(mu = mu),
    pmf = pmf_fun,
    tail = tail_fun,
    rng = rng_fun,
    moments = moment_fun
  )
}
