predcheck_pointpred_brms_hurdle_negbinomial <- function(fit,
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

  hu <- brms::posterior_epred(
    fit,
    newdata = data,
    dpar = "hu",
    draw_ids = draw_ids,
    transform = TRUE
  )
  hu <- pmin(pmax(hu, eps), 1 - eps)

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
      out <- matrix(NA_real_, nrow = nd, ncol = n)
      for (i in seq_len(nd)) {
        out[i, ] <- dhurdlenb(y_new, mu = mu[i, ], size = shape[i, ], pi = hu[i, ])
      }
      return(out)
    }
    dhurdlenb(y_new, mu = mu[draw, ], size = shape[draw, ], pi = hu[draw, ])
  }

  tail_fun <- function(y_new, draw = NULL) {
    y_new <- as.integer(y_new)
    if (is.null(draw)) {
      out <- matrix(NA_real_, nrow = nd, ncol = n)
      for (i in seq_len(nd)) {
        out[i, ] <- phurdlenb(y_new, mu = mu[i, ], size = shape[i, ], pi = hu[i, ], lower.tail = FALSE)
      }
      return(out)
    }
    phurdlenb(y_new, mu = mu[draw, ], size = shape[draw, ], pi = hu[draw, ], lower.tail = FALSE)
  }

  rng_fun <- function() {
    brms::posterior_predict(
      fit,
      newdata = data,
      draw_ids = draw_ids
    )
  }

  moment_fun <- function(draw = NULL) {
    calc_one <- function(mu_vec, hu_vec, shape_vec) {
      p0 <- (shape_vec / (mu_vec + shape_vec))^shape_vec
      p0 <- pmin(pmax(p0, eps), 1 - eps)

      mean_zt <- mu_vec / (1 - p0)
      var_nb <- mu_vec + mu_vec^2 / shape_vec
      ey2_nb <- var_nb + mu_vec^2
      ey2_zt <- ey2_nb / (1 - p0)
      var_zt <- ey2_zt - mean_zt^2

      mean_vec <- (1 - hu_vec) * mean_zt
      var_vec <- (1 - hu_vec) * var_zt + hu_vec * (1 - hu_vec) * mean_zt^2
      var_vec[!is.finite(var_vec) | var_vec <= 0] <- 1e-6

      list(mean = mean_vec, var = var_vec)
    }

    if (is.null(draw)) {
      out_mean <- matrix(NA_real_, nrow = nd, ncol = n)
      out_var <- matrix(NA_real_, nrow = nd, ncol = n)
      for (i in seq_len(nd)) {
        tmp <- calc_one(mu[i, ], hu[i, ], shape[i, ])
        out_mean[i, ] <- tmp$mean
        out_var[i, ] <- tmp$var
      }
      return(list(mean = out_mean, var = out_var))
    }

    calc_one(mu[draw, ], hu[draw, ], shape[draw, ])
  }

  list(
    support = "discrete",
    family = "hurdle_negbinomial",
    y = y,
    n = n,
    ndraws = nd,
    params = list(mu = mu, hu = hu, shape = shape),
    pmf = pmf_fun,
    tail = tail_fun,
    rng = rng_fun,
    moments = moment_fun
  )
}
