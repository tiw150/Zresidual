# =========================================================
# Helpers for predictive checks
# =========================================================

.predcheck_clip_prob <- function(p, eps = .Machine$double.eps) {
  pmin(pmax(p, eps), 1 - eps)
}

.predcheck_clip_logprob <- function(logp, eps = 1e-12) {
  lo <- log(eps)
  hi <- log1p(-eps)
  pmin(pmax(logp, lo), hi)
}

.predcheck_log_add_exp2 <- function(a, b) {
  both_ninf <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[both_ninf] <- -Inf

  a_ninf <- is.infinite(a) & a < 0 & !both_ninf
  b_ninf <- is.infinite(b) & b < 0 & !both_ninf
  out[a_ninf] <- b[a_ninf]
  out[b_ninf] <- a[b_ninf]

  out
}

.predcheck_as_matrix <- function(x, name, nrow_default = 1L) {
  if (is.null(x)) {
    stop(sprintf("predcheck: missing `%s`.", name), call. = FALSE)
  }

  if (is.matrix(x)) return(x)
  if (is.data.frame(x)) return(as.matrix(x))
  if (is.vector(x) || is.factor(x)) return(matrix(as.vector(x), nrow = nrow_default))

  stop(sprintf("predcheck: `%s` must be a vector or matrix.", name), call. = FALSE)
}

.predcheck_get_draw_ids <- function(fit, ndraws = NULL) {
  if (!inherits(fit, "brmsfit")) {
    if (is.null(ndraws)) return(NULL)

    ndraws <- suppressWarnings(as.integer(ndraws[1L]))
    if (!is.finite(ndraws) || ndraws < 1L) {
      stop("predcheck: `ndraws` must be a positive integer.", call. = FALSE)
    }

    return(seq_len(ndraws))
  }

  nd_all <- tryCatch(brms::ndraws(fit), error = function(e) nrow(as.matrix(fit)))
  if (!is.finite(nd_all) || nd_all < 1L) {
    stop("predcheck: no posterior draws found in `fit`.", call. = FALSE)
  }

  if (is.null(ndraws)) return(seq_len(nd_all))

  ndraws <- suppressWarnings(as.integer(ndraws[1L]))
  if (!is.finite(ndraws) || ndraws < 1L) {
    stop("predcheck: `ndraws` must be a positive integer.", call. = FALSE)
  }

  sort(sample.int(nd_all, size = min(ndraws, nd_all), replace = FALSE))
}

.predcheck_resolve_log_pointpred <- function(log_pointpred = NULL,
                                             fit,
                                             pred_method = "analytic",
                                             type = NULL) {
  if (!is.null(log_pointpred)) {
    if (is.function(log_pointpred)) return(log_pointpred)
    if (is.character(log_pointpred) && length(log_pointpred) == 1L) return(match.fun(log_pointpred))
    stop("predcheck: `log_pointpred` must be NULL, a function, or a function name.", call. = FALSE)
  }

  fn <- required_log_pointpred(
    fit = fit,
    pred_method = pred_method,
    type = type,
    show_names = FALSE
  )

  if (is.character(fn)) stop(fn, call. = FALSE)
  fn
}

.predcheck_call_log_pointpred <- function(fn,
                                          fit,
                                          data,
                                          type = NULL,
                                          draw_ids = NULL,
                                          ...) {
  fmls <- tryCatch(names(formals(fn)), error = function(e) character(0))
  args <- list(fit = fit, data = data)

  if ("type" %in% fmls || "..." %in% fmls) args$type <- type
  if ("draw_ids" %in% fmls || "..." %in% fmls) args$draw_ids <- draw_ids

  do.call(fn, c(args, list(...)))
}

.predcheck_resolve_point_Zcov <- function(point_Zcov = NULL) {
  if (!is.null(point_Zcov)) {
    if (is.function(point_Zcov)) return(point_Zcov)
    if (is.character(point_Zcov) && length(point_Zcov) == 1L) return(match.fun(point_Zcov))
    stop("predcheck: `point_Zcov` must be NULL, a function, or a function name.", call. = FALSE)
  }

  get("point_Zcov", mode = "function", inherits = TRUE)
}

.predcheck_resolve_postpred_simulate <- function(postpred_simulate = NULL) {
  if (!is.null(postpred_simulate)) {
    if (is.function(postpred_simulate)) return(postpred_simulate)
    if (is.character(postpred_simulate) && length(postpred_simulate) == 1L) return(match.fun(postpred_simulate))
    stop("predcheck: `postpred_simulate` must be NULL, a function, or a function name.", call. = FALSE)
  }

  get("postpred_simulate", mode = "function", inherits = TRUE)
}

.predcheck_validate_log_pointpred <- function(pp) {
  if (!is.list(pp) || is.null(pp$log_surv) || is.null(pp$log_like) || is.null(pp$is_discrete)) {
    stop("predcheck: `log_pointpred` must return `log_surv`, `log_like`, and `is_discrete`.", call. = FALSE)
  }

  log_surv <- .predcheck_as_matrix(pp$log_surv, "log_surv")
  log_like <- .predcheck_as_matrix(pp$log_like, "log_like")
  is_disc  <- .predcheck_as_matrix(pp$is_discrete, "is_discrete")

  if (!identical(dim(log_surv), dim(log_like))) {
    stop("predcheck: `log_surv` and `log_like` must have identical dimensions.", call. = FALSE)
  }

  if (nrow(is_disc) != 1L || ncol(is_disc) != ncol(log_surv)) {
    stop("predcheck: `is_discrete` must be a length-N vector or a 1 x N matrix.", call. = FALSE)
  }

  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = as.integer(is_disc[1L, ])
  )
}

.predcheck_compute_z_matrix <- function(log_pred,
                                        randomized = TRUE,
                                        eps = 1e-12) {
  pp <- .predcheck_validate_log_pointpred(log_pred)

  log_surv <- pp$log_surv
  log_like <- pp$log_like
  is_disc  <- pp$is_discrete

  M <- nrow(log_surv)
  N <- ncol(log_surv)

  log_rsp <- log_surv
  idx_disc <- which(is_disc == 1L)

  if (length(idx_disc) > 0L) {
    u <- if (isTRUE(randomized)) {
      matrix(stats::runif(M * length(idx_disc)), nrow = M, ncol = length(idx_disc))
    } else {
      matrix(0.5, nrow = M, ncol = length(idx_disc))
    }

    log_rsp[, idx_disc] <- .predcheck_log_add_exp2(
      log_surv[, idx_disc, drop = FALSE],
      log_like[, idx_disc, drop = FALSE] + log(pmax(u, .Machine$double.xmin))
    )
  }

  log_rsp <- .predcheck_clip_logprob(log_rsp, eps = eps)
  z <- -stats::qnorm(log_rsp, log.p = TRUE)

  if (!is.matrix(z) || !identical(dim(z), c(M, N))) {
    stop("predcheck: failed to compute a valid Z-residual matrix.", call. = FALSE)
  }

  t(z)
}

.predcheck_expand_x_specs <- function(base_name, x, k_anova = 10, k_bl = 10) {
  if (is.null(x)) return(list())

  # Character vectors are interpreted as multiple covariate names, e.g.
  # x = c("lp", "x1") creates two separate AOV tests.
  if (is.character(x) && length(x) > 1L) {
    return(lapply(x, function(xi) {
      spec <- list(name = base_name, x = xi)
      if (identical(base_name, "aov")) spec$k.anova <- k_anova
      if (identical(base_name, "bartlett")) spec$k.bl <- k_bl
      spec
    }))
  }

  spec <- list(name = base_name, x = x)
  if (identical(base_name, "aov")) spec$k.anova <- k_anova
  if (identical(base_name, "bartlett")) spec$k.bl <- k_bl
  list(spec)
}

.predcheck_test_label <- function(spec) {
  nm <- tolower(as.character(spec$name)[1L])
  base <- switch(nm, aov = "AOV", bartlett = "BL", sw = "SW", toupper(nm))

  if (nm %in% c("aov", "bartlett") && !is.null(spec$x)) {
    x <- spec$x
    if (is.character(x) && length(x) == 1L) {
      return(paste0(base, "[", x, "]"))
    }
  }

  base
}

.predcheck_normalize_tests <- function(test = NULL,
                                       x = NULL,
                                       k_anova = 10,
                                       k_bl = 10) {
  if (is.null(test)) {
    test <- list(list(name = "sw"))
    if (!is.null(x)) {
      test <- c(test, .predcheck_expand_x_specs("aov", x, k_anova = k_anova, k_bl = k_bl))
    }
    return(test)
  }

  if (is.list(test) && !is.null(test$name)) test <- list(test)
  if (is.list(test) && !is.null(test$func)) test <- list(test)
  if (!is.list(test)) stop("predcheck: `test` must be a list of test specifications.", call. = FALSE)

  out <- list()

  for (spec in test) {
    if (!is.list(spec)) stop("predcheck: each test specification must be a list.", call. = FALSE)
    if (is.null(spec$name) && !is.null(spec$func)) spec$name <- spec$func
    if (is.null(spec$name)) stop("predcheck: each test specification needs `name` or `func`.", call. = FALSE)
    if (is.function(spec$name)) {
      spec$func <- spec$name
      spec$name <- "custom"
    }

    if (identical(spec$name, "anova")) spec$name <- "aov"
    if (identical(spec$name, "bl")) spec$name <- "bartlett"

    if (identical(spec$name, "aov") && is.null(spec$k.anova)) spec$k.anova <- k_anova
    if (identical(spec$name, "bartlett") && is.null(spec$k.bl)) spec$k.bl <- k_bl

    # Also support explicit test specs such as
    # test = list(list(name = "aov", x = c("lp", "x1"))).
    if (identical(spec$name, "aov") && is.character(spec$x) && length(spec$x) > 1L) {
      out <- c(out, lapply(spec$x, function(xi) {
        si <- spec
        si$x <- xi
        si
      }))
    } else if (identical(spec$name, "bartlett") && is.character(spec$x) && length(spec$x) > 1L) {
      out <- c(out, lapply(spec$x, function(xi) {
        si <- spec
        si$x <- xi
        si
      }))
    } else {
      out <- c(out, list(spec))
    }
  }

  out
}

.predcheck_vector_from_zcov <- function(zcov_out,
                                        name,
                                        draw = NULL,
                                        n_draws = NULL,
                                        n_obs = NULL) {
  if (is.null(name)) return(NULL)

  nm <- as.character(name)[1L]
  nm_low <- tolower(nm)

  if (nm_low %in% c("index", "obs", "observation")) {
    if (is.null(n_obs)) return(NULL)
    return(seq_len(n_obs))
  }

  if (nm_low %in% c("linear_pred", "linear.pred")) nm_low <- "lp"
  if (nm_low %in% c("var")) nm_low <- "variance"

  val <- NULL

  if (!is.null(zcov_out[[nm_low]])) {
    val <- zcov_out[[nm_low]]
  } else if (identical(nm_low, "lp") && !is.null(zcov_out$linear_pred)) {
    val <- zcov_out$linear_pred
  }

  if (!is.null(val)) {
    if (is.matrix(val)) {
      if (!is.null(draw) && !is.null(n_draws) && nrow(val) == n_draws) return(as.vector(val[draw, ]))
      if (!is.null(draw) && !is.null(n_draws) && ncol(val) == n_draws) return(as.vector(val[, draw]))
      if (nrow(val) == 1L) return(as.vector(val[1L, ]))
      if (ncol(val) == 1L) return(as.vector(val[, 1L]))
      return(as.vector(val))
    }
    return(as.vector(val))
  }

  cov_mat <- zcov_out$covariate
  if (!is.null(cov_mat)) {
    cov_mat <- as.matrix(cov_mat)
    cn <- colnames(cov_mat)
    if (!is.null(cn) && nm %in% cn) return(as.vector(cov_mat[, nm]))
    if (!is.null(cn) && nm_low %in% tolower(cn)) return(as.vector(cov_mat[, which(tolower(cn) == nm_low)[1L]]))
  }

  cov_df <- zcov_out$covariates
  if (!is.null(cov_df)) {
    cov_df <- as.data.frame(cov_df)
    cn <- names(cov_df)
    if (!is.null(cn) && nm %in% cn) return(cov_df[[nm]])
    if (!is.null(cn) && nm_low %in% tolower(cn)) return(cov_df[[which(tolower(cn) == nm_low)[1L]]])
  }

  NULL
}



.predcheck_eval_one_test <- function(z,
                                     zcov_out,
                                     spec,
                                     draw,
                                     n_draws) {
  if (is.function(spec$func) && identical(as.character(spec$name)[1L], "custom")) {
    args <- spec[names(spec) %in% names(formals(spec$func))]
    args$name <- NULL
    args$func <- NULL
    return(as.numeric(do.call(spec$func, c(list(z = z), args)))[1L])
  }

  name <- tolower(as.character(spec$name)[1L])

  if (name %in% c("sw", "shapiro", "shapiro_wilk")) {
    return(as.numeric(sw.test.zresid(z))[1L])
  }

  if (name %in% c("aov", "anova", "nl_aov", "nonlinear_anova")) {
    x_name <- spec$x %||% spec$X %||% "lp"
    xv <- .predcheck_vector_from_zcov(
      zcov_out = zcov_out,
      name = x_name,
      draw = draw,
      n_draws = n_draws,
      n_obs = length(z)
    )
    if (is.null(xv)) {
      if (length(x_name) == length(z)) xv <- x_name
      else stop(sprintf("predcheck: cannot find x='%s' in point_Zcov output.", as.character(x_name)[1L]), call. = FALSE)
    }
    return(test.nl.aov(z, xv, k.anova = spec$k.anova %||% 10))
  }

  if (name %in% c("bartlett", "bl", "variance")) {
    x_name <- spec$x %||% spec$X %||% "lp"
    xv <- .predcheck_vector_from_zcov(
      zcov_out = zcov_out,
      name = x_name,
      draw = draw,
      n_draws = n_draws,
      n_obs = length(z)
    )
    if (is.null(xv)) {
      if (length(x_name) == length(z)) xv <- x_name
      else stop(sprintf("predcheck: cannot find x='%s' in point_Zcov output.", as.character(x_name)[1L]), call. = FALSE)
    }
    return(test.var.bartl(z, xv, k.bl = spec$k.bl %||% 10))
  }

  if (is.function(spec$func)) {
    args <- spec[names(spec) %in% names(formals(spec$func))]
    args$name <- NULL
    args$func <- NULL
    return(as.numeric(do.call(spec$func, c(list(z = z), args)))[1L])
  }

  stop(sprintf("predcheck: unknown test name '%s'.", name), call. = FALSE)
}

.predcheck_matrix_draw <- function(mat, draw, n_draws, name) {
  mat <- .predcheck_as_matrix(mat, name)

  if (nrow(mat) == n_draws) return(as.numeric(mat[draw, ]))
  if (ncol(mat) == n_draws) return(as.numeric(mat[, draw]))
  if (nrow(mat) == 1L) return(as.numeric(mat[1L, ]))

  stop(sprintf("predcheck: `%s` has incompatible dimensions.", name), call. = FALSE)
}

.predcheck_recycle_matrix_to_draws <- function(mat, n_draws, name) {
  mat <- .predcheck_as_matrix(mat, name)
  if (nrow(mat) == n_draws) return(mat)
  if (nrow(mat) == 1L) {
    return(matrix(rep(mat[1L, ], each = n_draws), nrow = n_draws, ncol = ncol(mat), byrow = FALSE))
  }
  stop(sprintf("predcheck: `%s` has incompatible draw dimension.", name), call. = FALSE)
}



# Fast replicated Z-residual path for selected brms families.
# This avoids calling log_pointpred() once per posterior draw in predcheck_zresid().
# Supported families currently include:
#   - hurdle_poisson / hurdle_negbinomial, with type = hurdle/count/zero
#   - poisson
#   - negbinomial
#   - bernoulli
#   - gaussian
#   - inverse.gaussian
.predcheck_fast_rep_z_matrix <- function(fit,
                                         data,
                                         zcov_out,
                                         y_rep,
                                         type = NULL,
                                         draw_ids = NULL,
                                         randomized = TRUE,
                                         eps = 1e-12,
                                         ...) {
  if (!inherits(fit, "brmsfit")) return(NULL)
  if (!requireNamespace("brms", quietly = TRUE)) return(NULL)

  fam <- tryCatch(fit$family$family, error = function(e) NULL)
  if (is.null(fam)) return(NULL)

  y_rep <- .predcheck_as_matrix(y_rep, "y_rep")
  M <- nrow(y_rep)
  N <- ncol(y_rep)

  obs_id <- zcov_out$obs_id
  if (is.null(obs_id)) obs_id <- seq_len(nrow(as.data.frame(data)))
  obs_id <- as.integer(obs_id)

  if (length(obs_id) != N) {
    return(NULL)
  }

  data_used <- as.data.frame(data)[obs_id, , drop = FALSE]

  get_par <- function(dpar) {
    out <- posterior.pred(
      fit = fit,
      dpar = dpar,
      data = data_used,
      count.only = FALSE,
      draw_ids = draw_ids,
      ...
    )
    out <- as.matrix(out)
    if (!identical(dim(out), c(M, N))) {
      stop(
        sprintf(
          "predcheck_zresid fast path: parameter '%s' has dimension %s but expected %s.",
          dpar,
          paste(dim(out), collapse = " x "),
          paste(c(M, N), collapse = " x ")
        ),
        call. = FALSE
      )
    }
    out
  }

  y_vec <- as.numeric(y_rep)
  if (any(!is.finite(y_vec))) {
    stop("predcheck_zresid fast path: `y_rep` contains non-finite values.", call. = FALSE)
  }

  log_like <- NULL
  log_surv <- NULL
  is_discrete <- 1L

  type_req <- zcov_out$type %||% type

  if (grepl("^hurdle", fam)) {
    if (!grepl("poisson|negbinomial", fam)) return(NULL)
    if (is.null(type_req)) type_req <- "hurdle"

    if (identical(type_req, "zero")) {
      hu <- get_par("zero")
      o_vec <- as.integer(y_vec >= 0.5)

      log_like <- matrix(
        stats::dbinom(
          x = o_vec,
          size = 1L,
          prob = as.numeric(hu),
          log = TRUE
        ),
        nrow = M,
        ncol = N
      )

      log_surv <- matrix(
        stats::pbinom(
          q = o_vec,
          size = 1L,
          prob = as.numeric(hu),
          lower.tail = FALSE,
          log.p = TRUE
        ),
        nrow = M,
        ncol = N
      )
    } else if (identical(type_req, "count")) {
      mu <- get_par("mu")

      if (grepl("negbinomial", fam)) {
        shape <- get_par("shape")
        log_like <- matrix(
          dtruncnb(
            y = y_vec,
            mu = as.numeric(mu),
            size = as.numeric(shape),
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
        log_surv <- matrix(
          ptruncnb(
            y = y_vec,
            mu = as.numeric(mu),
            size = as.numeric(shape),
            lower.tail = FALSE,
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
      } else {
        log_like <- matrix(
          dtruncpois(
            y = y_vec,
            lambda = as.numeric(mu),
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
        log_surv <- matrix(
          ptruncpois(
            y = y_vec,
            lambda = as.numeric(mu),
            lower.tail = FALSE,
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
      }
    } else {
      hu <- get_par("zero")
      mu <- get_par("mu")

      if (grepl("negbinomial", fam)) {
        shape <- get_par("shape")
        log_like <- matrix(
          dhurdlenb(
            y = y_vec,
            mu = as.numeric(mu),
            size = as.numeric(shape),
            pi = as.numeric(hu),
            log = TRUE
          ),
          nrow = M,
          ncol = N
        )
        log_surv <- matrix(
          phurdlenb(
            y = y_vec,
            mu = as.numeric(mu),
            size = as.numeric(shape),
            pi = as.numeric(hu),
            lower.tail = FALSE,
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
      } else {
        log_like <- matrix(
          dhurdlepois(
            y = y_vec,
            lambda = as.numeric(mu),
            pi = as.numeric(hu),
            log = TRUE
          ),
          nrow = M,
          ncol = N
        )
        log_surv <- matrix(
          phurdlepois(
            y = y_vec,
            lambda = as.numeric(mu),
            pi = as.numeric(hu),
            lower.tail = FALSE,
            log.p = TRUE
          ),
          nrow = M,
          ncol = N
        )
      }
    }
  } else if (identical(fam, "negbinomial")) {
    mu <- get_par("mu")
    shape <- get_par("shape")

    log_like <- matrix(
      stats::dnbinom(
        x = y_vec,
        size = as.numeric(shape),
        mu = as.numeric(mu),
        log = TRUE
      ),
      nrow = M,
      ncol = N
    )

    log_surv <- matrix(
      stats::pnbinom(
        q = y_vec,
        size = as.numeric(shape),
        mu = as.numeric(mu),
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = M,
      ncol = N
    )
  } else if (identical(fam, "poisson")) {
    mu <- get_par("mu")

    log_like <- matrix(
      stats::dpois(
        x = y_vec,
        lambda = as.numeric(mu),
        log = TRUE
      ),
      nrow = M,
      ncol = N
    )

    log_surv <- matrix(
      stats::ppois(
        q = y_vec,
        lambda = as.numeric(mu),
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = M,
      ncol = N
    )
  } else if (identical(fam, "bernoulli")) {
    p <- get_par("mu")
    y01 <- as.integer(y_vec >= 0.5)

    log_like <- matrix(
      stats::dbinom(
        x = y01,
        size = 1L,
        prob = as.numeric(p),
        log = TRUE
      ),
      nrow = M,
      ncol = N
    )

    log_surv <- matrix(
      stats::pbinom(
        q = y01,
        size = 1L,
        prob = as.numeric(p),
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = M,
      ncol = N
    )
  } else if (identical(fam, "gaussian")) {
    mu <- get_par("mu")
    sigma <- get_par("sigma")
    is_discrete <- 0L

    log_like <- matrix(
      stats::dnorm(
        x = y_vec,
        mean = as.numeric(mu),
        sd = as.numeric(sigma),
        log = TRUE
      ),
      nrow = M,
      ncol = N
    )

    log_surv <- matrix(
      stats::pnorm(
        q = y_vec,
        mean = as.numeric(mu),
        sd = as.numeric(sigma),
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = M,
      ncol = N
    )
  } else if (identical(fam, "inverse.gaussian")) {
    mu <- get_par("mu")
    shape <- get_par("shape")
    is_discrete <- 0L

    if (any(y_vec <= 0)) {
      stop("predcheck_zresid fast path: inverse Gaussian replicated responses must be positive.", call. = FALSE)
    }

    log_like <- matrix(
      dinvgauss(
        x = y_vec,
        mu = as.numeric(mu),
        shape = as.numeric(shape),
        log = TRUE
      ),
      nrow = M,
      ncol = N
    )

    log_surv <- matrix(
      pinvgauss(
        q = y_vec,
        mu = as.numeric(mu),
        shape = as.numeric(shape),
        lower.tail = FALSE,
        log.p = TRUE
      ),
      nrow = M,
      ncol = N
    )
  } else {
    return(NULL)
  }

  if (!identical(dim(log_like), c(M, N)) || !identical(dim(log_surv), c(M, N))) {
    stop("predcheck_zresid fast path: log matrix dimension mismatch.", call. = FALSE)
  }

  if (anyNA(log_like) || anyNA(log_surv)) {
    stop("predcheck_zresid fast path: `log_like` or `log_surv` contains NA.", call. = FALSE)
  }

  .predcheck_compute_z_matrix(
    log_pred = list(
      log_surv = log_surv,
      log_like = log_like,
      is_discrete = matrix(is_discrete, nrow = 1L, ncol = N)
    ),
    randomized = randomized,
    eps = eps
  )
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
