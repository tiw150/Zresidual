log_pointpred_brmsfit_gaussian <- function(fit, data, type = NULL, ...) {
  if (!inherits(fit, "brmsfit")) {
    stop("log_pointpred_brmsfit_gaussian: `fit` must be a brmsfit object.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_gaussian: `data` must be provided.", call. = FALSE)
  }

  fam <- tryCatch(fit$family$family, error = function(e) NULL)
  if (!is.null(fam) && !identical(fam, "gaussian")) {
    stop("log_pointpred_brmsfit_gaussian: fitted brms family must be 'gaussian'.", call. = FALSE)
  }

  data <- as.data.frame(data)
  n <- nrow(data)

  if (n < 1L) {
    stop("log_pointpred_brmsfit_gaussian: `data` has zero rows.", call. = FALSE)
  }

  formula_obj <- tryCatch(fit$formula$formula, error = function(e) NULL)
  if (is.null(formula_obj)) {
    formula_obj <- fit$formula
  }

  model.data <- stats::model.frame(
    formula_obj,
    data = data,
    na.action = stats::na.pass
  )

  response <- tryCatch(fit$formula$resp, error = function(e) NULL)
  if (is.null(response) || !(response %in% names(model.data))) {
    response <- all.vars(formula_obj)[1L]
  }

  if (is.null(response) || !(response %in% names(model.data))) {
    stop("log_pointpred_brmsfit_gaussian: response variable not found in `data`.", call. = FALSE)
  }

  y <- as.numeric(model.data[[response]])

  if (length(y) != n) {
    stop("log_pointpred_brmsfit_gaussian: response length does not match nrow(data).", call. = FALSE)
  }

  if (anyNA(y) || any(!is.finite(y))) {
    stop("log_pointpred_brmsfit_gaussian: response contains NA or non-finite values.", call. = FALSE)
  }

  as_param_matrix <- function(x, name, n_obs, ndraws = NULL) {
    if (is.null(dim(x))) {
      x <- matrix(as.numeric(x), nrow = 1L)
    } else {
      x <- as.matrix(x)
    }

    if (nrow(x) < 1L || ncol(x) < 1L) {
      stop(sprintf("log_pointpred_brmsfit_gaussian: `%s` is empty.", name), call. = FALSE)
    }

    if (!is.null(ndraws) && nrow(x) != ndraws) {
      stop(sprintf("log_pointpred_brmsfit_gaussian: `%s` draw dimension mismatch.", name), call. = FALSE)
    }

    if (ncol(x) == 1L && n_obs > 1L) {
      x <- matrix(rep(x[, 1L], times = n_obs), nrow = nrow(x), ncol = n_obs)
    }

    if (ncol(x) != n_obs) {
      stop(sprintf("log_pointpred_brmsfit_gaussian: `%s` column dimension must match nrow(data).", name), call. = FALSE)
    }

    if (anyNA(x) || any(!is.finite(x))) {
      stop(sprintf("log_pointpred_brmsfit_gaussian: `%s` contains NA or non-finite values.", name), call. = FALSE)
    }

    x
  }

  mu <- brms::posterior_epred(
    fit,
    newdata = data,
    ...
  )
  mu <- as_param_matrix(mu, name = "mu", n_obs = n)

  ndraws <- nrow(mu)

  sigma <- tryCatch(
    brms::posterior_linpred(
      fit,
      newdata = data,
      dpar = "sigma",
      transform = TRUE,
      ...
    ),
    error = function(e) NULL
  )

  if (is.null(sigma)) {
    if (!requireNamespace("posterior", quietly = TRUE)) {
      stop("log_pointpred_brmsfit_gaussian: package 'posterior' is required to extract sigma.", call. = FALSE)
    }

    dd <- posterior::as_draws_df(fit)

    if (!("sigma" %in% names(dd))) {
      stop("log_pointpred_brmsfit_gaussian: cannot extract posterior draws for `sigma`.", call. = FALSE)
    }

    sigma_vec <- as.numeric(dd[["sigma"]])
    sigma <- matrix(sigma_vec, nrow = length(sigma_vec), ncol = n)
  }

  sigma <- as_param_matrix(sigma, name = "sigma", n_obs = n, ndraws = ndraws)

  if (any(sigma <= 0)) {
    stop("log_pointpred_brmsfit_gaussian: `sigma` must be strictly positive.", call. = FALSE)
  }

  y_mat <- matrix(
    rep(y, each = ndraws),
    nrow = ndraws,
    ncol = n,
    byrow = FALSE
  )

  log_like <- matrix(
    stats::dnorm(
      x = y_mat,
      mean = mu,
      sd = sigma,
      log = TRUE
    ),
    nrow = ndraws,
    ncol = n
  )

  log_surv <- matrix(
    stats::pnorm(
      q = y_mat,
      mean = mu,
      sd = sigma,
      lower.tail = FALSE,
      log.p = TRUE
    ),
    nrow = ndraws,
    ncol = n
  )

  if (!identical(dim(log_like), dim(mu))) {
    stop("log_pointpred_brmsfit_gaussian: `log_like` dimension mismatch.", call. = FALSE)
  }

  if (!identical(dim(log_surv), dim(mu))) {
    stop("log_pointpred_brmsfit_gaussian: `log_surv` dimension mismatch.", call. = FALSE)
  }

  if (anyNA(log_like)) {
    stop("log_pointpred_brmsfit_gaussian: `log_like` contains NA.", call. = FALSE)
  }

  if (anyNA(log_surv)) {
    stop("log_pointpred_brmsfit_gaussian: `log_surv` contains NA.", call. = FALSE)
  }

  is_discrete <- matrix(0L, nrow = 1L, ncol = n)

  list(
    log_surv = log_surv,
    log_like = log_like,
    is_discrete = is_discrete
  )
}
