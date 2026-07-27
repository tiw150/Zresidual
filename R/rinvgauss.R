# =========================================================
# Inverse Gaussian random generation
# =========================================================

# Generate inverse Gaussian random variables using the
# Michael--Schucany--Haas method.
#
# Parameterization:
#   mu    = mean parameter, mu > 0
#   shape = shape parameter, shape > 0
#
# This function is used by postpred_simulate.brmsfit() for
# brms inverse.gaussian models instead of brms::posterior_predict(),
# because the latter can occasionally return non-positive replicated
# values in extreme numerical cases.
rinvgauss <- function(n = NULL, mu, shape) {
  if (missing(mu) || missing(shape)) {
    stop("rinvgauss: `mu` and `shape` must be supplied.", call. = FALSE)
  }

  mu_dim <- dim(mu)
  shape_dim <- dim(shape)

  if (!is.null(mu_dim) && !is.null(shape_dim) && !identical(mu_dim, shape_dim)) {
    stop("rinvgauss: `mu` and `shape` must have identical dimensions.", call. = FALSE)
  }

  mu_vec <- as.numeric(mu)
  shape_vec <- as.numeric(shape)

  if (is.null(n)) {
    n_out <- max(length(mu_vec), length(shape_vec))
  } else {
    n_out <- as.integer(n[1L])
    if (!is.finite(n_out) || n_out < 0L) {
      stop("rinvgauss: `n` must be a non-negative integer.", call. = FALSE)
    }
  }

  if (n_out == 0L) {
    return(numeric(0L))
  }

  mu_vec <- rep(mu_vec, length.out = n_out)
  shape_vec <- rep(shape_vec, length.out = n_out)

  if (any(!is.finite(mu_vec)) || any(!is.finite(shape_vec))) {
    stop("rinvgauss: `mu` and `shape` must be finite.", call. = FALSE)
  }
  if (any(mu_vec <= 0) || any(shape_vec <= 0)) {
    stop("rinvgauss: `mu` and `shape` must be positive.", call. = FALSE)
  }

  v <- stats::rnorm(n_out)
  y <- v * v

  # Candidate root from the Michael--Schucany--Haas algorithm.
  x <- mu_vec + (mu_vec * mu_vec * y) / (2 * shape_vec) -
    (mu_vec / (2 * shape_vec)) *
      sqrt(pmax(0, 4 * mu_vec * shape_vec * y + mu_vec * mu_vec * y * y))

  # Numerical guard: x should be strictly positive mathematically, but can
  # underflow to 0 in extreme cases. Keep it inside the distribution support.
  x <- pmax(x, .Machine$double.xmin)

  u <- stats::runif(n_out)
  out <- ifelse(u <= mu_vec / (mu_vec + x), x, (mu_vec * mu_vec) / x)
  out <- pmax(out, .Machine$double.xmin)

  if (is.null(n) && !is.null(mu_dim)) {
    dim(out) <- mu_dim
  }

  out
}
