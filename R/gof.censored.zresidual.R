#' Goodness-of-fit test for censored Z-residuals
#'
#' @description
#' Applies a goodness-of-fit test to censored Z-residuals using
#' \code{EnvStats::gofTestCensored}. This is typically used to assess whether
#' censored Z-residuals are approximately standard normal under a fitted
#' survival model.
#'
#' @param censored.Zresidual Numeric vector or one-column matrix of censored
#'   Z-residuals. It must carry an attribute \code{"censored.status"}.
#' @param test Character string passed to \code{EnvStats::gofTestCensored()}.
#'   The default is \code{"sf"}.
#' @param ... Additional arguments passed to
#'   \code{EnvStats::gofTestCensored()}.
#'
#' @return A single numeric p-value.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE) &&
#'     requireNamespace("EnvStats", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.3 * x))
#'   t_cens  <- rexp(n, rate = 0.5)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   rz <- surv_residuals(fit, data = dat, residual.type = "censored Z-residual")
#'   gof.censored.zresidual(rz)
#' }
#'
#' @seealso \code{\link[EnvStats]{gofTestCensored}}
#'
#' @export
gof.censored.zresidual <- function(censored.Zresidual, test = "sf", ...) {
  censored.status <- attr(censored.Zresidual, "censored.status")
  if (is.null(censored.status)) {
    stop("`censored.status` attribute is missing.", call. = FALSE)
  }
  
  if (!requireNamespace("EnvStats", quietly = TRUE)) {
    stop("Package 'EnvStats' is required for gof.censored.zresidual().", call. = FALSE)
  }
  
  x <- censored.Zresidual
  if (is.matrix(x)) {
    if (ncol(x) != 1L) {
      stop("`censored.Zresidual` must be a numeric vector or a one-column matrix.",
           call. = FALSE)
    }
    x <- x[, 1L]
  }
  x <- as.numeric(x)
  
  if (length(censored.status) != length(x)) {
    stop("`censored.status` must have the same length as `censored.Zresidual`.",
         call. = FALSE)
  }
  
  ## In survival objects, status is usually 1 = event, 0 = censored.
  ## EnvStats::gofTestCensored expects TRUE/FALSE for whether each value is censored.
  censored <- if (is.logical(censored.status)) {
    censored.status
  } else {
    as.integer(censored.status) == 0L
  }
  
  out <- EnvStats::gofTestCensored(
    x = x,
    censored = censored,
    test = test,
    ...
  )
  
  unname(out$p.value)
}

