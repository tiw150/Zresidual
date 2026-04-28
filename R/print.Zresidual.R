#' Print methods for Z-residual objects
#'
#' @description
#' Compact print methods for \code{"zresid"} and \code{"cvzresid"} objects.
#'
#' @param x An object returned by \code{\link{Zresidual}} or
#'   \code{\link{CV.Zresidual}}.
#' @param ... Further arguments passed to \code{print()}.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 20
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.2 * x))
#'   t_cens  <- rexp(n, rate = 0.5)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'
#'   z <- Zresidual(fit, data=dat, nrep = 1, seed = 1)
#'   print(z)
#' }
#'
#' @name print.zresid
#' @export
print.zresid <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 2L) {
    cat(sprintf("<zresid> n=%d, reps=%d\n", d[1], d[2]))
    print(utils::head(unclass(x)), ...)
  } else {
    cat("<zresid>\n")
    print(unclass(x), ...)
  }
  invisible(x)
}

#' @rdname print.zresid
#' @export
print.cvzresid <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 2L) {
    cat(sprintf("<cvzresid> n=%d, reps=%d\n", d[1], d[2]))
    print(utils::head(unclass(x)), ...)
  } else {
    cat("<cvzresid>\n")
    print(unclass(x), ...)
  }
  invisible(x)
}
