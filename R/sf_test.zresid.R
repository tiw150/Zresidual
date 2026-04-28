#' Shapiro-Francia normality test for Z-residuals
#'
#' @description
#' Applies the Shapiro-Francia test column-wise to a matrix of Z-residuals.
#'
#' @param Zresidual A numeric vector or matrix of Z-residuals.
#'
#' @return A numeric vector of p-values, one per column.
#'
#' @examples
#' if (requireNamespace("nortest", quietly = TRUE)) {
#'   set.seed(1)
#'   z <- matrix(rnorm(60), ncol = 2)
#'   sf_test.zresid(z)
#' }
#'
#' @export sf_test.zresid
sf_test.zresid <- function(Zresidual) {
  if (!requireNamespace("nortest", quietly = TRUE)) {
    stop("Package 'nortest' is required for sf_test.zresid().", call. = FALSE)
  }
  
  Zresidual <- as.matrix(Zresidual)
  
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf   <- which(is.infinite(Zresidual) & Zresidual > 0)
  if (length(id.negtv.inf)) Zresidual[id.negtv.inf] <- -1e10
  if (length(id.pos.inf))   Zresidual[id.pos.inf]   <-  1e10
  
  id.nan <- which(is.nan(Zresidual))
  id.inf <- which(is.infinite(Zresidual))
  if (length(id.inf) > 0L) message("Non-finite Zresiduals exist! The model or fitting process may have a problem!")
  if (length(id.nan) > 0L) message("NaNs exist! The model or fitting process may have a problem!")
  
  sf.pv <- rep(NA_real_, ncol(Zresidual))
  for (i in seq_len(ncol(Zresidual))) {
    sf.pv[i] <- tryCatch(
      nortest::sf.test(Zresidual[, i])$p.value,
      error = function(e) NA_real_
    )
  }
  sf.pv
}
