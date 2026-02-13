#' Shapiro-Wilk Normality Test for Z-Residuals
#'
#' Performs the Shapiro-Wilk test for normality on each column of a matrix of Z-residuals.
#'
#' @param Zresidual A numeric matrix of Z-residuals, where each column represents
#'   a separate set of residuals (e.g., from different posterior predictive draws or variables).
#' @param ... Additional arguments (ignored unless named):
#'   \describe{
#'     \item{max_n}{Integer. Maximum sample size passed to \code{shapiro.test} (default 5000).}
#'     \item{seed}{Integer. Optional random seed used when subsampling is needed.}
#'   }
#'
#' @return A numeric vector of Shapiro-Wilk p-values, one for each column of \code{Zresidual}.
#'
#' @details
#' Non-finite values are handled by replacing \code{Inf/-Inf} with a large finite value
#' based on the maximum finite magnitude in the same column. \code{NA/NaN} are removed
#' before testing. If the column has fewer than 3 finite values, or if the test fails,
#' \code{NA} is returned for that column.
#'
#' @export sw.test.zresid
sw.test.zresid <- function(Zresidual, ...) {
  args <- list(...)
  
  # Only honor named args to avoid breaking calls like sw.test.zresid(Z, X, k)
  max_n <- 5000L
  if (!is.null(args$max_n)) {
    max_n <- suppressWarnings(as.integer(args$max_n))
    if (!is.finite(max_n) || max_n < 3L) max_n <- 5000L
  }
  if (!is.null(args$seed)) {
    seed <- suppressWarnings(as.integer(args$seed))
    if (is.finite(seed)) set.seed(seed)
  }
  
  if (!is.matrix(Zresidual)) {
    Zresidual <- as.matrix(Zresidual)
  }
  
  sw.pv <- rep(NA_real_, ncol(Zresidual))
  
  for (i in seq_len(ncol(Zresidual))) {
    z <- Zresidual[, i]
    
    # Replace +/-Inf with large finite values tied to column scale
    inf_idx <- is.infinite(z)
    if (any(inf_idx, na.rm = TRUE)) {
      finite_abs_max <- suppressWarnings(max(abs(z[is.finite(z)]), na.rm = TRUE))
      if (!is.finite(finite_abs_max)) finite_abs_max <- 0
      z[inf_idx] <- sign(z[inf_idx]) * (finite_abs_max + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    
    nan_idx <- is.nan(z)
    if (any(nan_idx, na.rm = TRUE)) {
      message("NaNs exist! The model or the fitting process has a problem!")
    }
    
    # Drop NA/NaN and any remaining non-finite
    z <- z[is.finite(z)]
    
    # Shapiro-Wilk requires 3..5000 and not all identical
    if (length(z) < 3L) {
      sw.pv[i] <- NA_real_
      next
    }
    
    if (length(z) > max_n) {
      z <- sample(z, max_n)
    }
    
    sw.pv[i] <- tryCatch(
      stats::shapiro.test(z)$p.value,
      error = function(e) NA_real_
    )
  }
  
  sw.pv
}
