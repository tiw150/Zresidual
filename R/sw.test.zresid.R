#' Shapiro-Wilk normality test for Z-residuals
#'
#' @description
#' Applies the Shapiro-Wilk test column-wise to a matrix of Z-residuals.
#'
#' @param Zresidual A numeric vector or matrix of Z-residuals.
#' @param ... Optional named arguments. Supported arguments include
#'   \code{max_n} and \code{seed}.
#'
#' @return A numeric vector of p-values, one per column.
#'
#' @examples
#' set.seed(1)
#' z <- matrix(rnorm(60), ncol = 2)
#' sw.test.zresid(z)
#'
#' @export
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
