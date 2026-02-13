#' A function to calculate Shapiro-Francia test of Zresidual
#' @importFrom nortest sf.test
#' @param Zresidual A matrix (n x p) of Z-residuals.
#' @export sf_test.zresid
sf_test.zresid <- function(Zresidual) {
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
