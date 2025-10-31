#' Shapiro-Wilk Normality Test for Z-Residuals
#'
#' Performs the Shapiro-Wilk test for normality on each column of a matrix of Z-residuals.
#'
#' @param Zresidual A numeric matrix of Z-residuals, where each column represents
#'   a separate set of residuals (e.g., from different posterior predictive draws or variables).
#' @param ... Additional arguments.
#'
#' @return A numeric vector of Shapiro-Wilk p-values, one for each column of \code{Zresidual}.
#'
#' @details
#' Infinite or non-finite values are handled by replacing negative infinity with -1e10
#' and positive infinity with 1e10. Any NaN or remaining infinite values are reported via messages.
#'
#' @examples
#' \dontrun{
#' Zres <- matrix(rnorm(100), ncol=5)
#' sw.pvals <- sw.test.zresid(Zres)
#' }
#'
#' @seealso
#' \code{\link[stats]{shapiro.test}}, \code{\link{Zresidual}}
#'
#' @export sw.test.zresid

sw.test.zresid <- function (Zresidual, ...)
{
  sw.pv <- rep(0,ncol(Zresidual))
  for(i in 1:ncol(Zresidual)){
    id.negtv.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] < 0)
    id.pos.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] > 0)
    Zresidual[,i][id.negtv.inf]<- -1e10
    Zresidual[,i][id.pos.inf]<- 1e10

    id.nan <- which(is.nan(Zresidual[,i]))
    id.infinity <- which (is.infinite(Zresidual[,i]))

    if (length(id.infinity) > 0L) message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    if(length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")


    sw.pv[i]<- shapiro.test(Zresidual[,i])$p.value
  }
  return(sw.pv)
}


