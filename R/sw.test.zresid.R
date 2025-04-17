#' A function to calculate Shapiro Wilk test of Zresidual
#'
#' @param Zresidual A Z-residual.
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