#' A function to calculate Bartlett of Zresidual
#'
#' @param Zresidual A Z-residual.
#' @param X Linear predictor or covariate
#' @param k.bl Number of bins if applicable
#' @export
#'

bartlett.test.zresid <- function (Zresidual, X = c("lp", "covariate"), k.bl=10)
{
  if (missing(X))
    X = "lp"

  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10

  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")
    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      bl.pv[j]<-test.var.bartl(Zresidual[,j], fitted.value, k.bl)
    }
    bl.pv
  }
  if (X != "lp") {
    fitted.value <- attr(Zresidual, "covariates")
    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". ")) }

    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      id.na <- which(is.na(Zresidual[,j]))
      count.id <- which(!is.na(Zresidual[,j]))
      new.Zresidual <- Zresidual[count.id, j]
      #if(length(id.na) > 0) message("NAs omitted.")
      bl.pv[j]<-test.var.bartl(new.Zresidual, fitted.value[,i][count.id], k.bl)
    }
    bl.pv
  }
  bl.pv
}
