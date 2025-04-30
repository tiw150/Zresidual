#' A function to calculate ANOVA of Zresidual
#'
#' @param Zresidual A z-residual.
#' @param X Linear predictor or covariate
#' @param k.anova Number of bins if applicable
#' @export

aov.test.zresid <- function (Zresidual,X = c("lp", "covariate"), k.anova=10)
{
  if (missing(X)) X = "lp"
  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")
    id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
    id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
    Zresidual[id.negtv.inf]<- -1e10
    Zresidual[id.pos.inf]<- 1e10
    aov.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      aov.pv[j]<- test.nl.aov(Zresidual[,j], fitted.value, k.anova)
    }
    aov.pv
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
    } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". "))}

    id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
    id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
    Zresidual[id.negtv.inf]<- -1e10
    Zresidual[id.pos.inf]<- 1e10
    aov.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      aov.pv[j]<- test.nl.aov(Zresidual[,j], fitted.value[,i], k.anova)
    }
    aov.pv
  }
  return(aov.pv)
}

