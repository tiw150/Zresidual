#' Z-residuals for parametric survival regression models (survival)
#'
#' @description
#' S3 method for [Zresidual()] when the fitted model is a
#' [survival::survreg()] object (internally tagged as `"survreg.survival"`).
#' This is a thin wrapper around the existing `Zresidual.survreg()` core
#' implementation: it simply passes the fitted object and optional `data`
#' to the core function and then adds the `"zresid"` class.
#'
#' Users are expected to call [Zresidual()] directly rather than calling
#' `Zresidual.survreg.survival()`.
#'
#' @param object A fitted [survival::survreg()] model.
#' @param nrep Integer; number of randomized Z-residual replicates to
#'   generate. Defaults to `1`.
#' @param data Optional data frame containing the survival response and
#'   covariates; if `NULL`, the original model frame is used.
#' @param type Optional character string controlling the residual type,
#'   interpreted by the underlying implementation (if used). For `survreg` models,
#'   this is typically set internally to \code{"survival"}.
#' @param method Character string specifying the residual calculation method
#'   (if applicable to the underlying worker function). Currently unused
#'   by the default implementation.
#' @param ... Further arguments passed to the underlying implementation
#'   functions. Currently unused.
#' @return
#' A numeric matrix of dimension \eqn{n \times} \code{nrep}, with
#' additional attributes as produced by `Zresidual_survreg_survival()`. The
#' returned object is given class `"zresid"` in addition to any
#' existing classes.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'   fit_surv <- survreg(Surv(time, status) ~ age + sex,
#'                       data = lung, dist = "weibull")
#'   z_surv <- Zresidual(fit_surv, nrep = 5, data = lung)
#' }
#'
#' @method Zresidual survreg.survival
#' @export
Zresidual.survreg.survival <- function(object,nrep = 1,data = NULL,type = "survival", method=NULL,
                                       ...) {

  out <- Zresidual_survreg_survival(
    fit_survreg = object,
    newdata     = data,
    n.rep       = nrep,
    ...
  )

  class(out) <- c("zresid", class(out))
  out
}


#' @keywords internal
Zresidual_survreg_survival<-function(fit_survreg,newdata,n.rep=1, ...)
{
  if(is.null(newdata)){
    mf<-model.frame.survreg(fit_survreg)
    mf_nc<- ncol(mf)
    fix_var<-model.matrix.survreg(fit_survreg)
  }

  if(!is.null(newdata)){
    mf <- model.frame(fit_survreg$terms, newdata)
    mf_nc<-ncol (mf)
    fix_var<-model.matrix(fit_survreg$terms, newdata,drop=FALSE)
  }

  y<-  mf[[1]]
  distr<-fit_survreg$dist
  parms<- as.numeric(fit_survreg[["parms"]])
  alpha_hat<-1/fit_survreg$scale

  if (distr %in% c("weibull","exponential","logistic","lognormal",
                   "loglogistic","gaussian", "loggaussian","rayleigh"))
  {
    SP<-1-(psurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%fit_survreg$coefficients),
                    scale=1/alpha_hat, distribution=distr))
    haz <- dsurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%fit_survreg$coefficients),
                    scale=1/alpha_hat, distribution=distr)/SP
  }else if  (distr %in% "t")
  {
    SP<-1-(psurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%fit_survreg$coefficients),
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms))
    haz <- dsurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%fit_survreg$coefficients),
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms)/SP
  }else stop ("The distribution is not supported")
  censored <- which(as.matrix(y)[,-1]==0)
  n.censored <- length(censored)
  # Z-residual
  Zresid<-matrix(0,length(SP),n.rep)
  col_name<-rep(0,n.rep)
  for(i in 1:n.rep){
    RSP <- SP
    RSP[censored] <- RSP[censored]*runif(n.censored)
    Zresid[,i] <- -qnorm(RSP)
    col_name[i]<-paste("Z-residual ",i,sep = "")
  }
  colnames(Zresid)<- col_name
  #####
  censored.status<- (as.matrix(y)[,-1])
  lp<-fix_var %*%fit_survreg$coefficients

  Zresid.value<-as.matrix(Zresid)

  attributes(Zresid.value) <- c(attributes(Zresid.value), list(
    Survival.Prob= SP,
    linear.pred = lp,
    covariates = mf[,-1],
    censored.status= censored.status,
    object.model.frame=mf,
    type = "survival"

  ))
  return(Zresid.value)

}


# Zresidual.survreg<-function(fit_survreg)
# {
#   distr<-fit_survreg$dist
#   y<- fit_survreg$y
#   m <- nrow (y)
#   parms<- as.numeric(fit_survreg[["parms"]])
#   alpha_hat<-1/fit_survreg$scale
#
#   if (distr %in% c("weibull","exponential","logistic","lognormal",
#                    "loglogistic","gaussian", "loggaussian","rayleigh"))
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr)/SP
#   }else if  (distr %in% "t")
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms)/SP
#   }else stop ("The distribution is not supported")
#   censored <- which(as.data.frame(as.matrix(y))[,-1]==0)
#   n.censored <- length(censored)
#   # Z-residual
#   RSP <- SP
#   RSP[censored] <- RSP[censored]*runif(n.censored)
#   Zresid <- -qnorm(RSP)
#   #Normalized unmodified SPs
#   USP<-SP
#   USP[USP==1] <- .999999999
#   censored.Zresid<- -qnorm(USP)
#   # Unmodified CS residual
#   ucs<- -log(SP)
#   # Modified CS residual
#   MSP<- SP
#   MSP[censored] <- SP[censored]/exp(1)
#   mcs <- -log(MSP)
#   #Martingale Residual
#   martg<- as.data.frame(as.matrix(fit_survreg$y))[,-1]- ucs
#   #Deviance Residual
#   dev<- sign(martg)* sqrt((-2)*(martg+as.data.frame(as.matrix(fit_survreg$y))[,-1]*
#                                   log(as.data.frame(as.matrix(fit_survreg$y))[,-1]-martg)))
#
#   list(Zresid=Zresid,censored.Zresid=censored.Zresid,SP=SP,
#        ucs=ucs, mcs=mcs,martg=martg, dev=dev)
# }
#
# Zresidual.survreg1<-function(fit_survreg)
# {
#   distr<-fit_survreg$dist
#   y<- fit_survreg$y
#   m <- nrow (y)
#   parms<- as.numeric(fit_survreg[["parms"]])
#   alpha_hat<-1/fit_survreg$scale
#
#   if (distr %in% c("weibull","exponential","logistic","lognormal",
#                    "loglogistic","gaussian", "loggaussian","rayleigh"))
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr)/SP
#   }else if  (distr %in% "t")
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     fit_survreg[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms)/SP
#   }else stop ("The distribution is not supported")
#   censored <- which(as.data.frame(as.matrix(y))[,-1]==0)
#   n.censored <- length(censored)
#   # Z-residual
#   RSP <- SP
#   RSP[censored] <- RSP[censored]*runif(n.censored)
#   Zresid <- -qnorm(RSP)
#   #Normalized unmodified SPs
#   USP<-SP
#   USP[USP==1] <- .999999999
#   censored.Zresid<- -qnorm(USP)
#   # Unmodified CS residual
#   ucs<- -log(SP)
#   # Modified CS residual
#   MSP<- SP
#   MSP[censored] <- SP[censored]/exp(1)
#   mcs <- -log(MSP)
#   #Martingale Residual
#   martg<- as.data.frame(as.matrix(fit_survreg$y))[,-1]- ucs
#   #Deviance Residual
#   dev<- sign(martg)* sqrt((-2)*(martg+as.data.frame(as.matrix(fit_survreg$y))[,-1]*
#                                   log(as.data.frame(as.matrix(fit_survreg$y))[,-1]-martg)))
#
#   list(Zresid=Zresid,censored.Zresid=censored.Zresid,SP=SP,
#        ucs=ucs, mcs=mcs,martg=martg, dev=dev)
# }


