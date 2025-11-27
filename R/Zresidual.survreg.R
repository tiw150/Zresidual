#' Z-residuals for Accelerated failure time model models
#'
#' This function calculates Z-residuals based on a fitted `survreg` object from the `survival` package.
#'
#'
#' @importFrom survival psurvreg dsurvreg
#' @importFrom stats qnorm runif
#'
#' @param fit_survreg A fitted object from the `survreg` function in the `survival` package.
#' @param newdata Optional \code{data.frame} containing the variables used in \code{fit_coxph$formula}. If \code{NULL} (the default), residuals are computed on the original data used to fit the model. If supplied, \code{newdata} must contain the survival response and all covariates appearing in the original model formula.
#' @param n.rep An integer specifying the number of random draws to use for calculating the Z-residuals for censored observations. Defaults to `nrep` (which should be defined in your environment, a common choice is 100).
#'
#' @export
#'
#' @return A numeric matrix of dimension \eqn{n \times} \code{n.rep},  where \eqn{n} is the number of observations in the (new) model frame.   Each column corresponds to one set of Z-residuals. The returned matrix has the following attributes attached:
#'   \itemize{
#'     \item \code{Survival.Prob}: vector of survival probabilities \eqn{S_i(t_i)}.
#'     \item \code{linear.pred}: vector of linear predictors \eqn{\eta_i}.
#'     \item \code{covariates}: data frame of covariates (model frame without the response).
#'     \item \code{censored.status}: event indicator (1 = event,0 = censored).
#'     \item \code{object.model.frame}: the \code{model.frame} used to compute the residuals.
#'     \item \code{type}: character string \code{"survival"}.
#'   }
#'
#'
#' @examples
#'
#' library(survival)
#'
#' # Fit a Weibull survival regression model
#' data(cancer, package="survival")
#' fit_weibull <- survreg(Surv(time, status) ~ age, data = lung, dist = "weibull")
#'
#' # Calculate Z-residuals for the fitted model
#' z_residuals <- Zresidual.survreg(fit_weibull,newdata=lung)
#' head(z_residuals)
#'
#'
Zresidual.survreg<-function(fit_survreg,newdata,n.rep=1)
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


