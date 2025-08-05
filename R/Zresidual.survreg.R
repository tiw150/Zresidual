#' Calculate Z-residuals for a fitted survival regression model
#'
#' This function calculates Z-residuals based on a fitted `survreg` object from the `survival` package.
#'
#'
#' @importFrom survival psurvreg dsurvreg
#' @importFrom stats qnorm runif
#'
#' @param fit_survreg A fitted object from the `survreg` function in the `survival` package.
#' @param newdata An optional data frame containing new observations for which to calculate the Z-residuals. If `NULL` (default), residuals are calculated for the data used to fit the model.
#' @param n.rep An integer specifying the number of random draws to use for calculating the Z-residuals for censored observations. Defaults to `nrep` (which should be defined in your environment, a common choice is 100).
#'
#' @export
#'
#' @return  \itemize{ A matrix of Z-residuals. Each column represents a set of Z-residuals based on different random draws for censored observations. The matrix has the following attributes:
#'   \item{Survival.Prob:} {The estimated survival probabilities.}
#'   \item{linear.pred:} {The linear predictors from the survival regression model.}
#'   \item{covariates:} {The covariate values used in the model.}
#'   \item{censored.status:} {The censoring status (0 for censored, 1 for event).}
#'   \item{object.model.frame:} {The model frame used for the analysis.}
#'}
#'
#' @examples
#'
#' library(survival)
#'
#' # Fit a Weibull survival regression model
#' fit_weibull <- survreg(Surv(time, status) ~ x, data = lung, dist = "weibull")
#'
#' # Calculate Z-residuals for the fitted model
#' z_residuals <- Zresidual.survreg(fit_weibull)
#' head(z_residuals)
#'
#' # Calculate Z-residuals for new data
#' new_data <- data.frame(x = c(1, 2, 0.5))
#' Z_residuals_new <- Zresidual.survreg(fit_weibull, newdata = new_data)
#' head(z_residuals_new)
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
    object.model.frame=mf

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


