#input: survreg_fit is a survreg object
residual.survreg<-function(survreg_fit,newdata,
                           residual.type=c("censored Z-residual", "Cox-Snell",
                                           "martingale", "deviance"))
{
  distr<-survreg_fit$dist
  parms<- as.numeric(survreg_fit[["parms"]])
  alpha_hat<-1/survreg_fit$scale

  mf <- model.frame(survreg_fit$terms, newdata)
  mf_nc<-ncol (mf)
  y<-  mf[[1]]
  fix_var<-model.matrix(survreg_fit$terms, newdata,drop=FALSE)

  if (distr %in% c("weibull","exponential","logistic","lognormal",
                   "loglogistic","gaussian", "loggaussian","rayleigh"))
  {
    SP<-1-(psurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%survreg_fit$coefficients),
                    scale=1/alpha_hat, distribution=distr))
    haz <- dsurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%survreg_fit$coefficients),
                    scale=1/alpha_hat, distribution=distr)/SP
  }else if  (distr %in% "t")
  {
    SP<-1-(psurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%survreg_fit$coefficients),
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms))
    haz <- dsurvreg(as.matrix(y)[,1],
                    mean =as.vector(fix_var %*%survreg_fit$coefficients),
                    scale=1/alpha_hat, distribution=distr,
                    parms = parms)/SP
  }else stop ("The distribution is not supported")
  censored <- which(as.matrix(y)[,-1]==0)
  n.censored <- length(censored)

  if (residual.type == "censored Z-residual"){
    #Normalized unmodified SPs
    USP<-SP
    USP[USP==1] <- .999999999
    resid<- -qnorm(USP)
  }
  if (residual.type == "Cox-Snell"){
    # Unmodified CS residual
    resid<- -log(SP)
  }

  if (residual.type == "martingale"){
    #Martingale Residual
    ucs<- -log(SP)
    resid<- y[,3] - ucs
  }

  if (residual.type == "deviance"){
    #Deviance Residual
    ucs<- -log(SP)
    martg<-y[,3] - ucs
    resid<- sign(martg)* sqrt((-2)*(martg+y[,3]*log(y[,3]-martg)))
  }

  censored.status<- (as.matrix(y)[,-1])[,2]
  lp<-fix_var %*%survreg_fit$coefficients

  resid.value<-as.matrix(resid)
  colnames(resid.value)[1] <- residual.type

  attributes(resid.value) <- c(attributes(resid.value), list(
    Survival.Prob= SP,
    linear.pred = lp,
    covariates = fix_var,
    censored.status= censored.status,
    object.model.frame=mf

  ))
  return(resid.value)


}


# Zresidual.survreg1<-function(survreg_fit)
# {
#   distr<-survreg_fit$dist
#   y<- survreg_fit$y
#   m <- nrow (y)
#   parms<- as.numeric(survreg_fit[["parms"]])
#   alpha_hat<-1/survreg_fit$scale
#
#   if (distr %in% c("weibull","exponential","logistic","lognormal",
#                    "loglogistic","gaussian", "loggaussian","rayleigh"))
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr)/SP
#   }else if  (distr %in% "t")
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
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
#   martg<- as.data.frame(as.matrix(survreg_fit$y))[,-1]- ucs
#   #Deviance Residual
#   dev<- sign(martg)* sqrt((-2)*(martg+as.data.frame(as.matrix(survreg_fit$y))[,-1]*
#                                   log(as.data.frame(as.matrix(survreg_fit$y))[,-1]-martg)))
#
#   list(Zresid=Zresid,censored.Zresid=censored.Zresid,SP=SP,
#        ucs=ucs, mcs=mcs,martg=martg, dev=dev)
# }
#
# Zresidual.survreg1<-function(survreg_fit)
# {
#   distr<-survreg_fit$dist
#   y<- survreg_fit$y
#   m <- nrow (y)
#   parms<- as.numeric(survreg_fit[["parms"]])
#   alpha_hat<-1/survreg_fit$scale
#
#   if (distr %in% c("weibull","exponential","logistic","lognormal",
#                    "loglogistic","gaussian", "loggaussian","rayleigh"))
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr)/SP
#   }else if  (distr %in% "t")
#   {
#     SP<-1-(psurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
#                     scale=1/alpha_hat, distribution=distr,
#                     parms = parms))
#     haz <- dsurvreg(as.data.frame(as.matrix(y))[,-2],
#                     survreg_fit[["linear.predictors"]],
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
#   martg<- as.data.frame(as.matrix(survreg_fit$y))[,-1]- ucs
#   #Deviance Residual
#   dev<- sign(martg)* sqrt((-2)*(martg+as.data.frame(as.matrix(survreg_fit$y))[,-1]*
#                                   log(as.data.frame(as.matrix(survreg_fit$y))[,-1]-martg)))
#
#   list(Zresid=Zresid,censored.Zresid=censored.Zresid,SP=SP,
#        ucs=ucs, mcs=mcs,martg=martg, dev=dev)
# }


