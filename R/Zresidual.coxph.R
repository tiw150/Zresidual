#input: coxfit_fit is a coxph object
Zresidual.coxph<-function (fit_coxph, newdata,n.rep=nrep)
{
  if(is.null(newdata)){
    mf_new<-model.frame.coxph(fit_coxph)
    mf_nc_new<- ncol (mf_new)
    mm_new<-model.matrix.coxph(fit_coxph)
    fix_var_new<-mm_new
  }

  if(!is.null(newdata)){
    mf_new<-model.frame(fit_coxph$formula,newdata)
    mf_nc_new<- ncol (mf_new)
    mm_new<- model.matrix.coxph(fit_coxph,newdata)
    fix_var_new<-mm_new[, ,drop=FALSE]
  }

  basecumhaz<-basehaz(fit_coxph,centered = F)
  t<-basecumhaz$time
  H_0<-basecumhaz$hazard
  f <- stepfun(t[-1], H_0)

  explp_new<-exp(fix_var_new %*% fit_coxph$coefficients)
  Y_new <- mf_new[[1]]
  if(!inherits(Y_new, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y_new) != 3) {
    # making it all in (tstart, tstop) format
    Y_new <- Surv(rep(0, nrow(Y_new)), Y_new[,1], Y_new[,2])
  }
  H0_new<-f(Y_new[,2])

  #Survival Function
  SP<- exp(-as.vector(explp_new)*H0_new)
  censored <- which(Y_new[,3]==0)
  n.censored <- length(censored)
  #Z-residual
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
  censored.status<- Y_new[,3]
  lp.new<-fix_var_new %*% fit_coxph$coefficients

  Zresid.value<-as.matrix(Zresid)

#  class(Zresid.value) <- c("zresid", class(Zresid.value))

  attributes(Zresid.value) <- c(attributes(Zresid.value), list(
      Survival.Prob= SP,
      linear.pred = lp.new,
      covariates = mf_new[,-1],
      censored.status= censored.status,
      object.model.frame=mf_new

  ))
  return(Zresid.value)
}


#input: coxfit_fit is a coxph object
# Zresidual.coxph<-function(coxfit_fit)
# {
#   y<- coxfit_fit$y
#   m <- nrow (y)
#   mre <- resid(coxfit_fit, type="martingale")
#   dre <- resid(coxfit_fit, type="deviance")
#   #Unmodified Cox-Snell residual
#   ucs <- as.data.frame(as.matrix(y))[,-1] - mre
#   #Survival Function
#   SP<- exp(-ucs)
#   censored <- which(as.data.frame(as.matrix(y))[,-1]==0)
#   n.censored <- length(censored)
#   #NMSP residual (normal-transformed modified survival prob)
#   MSP<- SP
#   MSP[censored] <- SP[censored]/exp(1)
#   mcs <- -log(MSP)
#   nmsp<- qnorm(MSP)
#   # Z-residual
#   RSP <- SP
#   RSP[censored] <- RSP[censored]*runif(n.censored)
#   Zresid <- -qnorm(RSP)
#   Zresid.sw.pvalue<-shapiro.test(Zresid)$p.value
#   list(Zresid=Zresid,RSP=RSP,Zresid.sw.pvalue=Zresid.sw.pvalue,
#        ucs=ucs,mcs=mcs)
# }
# outputs:
## RSP --- Randomized Survival Probabilities
## Zresid --- Z-residual
## Zresid.sw.pvalue --- GOF test p-values by applying SW test to Z-residual
## ucs --- unmodified CS residuals
## mcs --- modified CS residuals
## martg --- Martingale residuals
## dev --- Deviance residuals
## haz_fn --- hazard function of cs residuals


