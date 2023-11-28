#input: coxfit_fit is a coxph object
residual.coxph<-function (fit_coxph, newdata,
                          residual.type=c("censored Z-residual", "Cox-Snell",
                                          "martingale", "deviance"))
{
  basecumhaz<-basehaz(fit_coxph,centered = F)
  t<-basecumhaz$time
  H_0<-basecumhaz$hazard
  f <- stepfun(t[-1], H_0)

  mf_new<-model.frame(fit_coxph$formula,newdata)
  mf_nc_new<- ncol (mf_new)
  mm_new<-model.matrix(fit_coxph$formula,newdata)
  fix_var_new<-mm_new[,-1,drop=FALSE]

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
    resid<- Y_new[,3] - ucs
  }

  if (residual.type == "deviance"){
    #Deviance Residual
    ucs<- -log(SP)
    martg<-Y_new[,3] - ucs
    resid<- sign(martg)* sqrt((-2)*(martg+Y_new[,3]*log(Y_new[,3]-martg)))
  }

  censored.status<- (as.matrix(Y_new)[,-1])[,2]
  lp.new<-fix_var_new %*% fit_coxph$coefficients

  resid.value<-as.matrix(resid)
  colnames(resid.value)[1] <- residual.type

  attributes(resid.value) <- c(attributes(resid.value), list(
    Survival.Prob= SP,
    linear.pred = lp.new,
    covariates = fix_var_new,
    censored.status= censored.status,
    object.model.frame=mf_new

  ))
  return(resid.value)

}

