#'Z-residuals for Cox proportional hazards models
#'
#' This function calculates Z-residuals based on a fitted `coxph` object from the `survival` package.
#'
#' @importFrom survival
#' @importFrom stats qnorm runif
#'
#' @param fit_coxph A fitted object from the `coxph` function in the `survival` package.
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
#' \dontrun{
#'   library(survival)
#'   fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   # Calculate Z-residuals for the fitted model
#'   z <- Zresidual.coxph(fit)
#'   # multiple randomized replicates
#'   z_multi <- Zresidual.coxph(fit, n.rep = 10)
#'
#' # Calculate Z-residuals for new data
#' new_data <- data.frame(x = c(1, 2, 0.5))
#' Z_resid_new <- Zresidual.survreg(fit_weibull, newdata = new_data)
#' head(z_residuals_new)
#'}

Zresidual.coxph<-function (fit_coxph, newdata,n.rep=1)
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
 #     type = "coxph",
      Survival.Prob= SP,
      linear.pred = lp.new,
      covariates = mf_new[,-1],
      censored.status= censored.status,
      object.model.frame=mf_new,
      type = "survival"

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


