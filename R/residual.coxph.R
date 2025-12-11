#' Residuals for Cox proportional hazards models
#'
#' Compute several types of residuals (censored Z-residuals, Cox–Snell,
#' martingale, and deviance) for Cox proportional hazards models fitted with
#' \code{\link[survival]{coxph}}.
#'
#' @importFrom survival basehaz Surv
#' @importFrom stats qnorm
#'
#' @param fit_coxph A fitted \code{\link[survival]{coxph}} model object
#'   without a shared frailty term. The response should be a right-censored
#'   survival object, typically \code{Surv(time, status)} or
#'   \code{Surv(tstart, tstop, status)}.
#' @param newdata A \code{data.frame} containing the variables required by
#'   \code{fit_coxph$formula}, including the \code{Surv()} response and all
#'   covariates. Residuals are evaluated for the observations in \code{newdata}.
#' @param residual.type Character string specifying the type of residual to
#'   compute. Must be one of \code{"censored Z-residual"},
#'   \code{"Cox-Snell"}, \code{"martingale"}, or \code{"deviance"}. The
#'   default is the full vector \code{c("censored Z-residual", "Cox-Snell",
#'   "martingale", "deviance")}, but in typical use a single value should be
#'   supplied.
#'
#' @return A numeric matrix of dimension \eqn{n \times 1}, where \eqn{n} is
#'   the number of observations in \code{newdata}. The single column is named
#'   according to \code{residual.type}. Several attributes are attached:
#'   \itemize{
#'     \item \code{Survival.Prob}: vector of survival probabilities
#'       \eqn{S_i(t_i)}.
#'     \item \code{linear.pred}: vector of linear predictors \eqn{\eta_i}.
#'     \item \code{covariates}: model matrix of covariates (columns used in
#'       the linear predictor).
#'     \item \code{censored.status}: event indicator (1 = event, 0 = censored).
#'     \item \code{object.model.frame}: the \code{model.frame} constructed
#'       from \code{fit_coxph$formula} and \code{newdata}.
#'   }
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   data(lung)
#'
#'   ## Cox PH model
#'   fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'
#'   ## Censored Z-residuals
#'   r_z <- residual.coxph(fit_cox, newdata = lung,
#'                         residual.type = "censored Z-residual")
#'
#'   ## Cox–Snell residuals
#'   r_cs <- residual.coxph(fit_cox, newdata = lung,
#'                          residual.type = "Cox-Snell")
#'
#'   ## Martingale residuals
#'   r_m <- residual.coxph(fit_cox, newdata = lung,
#'                         residual.type = "martingale")
#'
#'   ## Deviance residuals
#'   r_d <- residual.coxph(fit_cox, newdata = lung,
#'                         residual.type = "deviance")
#' }


residual.coxph<-function (fit_coxph, newdata,
                          residual.type=c("censored Z-residual", "Cox-Snell",
                                          "martingale", "deviance"))
{
  residual.type <- match.arg(residual.type)

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

  resid.class <- switch(residual.type,
                        "Cox-Snell"            = "cs.residual",
                        "deviance"             = "dev.resid",
                        "martingale"           = "martg.resid",
                        "censored Z-residual"  = "cz.resid"
  )

  class(resid.value) <- resid.class

  return(resid.value)

}

