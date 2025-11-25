#' Residuals for shared frailty Cox proportional hazards models
#'
#' Compute several types of residuals (censored Z-residuals, Cox–Snell,
#' martingale, and deviance) for shared frailty Cox proportional hazards
#' models fitted with \code{\link[survival]{coxph}} using a multiplicative
#' frailty term (e.g., \code{frailty(group)}).
#'
#'
#' The function is designed for an out-of-sample / cross-validation setting:
#' \itemize{
#'   \item \code{traindata} is used to reconstruct the baseline cumulative
#'         hazard and the relationship between covariates and shared frailty.
#'   \item \code{newdata} is the dataset for which residuals are computed.
#' }
#' Both datasets must contain the survival response, all fixed-effect
#' covariates, and the same frailty grouping factor (with compatible factor
#' levels). The grouping variable is extracted from the \code{frailty()}
#' term in \code{fit_coxph$formula} and must be a factor in both
#' \code{traindata} and \code{newdata}. The internal implementation treats
#' models with many groups (\code{gpnumber > 5}) and few groups
#' (\code{gpnumber <= 5}) slightly differently, reflecting how frailty
#' coefficients are stored in the fitted \code{coxph} object.
#'
#' @importFrom survival Surv
#'
#' @param fit_coxph A fitted \code{\link[survival]{coxph}} object for a shared
#'   frailty Cox model, typically specified with a term such as
#'   \code{frailty(group)} in the model formula. The object must contain
#'   cluster-level frailty estimates in \code{fit_coxph$frail}.
#' @param traindata A \code{data.frame} used to reconstruct the baseline
#'   hazard and frailty effects. It must contain the survival response, all
#'   fixed-effect covariates, and the shared frailty grouping factor appearing
#'   in \code{fit_coxph$formula}.
#' @param newdata A \code{data.frame} for which residuals are to be computed.
#'   It must contain the survival response, covariates, and the shared frailty
#'   grouping factor used in the original model. The factor levels of the
#'   grouping variable must be compatible with those in \code{traindata}.
#' @param residual.type Character string specifying the type of residual to
#'   compute. Must be one of \code{"censored Z-residual"}, \code{"Cox-Snell"},
#'   \code{"martingale"}, or \code{"deviance"}. The default is the full vector
#'   \code{c("censored Z-residual", "Cox-Snell", "martingale", "deviance")},
#'   but in typical use a single value should be supplied.
#'
#' @return A numeric matrix of dimension \eqn{n \times 1}, where \eqn{n} is
#'   the number of observations in \code{newdata}. The single column is named
#'   according to \code{residual.type}. Several attributes are attached:
#'   \itemize{
#'     \item \code{Survival.Prob}: vector of survival probabilities
#'       \eqn{S_{ij}(t_i)} for each observation in \code{newdata}.
#'     \item \code{linear.pred}: vector of fixed-effect linear predictors
#'       \eqn{\eta_{ij}} (excluding the frailty term).
#'     \item \code{covariates}: model matrix of fixed-effect covariates
#'       used in the linear predictor.
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
#'   lung$inst <- factor(lung$inst)
#'
#'   ## Shared frailty Cox model
#'   set.seed(1)
#'   idx <- sample(seq_len(nrow(lung)), size = floor(0.7 * nrow(lung)))
#'   train_dat <- lung[idx, ]
#'   test_dat  <- lung[-idx, ]
#'
#'   fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
#'                      data = train_dat)
#'
#'   ## Censored Z-residuals on the test set
#'   r_z <- residual.coxph.frailty(fit_frail,
#'                                 traindata = train_dat,
#'                                 newdata   = test_dat,
#'                                 residual.type = "censored Z-residual")
#'
#'   ## Cox–Snell residuals
#'   r_cs <- residual.coxph.frailty(fit_frail,
#'                                  traindata = train_dat,
#'                                  newdata   = test_dat,
#'                                  residual.type = "Cox-Snell")
#' }

residual.coxph.frailty <- function (fit_coxph, traindata, newdata,
                                    residual.type=c("censored Z-residual", "Cox-Snell",
                                                    "martingale", "deviance"))
{
  form<-(fit_coxph$formula)[[3]]
  group_id_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[3]
  if(!is.factor(traindata[[group_id_name]])) stop("The group ID must be factor!")
  if(!is.factor(newdata[[group_id_name]])) stop("The group ID must be factor!")
  gpnumber<-length(levels(traindata[[group_id_name]]))

  mf <- model.frame(fit_coxph$formula, traindata)
  mf_nc<-ncol (mf)
  mm <- model.matrix(fit_coxph$formula, traindata)
  mm_nc <- ncol (mm)
  if(gpnumber>5){
    fix_var<-mm[,-c(1,mm_nc),drop=FALSE]
    explp<-exp(fix_var %*% fit_coxph$coefficients+fit_coxph$frail[mf[,mf_nc]])
    z_hat<-exp(fit_coxph$frail[mf[,mf_nc]])
  }
  if(gpnumber<=5){
    fix_var<-mm[,-c(1,mf_nc:mm_nc),drop=FALSE]
    coef_number<- ncol(fix_var)
    frailty<-fit_coxph$coefficients[coef_number+1:gpnumber]
    explp<-exp(fix_var %*% fit_coxph$coefficients[1:coef_number]+
                 frailty[mf[,mf_nc]])
    z_hat<-exp(frailty[mf[,mf_nc]])
  }

  Y <- mf[[1]]
  if(!inherits(Y, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y) != 3) {
    # making it all in (tstart, tstop) format
    Y <- Surv(rep(0, nrow(Y)), Y[,1], Y[,2])
  }

  # this one gives the baseline cumulative hazard at all the time points;
  getchz <- function(Y, explp) {
    death <- (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
    dtime <- Y[, ncol(Y) - 1] # this is the tstop

    time <- sort(unique(dtime)) # unique tstops

    nevent <- as.vector(rowsum(1 * death, dtime))

    nrisk <- rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
    delta <- min(diff(time))/2
    etime <- c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times

    indx <- approx(etime, 1:length(etime), time, method = "constant",
                   rule = 2, f = 1)$y

    esum <- rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
    nrisk <- nrisk - c(esum, 0)[indx]

    haz <- nevent/nrisk
    cumhaz <- cumsum(haz)

    out<-data.frame(time = time, haz = haz, cumhaz = cumhaz)
    return (out)
  }
  df<-getchz(Y=Y,explp=explp)
  t<-df$time
  haz<-df$haz
  H_0<-df$cumhaz
  f <- stepfun(t[-1], H_0)

  full_data<-rbind(newdata,traindata)
  mf_new<-model.frame(fit_coxph$formula,full_data)[c(1:nrow(newdata)),]
  mf_nc_new<- ncol (mf_new)
  mm_new <-model.matrix(fit_coxph$formula,full_data)[c(1:nrow(newdata)),,drop=FALSE]
  mm_nc_new <- ncol (mm_new)

  if(gpnumber>5){
    fix_var_new<-mm_new[,-c(1,mm_nc_new),drop=FALSE]
    z_hat_new<-exp(fit_coxph$frail[mf_new[,mf_nc_new]])
    lp_new<-fix_var_new %*% fit_coxph$coefficients+fit_coxph$frail[mf_new[,mf_nc_new]]
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients)
  }
  if(gpnumber<=5){
    fix_var_new<-mm_new[,-c(1,mf_nc_new:mm_nc_new),drop=FALSE]
    coef_number_new<- ncol(fix_var_new)
    frailty_new<-fit_coxph$coefficients[coef_number_new+1:gpnumber]
    lp_new<-fix_var_new %*% fit_coxph$coefficients[1:coef_number_new]+frailty_new[mf_new[,mf_nc_new]]
    explp_new<-exp(fix_var_new %*% fit_coxph$coefficients[1:coef_number_new])
    z_hat_new<-as.numeric(exp(frailty_new[mf_new[,mf_nc_new]]))
  }

  Y_new <- mf_new[[1]]
  if(!inherits(Y_new, "Surv")) stop("left hand side not a survival object")
  if(ncol(Y_new) != 3) {
    # making it all in (tstart, tstop) format
    Y_new <- Surv(rep(0, nrow(Y_new)), Y_new[,1], Y_new[,2])
  }
  H0_new<-f(Y_new[,2])
  #Survival Function
  SP<- exp(-z_hat_new*as.vector(explp_new)*H0_new)
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

  resid.value<-as.matrix(resid)
  colnames(resid.value)[1] <- residual.type

  attributes(resid.value) <- c(attributes(resid.value), list(
    Survival.Prob= SP,
    linear.pred = lp_new,
    covariates = fix_var_new,
    censored.status= censored.status,
    object.model.frame=mf_new

  ))
  return(resid.value)
}





