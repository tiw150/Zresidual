#' Z-residuals for Cox proportional hazards models (survival package)
#'
#' @description
#' `Zresidual.coxph.survival()` computes randomized Z-residuals for Cox
#' proportional hazards models fitted with [survival::coxph()], with or
#' without shared frailty terms. It is the S3 method for [Zresidual()]
#' when the fitted object is internally tagged with the class
#' `"coxph.survival"`, and is normally called via [Zresidual()] rather
#' than directly.
#'
#' The function automatically detects whether the model formula contains
#' a frailty term (e.g. `frailty(group)`) and dispatches to separate
#' internal implementations for standard Cox models and shared frailty
#' Cox models. In both cases, Z-residuals are computed for a single
#' data set (typically the data used to fit the model) and are intended
#' for in-sample diagnostics.
#'
#' @details
#' There are two main usage patterns:
#'
#' \itemize{
#'   \item \strong{Standard Cox models (no frailty term)}:
#'
#'   When the model formula does not contain a frailty term, the method
#'   computes Z-residuals for a fitted `coxph` object using the survival
#'   response and covariates in \code{data}. If \code{data} is \code{NULL},
#'   the model frame is reconstructed from the fitted object and residuals
#'   are computed on the original data used to fit the model.
#'
#'   Internally, this branch delegates to an implementation function such
#'   as \code{Zresidual_coxph_survival()}.
#'
#'   \item \strong{Shared frailty Cox models}:
#'
#'   If the model includes a multiplicative frailty term (e.g.
#'   \code{frailty(group)}), the method computes randomized Z-residuals
#'   that account for the cluster-level frailty. The same data set
#'   (either \code{data} if supplied, or the original model frame) is used
#'   both to reconstruct the baseline hazard / frailty effects and to
#'   evaluate residuals.
#'
#'   Internally, this branch calls a dedicated frailty implementation,
#'   such as \code{Zresidual_coxph_frailty_survival()}.
#' }
#'
#' @importFrom survival Surv
#'
#' @param object A fitted [survival::coxph()] model. The function supports
#'  both standard Cox models and shared frailty Cox models specified with
#'  a term such as \code{frailty(group)} in the formula.
#' @param data Optional \code{data.frame} containing the survival response
#'  and covariates used in \code{object$terms}. When \code{NULL} (default),
#'  the model frame is reconstructed from \code{object} and residuals are
#'  computed on the original data.
#' @param nrep Integer; number of independent randomized Z-residual
#'  replicates to generate. Defaults to \code{1}. Each replicate corresponds
#'  to a different randomization of censored observations.
#' @param type Optional character string controlling the residual type,
#'   interpreted by the underlying implementation (if used). For Cox models,
#'   this is typically set internally to \code{"survival"}.
#' @param method Character string specifying the residual calculation method
#'   (if applicable to the underlying worker function). Currently unused
#'   by the Cox PH implementation.
#' @param ... Further arguments passed to the underlying implementation
#'  functions. Currently unused.
#' @return
#' A numeric matrix of dimension \eqn{n \times} \code{nrep}, where \eqn{n} is
#' the number of observations in the data set on which residuals are
#' evaluated (either \code{data} if supplied, or the original model frame).
#' Each column corresponds to one set of Z-residuals. The returned matrix has
#' several attributes attached, including:
#' \itemize{
#'   \item \code{Survival.Prob}: vector of survival probabilities
#'     \eqn{S_i(t_i)} (or \eqn{S_{ij}(t_i)} in the frailty case).
#'   \item \code{linear.pred}: vector of linear predictors \eqn{\eta_i}
#'     (fixed effects; excluding the frailty term in shared frailty models).
#'   \item \code{covariates}: data frame of covariates (model frame without
#'     the survival response and grouping factor).
#'   \item \code{censored.status}: event indicator (1 = event, 0 = censored).
#'   \item \code{object.model.frame}: the \code{model.frame} used to compute
#'     the residuals.
#'   \item \code{type}: character string, typically \code{"survival"}.
#' }
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   ## Standard Cox model (no frailty term)
#'   fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   z_cox <- Zresidual(fit_cox, nrep = 10, data = lung)
#'
#'   ## Shared frailty Cox model (in-sample residuals)
#'   lung$inst <- factor(lung$inst)
#'   fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst),
#'                      data = lung)
#'   z_in <- Zresidual(fit_frail, nrep = 5)
#' }
#'
#' @seealso
#' Zresidual(),
#' Zresidual_coxph_survival(),
#' Zresidual_coxph_frailty_survival()
#'
#' @method Zresidual coxph.survival
#' @export
Zresidual.coxph.survival <- function(object,
                                     nrep = 1,
                                     data = NULL,
                                     type   = NULL,
                                     method = NULL,...) {

  frailty_terms <- attr(object$terms, "specials")$frailty

  if (!is.null(frailty_terms)) {
    out <- Zresidual_coxph_frailty_survival(
      fit_coxph = object,
      traindata = data,
      newdata   = data,
      n.rep     = nrep, ...
    )
  } else {
    out <- Zresidual_coxph_survival(
      fit_coxph = object,
      newdata   = data,
      n.rep     = nrep, ...
    )
  }

  class(out) <- c("zresid", class(out))
  out
}


#' @keywords internal
Zresidual_coxph_survival<-function (fit_coxph, newdata,n.rep=1, ...)
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


#' @keywords internal
Zresidual_coxph_frailty_survival <-
  function (fit_coxph, traindata, newdata, n.rep = 1, ...)
  {
    if (is.null(traindata) & is.null(newdata)) {
      mf <- model.frame.coxph(fit_coxph)
      mf_nc <- ncol (mf)
      groupid <- as.factor(mf[, ncol(mf)])
      if (!is.factor(groupid))
        stop("The group ID must be factor!")
      gpnumber <- length(levels(groupid))
      mm <- model.matrix.coxph(fit_coxph)
      mm_nc <- ncol (mm)

      if (gpnumber > 5) {
        fix_var <- mm[, -mm_nc, drop = FALSE]
        lp <- fix_var %*% fit_coxph$coefficients
        explp <- exp(fix_var %*% fit_coxph$coefficients +
                       fit_coxph$frail[as.vector(mm[, mm_nc])])
        z_hat <- exp(fit_coxph$frail[as.vector(mm[, mm_nc])])
      }
      if (gpnumber <= 5) {
        fix_var <- mm[, -c((mm_nc - gpnumber + 1):mm_nc), drop = FALSE]
        coef_number <- ncol(fix_var)
        frailty <- fit_coxph$coefficients[coef_number + 1:gpnumber]
        lp <- fix_var %*% fit_coxph$coefficients[1:coef_number]
        explp <- exp(fix_var %*% fit_coxph$coefficients[1:coef_number] +
                       frailty[mf[, mf_nc]])
        z_hat <- exp(frailty[mf[, mf_nc]])
      }
      Y <- mf[[1]]
      if (!inherits(Y, "Surv"))
        stop("left hand side not a survival object")
      if (ncol(Y) != 3) {
        # making it all in (tstart, tstop) format
        Y <- Surv(rep(0, nrow(Y)), Y[, 1], Y[, 2])
      }
      # this one gives the baseline cumulative hazard at all the time points;
      getchz <- function(Y, explp) {
        death <-
          (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
        dtime <- Y[, ncol(Y) - 1] # this is the tstop

        time <- sort(unique(dtime)) # unique tstops

        nevent <- as.vector(rowsum(1 * death, dtime))

        nrisk <-
          rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
        delta <- min(diff(time)) / 2
        etime <-
          c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times

        indx <-
          approx(
            etime,
            1:length(etime),
            time,
            method = "constant",
            rule = 2,
            f = 1
          )$y

        esum <-
          rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
        nrisk <- nrisk - c(esum, 0)[indx]

        haz <- nevent / nrisk
        cumhaz <- cumsum(haz)

        out <- data.frame(time = time,
                          haz = haz,
                          cumhaz = cumhaz)
        return (out)
      }
      df <- getchz(Y = Y, explp = explp)
      t <- df$time
      haz <- df$haz
      H_0 <- df$cumhaz
      f <- stepfun(t[-1], H_0)
      H0 <- f(Y[, 2])
      #Survival Function
      SP <- exp(-z_hat * as.vector(exp(lp)) * H0)
      censored <- which(Y[, 3] == 0)
      n.censored <- length(censored)
      #Z-residual
      Zresid <- matrix(0, length(SP), n.rep)
      col_name <- rep(0, n.rep)
      for (i in 1:n.rep) {
        RSP <- SP
        RSP[censored] <- RSP[censored] * runif(n.censored)
        Zresid[, i] <- -qnorm(RSP)
        col_name[i] <- paste("Z-residual ", i, sep = "")
      }
      colnames(Zresid) <- col_name
      #####
      censored.status <- Y[, 3]

      Zresid.value <- as.matrix(Zresid)


      attributes(Zresid.value) <- c(
        attributes(Zresid.value),
        list(
          Survival.Prob = SP,
          linear.pred = lp,
          covariates = fix_var,
          censored.status = censored.status,
          object.model.frame = mf
        )
      )
    }

    ##########serve for cv function.##################
    if (!is.null(traindata) & !is.null(newdata)) {
      form <- (fit_coxph$formula)[[3]]
      group_id_name <- gsub(".*[(]([^.]+)[,].*", "\\1", form)[3]

      if (!is.factor(traindata[[group_id_name]]))
        stop("The group ID must be factor!")
      if (!is.factor(newdata[[group_id_name]]))
        stop("The group ID must be factor!")

      gpnumber <- length(levels(traindata[[group_id_name]]))
      mf <- model.frame(fit_coxph$formula, traindata)
      mf_nc <- ncol (mf)
      mm <- model.matrix(fit_coxph$formula, traindata)
      mm_nc <- ncol (mm)

      if (gpnumber > 5) {
        fix_var <- mm[, -c(1, mm_nc), drop = FALSE]
        explp <-
          exp(fix_var %*% fit_coxph$coefficients + fit_coxph$frail[mf[, mf_nc]])
        z_hat <- exp(fit_coxph$frail[mf[, mf_nc]])
      }
      if (gpnumber <= 5) {
        fix_var <- mm[, -c(1, mf_nc:mm_nc), drop = FALSE]
        coef_number <- ncol(fix_var)
        frailty <- fit_coxph$coefficients[coef_number + 1:gpnumber]
        explp <- exp(fix_var %*% fit_coxph$coefficients[1:coef_number] +
                       frailty[mf[, mf_nc]])
        z_hat <- exp(frailty[mf[, mf_nc]])
      }
      Y <- mf[[1]]

      if (!inherits(Y, "Surv"))
        stop("left hand side not a survival object")
      if (ncol(Y) != 3) {
        # making it all in (tstart, tstop) format
        Y <- Surv(rep(0, nrow(Y)), Y[, 1], Y[, 2])
      }

      # this one gives the baseline cumulative hazard at all the time points;
      getchz <- function(Y, explp) {
        death <-
          (Y[, ncol(Y)] == 1) # this is a TRUE FALSE for the status column
        dtime <- Y[, ncol(Y) - 1] # this is the tstop

        time <- sort(unique(dtime)) # unique tstops

        nevent <- as.vector(rowsum(1 * death, dtime))

        nrisk <-
          rev(cumsum(rev(rowsum(explp, dtime)))) # This gives the sum
        delta <- min(diff(time)) / 2
        etime <-
          c(sort(unique(Y[, 1])), max(Y[, 1]) + delta)  #unique entry times

        indx <-
          approx(
            etime,
            1:length(etime),
            time,
            method = "constant",
            rule = 2,
            f = 1
          )$y

        esum <-
          rev(cumsum(rev(rowsum(explp, Y[, 1]))))  #not yet entered
        nrisk <- nrisk - c(esum, 0)[indx]

        haz <- nevent / nrisk
        cumhaz <- cumsum(haz)

        out <- data.frame(time = time,
                          haz = haz,
                          cumhaz = cumhaz)
        return (out)
      }
      df <- getchz(Y = Y, explp = explp)
      t <- df$time
      haz <- df$haz
      H_0 <- df$cumhaz
      f <- stepfun(t[-1], H_0)

      full_data <- rbind(newdata, traindata)
      mf_new <-
        model.frame(fit_coxph$formula, full_data)[c(1:nrow(newdata)), ]
      mf_nc_new <- ncol (mf_new)
      mm_new <-
        model.matrix(fit_coxph$formula, full_data)[c(1:nrow(newdata)), , drop =
                                                     FALSE]
      mm_nc_new <- ncol (mm_new)

      if (gpnumber > 5) {
        fix_var_new <- mm_new[, -c(1, mm_nc_new), drop = FALSE]
        z_hat_new <- exp(fit_coxph$frail[mf_new[, mf_nc_new]])
        lp_new <-
          fix_var_new %*% fit_coxph$coefficients + fit_coxph$frail[mf_new[, mf_nc_new]]
        explp_new <- exp(fix_var_new %*% fit_coxph$coefficients)
      }
      if (gpnumber <= 5) {
        fix_var_new <- mm_new[, -c(1, mf_nc_new:mm_nc_new), drop = FALSE]
        coef_number_new <- ncol(fix_var_new)
        frailty_new <-
          fit_coxph$coefficients[coef_number_new + 1:gpnumber]
        lp_new <-
          fix_var_new %*% fit_coxph$coefficients[1:coef_number_new] + frailty_new[mf_new[, mf_nc_new]]
        explp_new <-
          exp(fix_var_new %*% fit_coxph$coefficients[1:coef_number_new])
        z_hat_new <- as.numeric(exp(frailty_new[mf_new[, mf_nc_new]]))
      }
      Y_new <- mf_new[[1]]

      if (!inherits(Y_new, "Surv"))
        stop("left hand side not a survival object")
      if (ncol(Y_new) != 3) {
        # making it all in (tstart, tstop) format
        Y_new <- Surv(rep(0, nrow(Y_new)), Y_new[, 1], Y_new[, 2])
      }
      H0_new <- f(Y_new[, 2])
      #Survival Function
      SP <- exp(-z_hat_new * as.vector(explp_new) * H0_new)
      censored <- which(Y_new[, 3] == 0)
      n.censored <- length(censored)

      #Z-residual
      Zresid <- matrix(0, length(SP), n.rep)
      col_name <- rep(0, n.rep)
      for (i in 1:n.rep) {
        RSP <- SP
        RSP[censored] <- RSP[censored] * runif(n.censored)
        Zresid[, i] <- -qnorm(RSP)
        col_name[i] <- paste("Z-residual ", i, sep = "")
      }
      colnames(Zresid) <- col_name
      #####S
      censored.status <- Y_new[, 3]

      Zresid.value <- as.matrix(Zresid)

      attributes(Zresid.value) <- c(
        attributes(Zresid.value),
        list(
          Survival.Prob = SP,
          linear.pred = lp_new,
          covariates = fix_var_new,
          censored.status = censored.status,
          object.model.frame = mf_new,
          type = "survival"

        )
      )

    }

    return(Zresid.value)

  }
