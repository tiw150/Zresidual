#' Z-residuals for Cox proportional hazards models (survival package)
#'
#' @description
#' `Zresidual.coxph.survival()` computes randomized Z-residuals for Cox
#' proportional hazards models fitted with \code{\link[survival]{coxph}},
#' supporting both standard and shared frailty models. This S3 method is
#' designed to be called via the generic function \code{\link{Zresidual}}.
#'
#' The function automatically detects the presence of a frailty term (e.g.,
#' \code{frailty(group)}) and dispatches the calculation to one of two
#' internal implementations. These residuals are intended for **in-sample**
#' diagnostics and model assessment.
#'
#' @details
#' This method dispatches work based on the model formula:
#'
#' \itemize{
#'   \item \strong{Standard Cox models (No Frailty)}:
#'
#'   The function calls \code{\link{Zresidual_coxph_survival}} to compute
#'   Z-residuals using the fixed effects (\eqn{\mathbf{x}\hat{\mathbf{\beta}}}{x * beta-hat}) and the
#'   estimated baseline cumulative hazard function \eqn{\hat{H}_0(t)}{H\_0\_hat(t)}.
#'
#'   \item \strong{Shared Frailty Cox models}:
#'
#'   The function calls \code{\link{Zresidual_coxph_frailty_survival}}. This
#'   implementation computes residuals accounting for the cluster-level frailty
#'   effect (\eqn{\hat{z}_{\text{group}}}{z\_group\_hat}). It requires the data used for fitting
#'   (\code{traindata}) to reconstruct the baseline hazard and estimate the
#'   frailty term.
#' }
#'
#' **Randomization for Censored Observations:**
#' Since the true survival probability for a censored observation $i$ is known
#' only to be greater than \eqn{S_i(t_i)}{S\_i(t\_i)}, the Z-residual uses a randomized survival
#' probability: \eqn{S_{i, \text{rand}}(t_i) = S_i(t_i) \cdot U}{S\_i\_rand(t\_i) = S\_i(t\_i) \* U}, where \eqn{U \sim \text{Unif}(0, 1)}{U \~ Unif(0, 1)}.
#' This randomization is repeated \code{nrep} times.
#'
#' @importFrom survival Surv
#'
#' @param object A fitted \code{\link[survival]{coxph}} model. Supports
#'   both standard Cox models and shared frailty models.
#' @param data Optional \code{data.frame} containing the survival response
#'   and covariates. When \code{NULL} (default), the residuals are computed
#'   on the data used to fit the \code{object}. This parameter is often aliased
#'   as \code{newdata} in the internal worker functions.
#' @param nrep Integer; number of independent randomized Z-residual
#'   replicates to generate. Defaults to \code{1}.
#' @param type Optional character string controlling the residual type.
#'   Set internally to \code{"survival"} for Cox models.
#' @param method Character string specifying the residual calculation method.
#'   Currently unused.
#' @param ... Further arguments passed to the underlying implementation
#'   functions.
#' @return
#' A numeric matrix of class \code{"zresid"} with dimension \eqn{N \times nrep}{N x nrep}.
#' Each column is an independent set of Z-residuals. The following diagnostic
#' attributes are attached:
#' \itemize{
#'   \item \code{Survival.Prob}: Vector of predicted survival probabilities \eqn{S_i(t_i)}{S\_i(t\_i)}.
#'   \item \code{linear.pred}: Vector of linear predictors \eqn{\eta_i = \mathbf{x}_i \mathbf{\hat{\beta}}}{eta\_i = x\_i * beta-hat}.
#'   \item \code{covariates}: Data frame of covariates used in the model.
#'   \item \code{censored.status}: Event indicator (1 = event, 0 = censored).
#'   \item \code{object.model.frame}: The \code{model.frame} used for computation.
#'   \item \code{type}: Character string, always \code{"survival"}.
#' }
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   ## Standard Cox model (no frailty term)
#'   fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   # Note: The internal class 'coxph.survival' is usually added by a wrapper,
#'   # but Zresidual() handles dispatch automatically.
#'   z_cox <- Zresidual(fit_cox, nrep = 10, data = lung)
#'
#'   ## Shared frailty Cox model (in-sample residuals)
#'   # Note: 'inst' must be a grouping factor.
#'   lung$inst_f <- factor(lung$inst)
#'   fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(inst_f),
#'                      data = lung)
#'   z_in <- Zresidual(fit_frail, nrep = 5)
#' }
#'
#' @seealso
#' \code{\link{Zresidual}},
#' \code{\link{Zresidual_coxph_survival}},
#' \code{\link{Zresidual_coxph_frailty_survival}}
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

#' @title Z-residuals for Standard Cox Models (Internal Worker)
#' @name Zresidual_coxph_survival
#' @keywords internal
#' @description Internal function to compute randomized Z-residuals for a **standard**
#' Cox proportional hazards model (without frailty).
#'
#' @param fit_coxph A fitted \code{survival::coxph} model object.
#' @param newdata Optional data frame on which to compute the residuals. If \code{NULL},
#'   the original model frame is used.
#' @param n.rep Integer. Number of randomized residual replicates to generate. Default is 1.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The Z-residual for an observation $i$ is calculated as
#' \deqn{Z_i = -\Phi^{-1}(\hat{S}_i(t_i, \text{rand}))}{Z\_i = -qnorm(S\_hat\_i(t\_i, rand))}
#' where \eqn{\Phi^{-1}}{Phi-inverse} is the inverse standard normal CDF, and \eqn{\hat{S}_i(t_i, \text{rand})}{S\_hat\_i(t\_i, rand)}
#' is the predicted survival probability at time \eqn{t_i}{t\_i}.
#'
#' For uncensored observations (\eqn{t_i}{t\_i} is an event time), \eqn{\hat{S}_i(t_i, \text{rand}) = \hat{S}_i(t_i)}{S\_hat\_i(t\_i, rand) = S\_hat\_i(t\_i)}.
#' For censored observations,
#' \eqn{\hat{S}_i(t_i, \text{rand}) = \hat{S}_i(t_i) \cdot U}{S\_hat\_i(t\_i, rand) = S\_hat\_i(t\_i) \* U}, where \eqn{U \sim \text{Unif}(0, 1)}{U \~ Unif(0, 1)}.
#'
#' The predicted survival is calculated as
#' \deqn{\hat{S}_i(t_i) = \exp(-\exp(\mathbf{x}_i \mathbf{\hat{\beta}}) \hat{H}_0(t_i))}{S\_hat\_i(t\_i) = exp(-exp(x\_i * beta-hat) * H\_0\_hat(t\_i))}.
#'
#' @return A matrix containing the Z-residuals (\eqn{N \times nrep}{N x nrep}) with diagnostic attributes:
#' \code{Survival.Prob}, \code{linear.pred}, \code{covariates}, \code{censored.status},
#' \code{object.model.frame}, and \code{type = "survival"}.
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

#' @title Z-residuals for Standard Cox Models with Frailty terms (Internal Worker)
#' @name Zresidual_coxph_frailty_survival
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
