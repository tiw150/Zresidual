#' Cross-Validated Z-Residual Diagnostics for Survival Models
#'
#' @description
#' Computes cross-validated Z-residuals for fitted survival models, including
#' Cox proportional hazards models (`coxph`) with or without frailty terms,
#' and parametric survival regression models (`survreg`).
#' The function automatically detects the model type from the fitted model
#' object and applies the appropriate Z-residual cross-validation method.
#'
#' @param fit.object A fitted survival model object of class `coxph` or `survreg`.
#' @param nfolds Integer. Number of folds for cross-validation.
#' @param foldlist Optional list specifying custom fold assignments. If `NULL`,
#'   folds are generated internally.
#' @param data Optional dataset used to refit the model during cross-validation.
#'   Required when `foldlist` is provided or when the original model call
#'   does not contain the data explicitly.
#' @param nrep Integer. Number of repeated cross-validations to perform.
#'   Default is 1.
#'
#' @details
#' The function identifies whether the fitted model is:
#' - a Cox model (`coxph`) with frailty terms,
#' - a Cox model without frailty,
#' - or a parametric survival model (`survreg`),
#' and dispatches to the appropriate internal cross-validation function:
#' `CV.Zresidual.coxph.frailty()`, `CV.Zresidual.coxph()`, or
#' `CV.Zresidual.survreg()`.
#'
#' All required packages are loaded via `pacman::p_load()`.
#'
#' @return
#' An object of class `"cvzresid"` containing the cross-validated Z-residual
#' results and any model-specific diagnostic information.
#'
#' @seealso
#' `CV.Zresidual.coxph()`, `CV.Zresidual.coxph.frailty()`,
#' `CV.Zresidual.survreg()`
#'
#' @examples
#' \dontrun{
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' out <- CV.Zresidual(fit, nfolds = 5)
#' }
#'
#' @export

CV.Zresidual<- function(fit.object, nfolds, foldlist=NULL, data=NULL,nrep=1)
{
  # Required packages:
  if (!requireNamespace("pacman")) {install.packages("pacman")}
  pacman::p_load("survival", "stringr")
  form<-fit.object$call
  get_object_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[1]

  if( get_object_name=="coxph") {
    frailty_terms<- attr(fit.object$terms, "specials")$frailty

    if (!is.null(frailty_terms)){
      CV.Zresid.fun <-CV.Zresidual.coxph.frailty(fit.coxph=fit.object, data=data,
                                                 nfolds=nfolds,foldlist=foldlist,
                                                 n.rep=nrep)

    } else CV.Zresid.fun <- CV.Zresidual.coxph(fit.coxph=fit.object, data=data,
                                               nfolds=nfolds,foldlist=foldlist,
                                               n.rep=nrep)
  }

  if (get_object_name=="survreg") {
    CV.Zresid.fun<-CV.Zresidual.survreg(fit.survreg=fit.object,data=data,
                                        nfolds=nfolds,foldlist=foldlist,
                                        n.rep=nrep)
    }
  class(CV.Zresid.fun) <- c("cvzresid", class(CV.Zresid.fun))
  CV.Zresid.fun
}



