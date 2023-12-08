#' A title
#'
#' Some descriptions
#'
#' @param fit.object The fit object are one of 'coxph' and 'survreg'.
#'
#' @param data Data that used for fitting the survival model.
#'
#' @export
#'
#' @return \itemize{
#'  \item{Zresid}{CV Z-residual}
#'  \item{SP}{Survival Probabilities}

#' }
#'

CV.Zresidual<- function(fit.object, foldlist=NULL, data=NULL, nfolds=NULL,nrep=1)
{
  # Required packages:
  if (!requireNamespace("pacman")) {install.packages("pacman")}
  pacman::p_load("survival", "stringr")
  form<-fit.object$call
  get_object_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[1]

  if( get_object_name=="coxph") {
    frailty_terms<- attr(fit.object$terms, "specials")$frailty

    if (!is.null(frailty_terms)){
      cv.Zresid.fun <-CV.Zresidual.coxph.frailty(fit.coxph=fit.object, data=data,
                                                 nfolds=nfolds,foldlist=foldlist,
                                                 n.rep=nrep)

    } else cv.Zresid.fun <- CV.Zresidual.coxph(fit.coxph=fit.object, data=data,
                                               nfolds=nfolds,foldlist=foldlist,
                                               n.rep=nrep)
  }

  if (get_object_name=="survreg") {
    cv.Zresid.fun<-CV.Zresidual.survreg(fit.survreg=fit.object,data=data,
                                        nfolds=nfolds,foldlist=foldlist,
                                        n.rep=nrep)
    }

  cv.Zresid.fun
}



