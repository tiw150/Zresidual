cv.Zresidual<- function(fit.object, data, nfolds=NULL,foldlist=NULL)
{
  # Required packages:
  if (!requireNamespace("pacman")) {install.packages("pacman")}
  pacman::p_load("survival", "stringr")
  form<-fit.object$call
  get_object_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[1]

  if( get_object_name=="coxph") {
    frailty_terms<- attr(fit.object$terms, "specials")$frailty

    if (!is.null(frailty_terms)){
      cv.Zresid.fun <-cv.zresidual.coxph.frailty(fit=fit.object, data=data,
                                                 nfolds=nfolds,foldlist=foldlist)

    } else cv.Zresid.fun <- cv.zresidual.coxph(fit=fit.object, data=data,
                                               nfolds=nfolds,foldlist=foldlist)
  }

  if (get_object_name=="survreg") {
    cv.Zresid.fun<-cv.zresidual.survreg(fit=fit.object,data=data,
                                        nfolds=nfolds,foldlist=foldlist)
    }

  cv.Zresid.fun
}



