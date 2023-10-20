Zresidual<- function(fit.object, data, fit.object2=NULL)
{
  # Required packages:
  if (!requireNamespace("pacman")) {install.packages("pacman")}
  pacman::p_load("survival", "stringr", "glmmTMB","actuar")
  form<- fit.object$call
  get_object_name<-gsub(".*[(]([^.]+)[,].*", "\\1", form)[1]

  # the coxph function in the survival package
  if( get_object_name=="coxph") {

    frailty_terms<- attr(fit.object$terms, "specials")$frailty

    if (!is.null(frailty_terms)){

      Zresid_fun <-Zresidual.coxph.frailty(fit_coxph=fit.object,traindata=data,newdata=data)

      } else Zresid_fun <- Zresidual.coxph(fit_coxph=fit.object,newdata=data)
    }

  # the survreg function in the survival package
  if (get_object_name=="survreg") {

    Zresid_fun <- Zresidual.survreg(survreg_fit=fit.object,newdata=data)
  }

  # the glmmTMB function in the glmmTMB package
  if (get_object_name=="glmmTMB") {

     distr<-family(fit.object)$family

     if(distr[1] %in% c("gaussian", "poisson","nbinom2")){

       Zresid_fun <- Zresidual.ZI(fit_ZI = fit.object)

     }else if  (distr %in% c("truncated_poisson", "truncated_nbinom2")) {

      if (is.null(fit.object2)) stop("a fit.object2 is required for modelling the probability of zero in zero-modified model.")
       Zresid_fun <- Zresidual.hurdle(model_count=fit.object, model_zero=fit.object2, data=data)

     }else stop ("The distribution is not supported")
   }


  Zresid_fun
}



