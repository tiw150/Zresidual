residuals<- function(fit.object, data,
                             residual.type=c("censored Z-residual", "Cox-Snell",
                                             "martingale", "deviance"))
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

      resid_fun <-residual.coxph.frailty(fit_coxph=fit.object,traindata=data,
                                         newdata=data, residual.type=residual.type)

    } else resid_fun <- residual.coxph(fit_coxph=fit.object,
                                       newdata=data,residual.type=residual.type)
  }

  # the survreg function in the survival package
  if (get_object_name=="survreg") {

    resid_fun <- residual.survreg(survreg_fit=fit.object,newdata=data,
                                  residual.type=residual.type)
  }
  resid_fun
}



