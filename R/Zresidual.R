#' A title
#'
#' Some descriptions
#'
#' @param fit.object The fit object are one of 'coxph', 'survreg' and 'brms'.
#'
#'
#' @export
#'
#' @return \itemize{
#'  \item{Zresid} {Z-residual}
#'
#' }
#'
Zresidual <- function(fit.object, nrep = 1,data = NULL,type=NULL,method = "iscv")
{
  get_object_name <- class(fit.object)

  if (inherits(fit.object, "coxph")) {
    frailty_terms <- attr(fit.object$terms, "specials")$frailty

    if (!is.null(frailty_terms)) {

      Zresid_fun <- Zresidual.coxph.frailty(fit_coxph = fit.object,traindata = data,
                                            newdata = data,n.rep = nrep)

    } else {

      Zresid_fun <-Zresidual.coxph(fit_coxph = fit.object,
                                   newdata = data,n.rep = nrep)

    }
  } else if (inherits(fit.object, "survreg")) {

    Zresid_fun <-Zresidual.survreg(fit_survreg = fit.object,
                                   newdata = data,n.rep = nrep)

  } else if (inherits(fit.object, "brmsfit")) {

    distr<- family(fit.object)$family

    if (distr =="hurdle_negbinomial") {

      Zresid_fun <- Zresidual.hurdle.negbinomial(fit = fit.object, type=type,
                                                 method = method,nrep = nrep)

    } else if (distr == "hurdle_poisson") {

      Zresid_fun <- Zresidual.hurdle.poisson(fit = fit.object,  type=type,
                                             method = method,nrep = nrep)

    } else if (distr == "negbinomial") {

      Zresid_fun <- Zresidual.negbinomial(fit = fit.object,  type = "NB",
                                          method = method,nrep = nrep)

    } else if (distr == "poisson") {

      Zresid_fun <- Zresidual.poisson(fit = fit.object, type = "pois",
                                      method = method,n.rep = nrep)

    } else if (distr == "bernoulli") {

      Zresid_fun <- Zresidual.bernoulli(fit=fit.object, method = method,nrep = nrep)

    } else {
      stop("The distribution '", distr, "' from the brmsfit object is not supported.")
    }

  }  else {
    stop("Objects of class '", paste(class(fit.object), collapse = "', '"), "' are not supported by this function.")
  }

  class(Zresid_fun) <- c("zresid", class(Zresid_fun))
  Zresid_fun
}


# else if (inherits(fit.object, "glmmTMB")) {
#
#   distr <- family(fit.object)$family
#
#   if (distr %in% c("gaussian", "poisson", "nbinom2")) {
#
#     Zresid_fun <- Zresidual.ZI(fit_ZI = fit.object,n.rep = nrep)
#
#   } else if (distr %in% c("truncated_poisson", "truncated_nbinom2")) {
#
#     Zresid_fun <- Zresidual.hurdle(fit_hd=fit.object,n.rep = nrep)
#
#   } else {
#     stop("The distribution '", distr, "' from the glmmTMB object is not supported.")
#   }
# }

