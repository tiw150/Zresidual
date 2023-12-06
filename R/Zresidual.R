#' A title
#'
#' Some descriptions
#'
#' @param fit.object The fit object are one of 'coxph', 'survreg' and 'glmmTMB'.
#'
#' @param data Data that used for fitting the survival model.
#'
#' @export
#'
#' @return \itemize{
#'  \item{Zresid}{Z-residual}
#'  \item{SP}{Survival Probabilities}

#' }
#'
Zresidual <-
  function(fit.object,
           nrep = 1,
           data = NULL)
  {
    form <- fit.object$call

    # if (is.na(form)) stop("a fit.object is required")

    get_object_name <- gsub(".*[(]([^.]+)[,].*", "\\1", form)[1]

    # the coxph function in the survival package
    if (get_object_name == "coxph") {
      frailty_terms <- attr(fit.object$terms, "specials")$frailty

      if (!is.null(frailty_terms)) {
        Zresid_fun <- Zresidual.coxph.frailty(
          fit_coxph = fit.object,
          traindata = data,
          newdata = data,
          n.rep = nrep
        )
      } else
        Zresid_fun <-
          Zresidual.coxph(fit_coxph = fit.object,
                          newdata = data,
                          n.rep = nrep)
    }

    # the survreg function in the survival package
    if (get_object_name == "survreg") {
      Zresid_fun <-
        Zresidual.survreg(survreg_fit = fit.object,
                          newdata = data,
                          n.rep = nrep)
    }

    # the glmmTMB function in the glmmTMB package
    if (get_object_name == "glmmTMB") {
      distr <- family(fit.object)$family

      if (distr[1] %in% c("gaussian", "poisson", "nbinom2")) {
        Zresid_fun <- Zresidual.ZI(fit_ZI = fit.object,n.rep = nrep)

      } else if (distr %in% c("truncated_poisson", "truncated_nbinom2")) {

        Zresid_fun <- Zresidual.hurdle(fit_hd=fit.object,n.rep = nrep)
      } else
        stop ("The distribution is not supported")
    }

    class(Zresid_fun) <- c("zresid", class(Zresid_fun))
    Zresid_fun
  }
