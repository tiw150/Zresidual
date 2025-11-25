#' Compute Z-Residuals for Survival and Bayesian Regression Models
#'
#' @description
#' Computes Z-residuals for a variety of model types, including Cox proportional
#' hazards models (`coxph`) with or without frailty terms, parametric survival
#' regression models (`survreg`), and Bayesian regression models fitted using
#' **brms** (`brmsfit`) for several supported distributions (Poisson,
#' negative binomial, hurdle Poisson, hurdle negative binomial, and Bernoulli).
#'
#' The function automatically detects the model class and dispatches the call
#' to the appropriate Z-residual calculation method.
#'
#' @param fit.object A fitted model object. Supported classes include:
#'   * `"coxph"` - Cox proportional hazards model (from **survival**)
#'   * `"survreg"` - Parametric survival regression
#'   * `"brmsfit"` - Bayesian regression model from **brms**
#'
#' @param nrep Integer. Number of repeated Z-residual samples to compute.
#'   Default is 1.
#'
#' @param data Optional dataset used for generating new model predictions
#'   and residuals. Required for some models, especially when frailty terms
#'   are present in `coxph`.
#'
#' @param type Optional character string. Residual type for **brms** hurdle
#'   and count models. The meaning of this argument depends on the chosen
#'   model family.
#'
#' @param method Character string indicating which predictive p-value
#'   method to use for Z-residual computation in `brmsfit` models.
#'   Options include `"iscv"` (importance-sampling cross-validated),
#'   `"loocv"`, and `"posterior"` (default `"iscv"`).
#'
#' @details
#' The function performs a class-specific dispatch:
#'
#' * **coxph**
#'   - Detects frailty terms using `attr(fit.object$terms, "specials")$frailty`
#'   - Calls `Zresidual.coxph()` or `Zresidual.coxph.frailty()`
#'
#' * **survreg**
#'   - Calls `Zresidual.survreg()`
#'
#' * **brmsfit**
#'   - Detects distribution via `family(fit.object)$family`
#'   - Supported families:
#'       `"hurdle_negbinomial"`, `"hurdle_poisson"`,
#'       `"negbinomial"`, `"poisson"`, `"bernoulli"`
#'   - Dispatches to the corresponding Z-residual function
#'
#' Unsupported model classes or unsupported **brms** families produce an error.
#'
#' @return
#' An object of class `"zresid"` containing computed Z-residuals and any
#' model-specific diagnostic output produced by the underlying Z-residual
#' function.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' fit1 <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' z1 <- Zresidual(fit1, nrep = 10, data = lung)
#'
#' library(brms)
#' fit2 <- brm(count ~ x, data = df, family = poisson)
#' z2 <- Zresidual(fit2, method = "posterior")
#' }
#'
#' @seealso
#' `CV.Zresidual()`,
#' `Zresidual.coxph()`,
#' `Zresidual.survreg()`,
#' model-specific Z-residual methods for **brmsfit** families.
#'
#' @export

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
                                                 method = method,n.rep = nrep)

    } else if (distr == "hurdle_poisson") {

      Zresid_fun <- Zresidual.hurdle.poisson(fit = fit.object,  type=type,
                                             method = method,n.rep = nrep)

    } else if (distr == "negbinomial") {

      Zresid_fun <- Zresidual.negbinomial(fit = fit.object,  type = "NB",
                                          method = method,n.rep = nrep)

    } else if (distr == "poisson") {

      Zresid_fun <- Zresidual.poisson(fit = fit.object, type = "pois",
                                      method = method,n.rep = nrep)

    } else if (distr == "bernoulli") {

      Zresid_fun <- Zresidual.bernoulli(fit=fit.object, method = method,n.rep = nrep)

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

