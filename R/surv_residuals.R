#' Compute residual diagnostics for supported survival models
#'
#' @description
#' Unified wrapper for residual diagnostics from supported survival models.
#' Depending on the fitted model class, the function dispatches internally to
#' the appropriate residual computation for Cox models, shared-frailty Cox
#' models, or parametric survival regression models.
#'
#' @param fit.object A fitted survival model object.
#' @param data A data frame containing the variables needed to evaluate
#'   residuals.
#' @param residual.type Character string specifying the residual type.
#'
#' @return
#' A residual object whose structure depends on \code{residual.type}. The result
#' is usually numeric with attributes used by downstream plotting functions.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.3 * x))
#'   t_cens  <- rexp(n, rate = 0.5)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   r <- surv_residuals(fit, data = dat, residual.type = "martingale")
#'   head(as.vector(r))
#' }
#'
#' @export
surv_residuals <- function(fit.object, data,
                           residual.type = c("censored Z-residual", "Cox-Snell",
                                             "martingale", "deviance")) {
  residual.type <- match.arg(residual.type)
  
  if (inherits(fit.object, "coxph")) {
    frailty_terms <- attr(fit.object$terms, "specials")$frailty
    
    if (!is.null(frailty_terms) && length(frailty_terms) > 0L) {
      resid_fun <- residual.coxph.frailty(
        fit_coxph = fit.object,
        traindata = data,
        newdata = data,
        residual.type = residual.type
      )
    } else {
      resid_fun <- residual.coxph(
        fit_coxph = fit.object,
        newdata = data,
        residual.type = residual.type
      )
    }
    
  } else if (inherits(fit.object, "survreg")) {
    resid_fun <- residual.survreg(
      survreg_fit = fit.object,
      newdata = data,
      residual.type = residual.type
    )
    
  } else {
    stop(
      "surv_residuals() currently supports objects of class 'coxph' and 'survreg'.",
      call. = FALSE
    )
  }
  
  class(resid_fun) <- c("survresid", class(resid_fun))
  resid_fun
}

