#' Cross-validated Z-residual diagnostics (generic)
#'
#' @description
#' Generic function for cross-validated Z-residual diagnostics.
#' Method dispatch is based on the class of \code{object}. Methods are
#' currently provided for survival models such as \code{coxph} and
#' \code{survreg}.
#'
#' @param object A fitted model object.
#' @param nfolds Integer. Number of folds for cross-validation.
#' @param foldlist Optional list specifying custom fold assignments. If
#'   \code{NULL}, folds are generated internally by the method.
#' @param data Optional data frame used to refit the model during
#'   cross-validation, when required by the method.
#' @param nrep Integer. Number of repeated cross-validations to perform.
#'   Default is 1.
#' @param ... Further arguments passed on to specific methods.
#'
#' @return
#' An object whose structure depends on the underlying method, typically
#' tagged with class \code{"cvzresid"} in addition to method-specific
#' classes.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'   fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   out <- CV.Zresidual(fit, nfolds = 5, data = lung)
#' }
#'
#' @export
CV.Zresidual <- function(object,
                         nfolds,
                         foldlist = NULL,
                         data     = NULL,
                         nrep     = 1, ...) {
  UseMethod("CV.Zresidual")
}
