#' Extract point-level quantities for predictive checks
#'
#' @description
#' Thin wrapper around \code{\link{Zcov}} that requests the standardized
#' point-level quantities required by predictive checks.
#'
#' @param fit A fitted model object.
#' @param data Evaluation data.
#' @param type Optional component selector.
#' @param ndraws Optional number of posterior draws.
#' @param draw_ids Optional posterior draw indices. When supplied, this is used
#'   to synchronize point quantities with \code{log_pointpred} and
#'   \code{postpred_simulate}.
#' @param ... Additional arguments passed to \code{Zcov}.
#'
#' @return A standardized list containing at least \code{response},
#'   \code{mean}, \code{variance}, \code{lp}, and \code{covariate} when
#'   supported by the model backend.
#'
#' @export
point_Zcov <- function(fit,
                       data,
                       type = NULL,
                       ndraws = NULL,
                       draw_ids = NULL,
                       ...) {
  Zcov(
    fit = fit,
    data = data,
    type = type,
    point_details = c("mean", "variance", "lp", "covariate"),
    ndraws = ndraws,
    draw_ids = draw_ids,
    ...
  )
}
