## R/Zresidual-poisson-brms.R

#' Z-residuals for Poisson models fitted with brms
#'
#' @description
#' `Zresidual.poisson.brms()` is the S3 method for [Zresidual()]
#' when applied to Poisson models fitted with [brms::brm()]
#' and `family = poisson()`. Objects are dispatched here when
#' the fitted object is a `"brmsfit"` with
#' `brms::family(object)$family == "poisson"` and has been
#' internally tagged with the class `"poisson.brms"` by
#' [Zresidual()].
#'
#' In normal use, users should call [Zresidual()] directly on the
#' `brmsfit` object (for example `Zresidual(fit)`), rather than calling
#' `Zresidual.poisson.brms()` explicitly.
#'
#' @param object A `brmsfit` object with Poisson family
#'   (`brms::family(object)$family == "poisson"`).
#' @param nrep Integer; the number of replicated Z-residual sets to
#'   generate. Default is `1`.
#' @param data Optional data frame used for prediction or residual computation.
#' If \code{NULL} (default), the data stored inside the \code{brmsfit} object are used.
#' @param type Optional character string controlling the residual type,
#' interpreted by the underlying implementation (if used).
#' @param method Character string specifying the residual calculation method:
#'   `"iscv"` for importance-sampled cross-validated randomized predictive
#'   p-values, `"rpost"` for posterior predictive p-values, or
#'   `"mpost"` for marginal posterior predictive p-values. Default is
#'   `"iscv"`.
#' @param ... Further arguments passed from [Zresidual()]. They are ignored
#'   by this method but are accepted for consistency with the generic.
#'
#' @return
#' A numeric matrix of Z-residuals (one column per replication) as returned
#' by [Zresidual_poisson_brms()], with the class `"zresid"` added to
#' its class vector.
#'
#' @examples
#' \dontrun{
#'   library(brms)
#'   fit_pois <- brm(y ~ x1 + x2, data = df,
#'                   family = poisson())
#'
#'   ## ISCV-based Z-residuals
#'   z_pois <- Zresidual(fit_pois, method = "iscv")
#'
#'   ## Posterior predictive Z-residuals with 2 replicates
#'   z_pois_post <- Zresidual(fit_pois, method = "rpost", nrep = 2)
#' }
#'
#' @method Zresidual poisson.brms
#' @export
Zresidual.poisson.brms <- function(object,
                                   nrep   = 1,data, type,
                                   method = "iscv", ...) {

  out <- Zresidual_poisson_brms(
    fit    = object,
    method = method,
    n.rep  = nrep, ...
  )

  class(out) <- c("zresid", class(out))
  out
}


#' Compute Z-residuals for Poisson brms models
#'
#' @description
#' Computes Z-residuals for fitted Bayesian Poisson models.
#' Z-residuals are useful for model diagnostics, including checking fit and
#' overdispersion, and can be calculated using posterior or cross-validated
#' predictive p-values.
#'
#' This is an internal workhorse for [Zresidual.poisson.brms()] and
#' is not intended to be called directly by end users.
#'
#' @param fit A fitted \pkg{brms} model object for a Poisson outcome.
#' @param method Character string specifying the residual calculation method:
#'   \code{"iscv"} for importance-sampled cross-validated randomized predictive
#'   p-values, \code{"rpost"} for posterior predictive p-values, or
#'   \code{"mpost"} for marginal posterior predictive p-values. Default is
#'   \code{"iscv"}.
#' @param n.rep Integer; the number of replicated Z-residual sets to generate.
#'   Default is 1.
#' @param ... Further arguments passed to lower-level helper functions.
#'
#' @details
#' The function typically performs the following steps:
#' \enumerate{
#'   \item Extracts the observed response vector from the model data.
#'   \item Computes the log-PMF and log-CDF for the Poisson model
#'         using \code{\link{log_pred_dist_pois}}.
#'   \item Generates randomized or posterior predictive p-values according to the
#'         specified \code{method}.
#'   \item Converts the p-values to Z-residuals via the negative quantile of the
#'         standard normal distribution.
#' }
#'
#' The output is a matrix of Z-residuals with one column per replication.
#'
#' @return A numeric matrix of Z-residuals with attributes such as:
#' \itemize{
#'   \item \code{zero_id}: Indices of zero outcomes.
#'   \item \code{log_pmf}: Log-probability mass function values.
#'   \item \code{log_cdf}: Log-cumulative distribution function values.
#'   \item \code{covariates}: Model covariates.
#'   \item \code{linear.pred}: Linear predictor values from the fitted model.
#' }
#' The S3 wrapper [Zresidual.poisson.brms()] will additionally attach
#' the class \code{"zresid"} to the returned object.
#'
#' @examples
#' \dontrun{
#'   # Compute Z-residuals for a Poisson model
#'   zres_pois <- Zresidual_poisson_brms(
#'     fit    = fit_pois,
#'     method = "iscv"
#'   )
#'
#'   # Compute Z-residuals with 2 replicates using posterior predictive p-values
#'   zres_pois_post <- Zresidual_poisson_brms(
#'     fit    = fit_pois,
#'     method = "rpost",
#'     n.rep  = 2
#'   )
#' }
#'
#' @seealso
#' \code{\link{log_pred_dist_pois}}, \code{\link{post_logrpp}},
#' \code{\link{iscv_logrpp}}, and the S3 wrapper
#' [Zresidual.poisson.brms()].
#'
#' @keywords internal
Zresidual_poisson_brms <- function(fit, method = "iscv", n.rep = 1, ...){

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  zero_id <- which(sim.y == 0)

  # Argument should be one of the element in type_list
  type_list <- c("pois")
  names(type_list) <- c("pois")

  ldist <- log_pred_dist_pois(fit)
  lpmf <- ldist$lpmf_hat
  lcdf <- ldist$lcdf_hat


  rpp_list <- c(iscv = "iscv_logrpp", post = "post_logrpp", "post_logmpp")
  names(rpp_list) <- c("iscv", "rpost", "mpost")

  z_res<- matrix(NA, ncol = n.rep, nrow = n)

  for (i in 1:n.rep) {
    rpp <- get(rpp_list[[method]])(lcdf, lpmf)
    z_res[, i] <- -qnorm(rpp, log.p = T)
  }

  colnames(z_res) <- paste0("Z-residual ", 1:n.rep)

  attributes(z_res) <- c(attributes(z_res),list(
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates =fit$data[, setdiff(names(fit$data), response), drop = FALSE],
    linear.pred = predict(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
