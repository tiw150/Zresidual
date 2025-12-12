#' Z-residuals for hurdle Poisson models fitted with brms
#'
#' @description
#' `Zresidual.hurdle_poisson.brms()` is the S3 method for [Zresidual()]
#' when applied to hurdle Poisson models fitted with [brms::brm()]
#' and `family = hurdle_poisson()`. Objects are dispatched here when
#' the fitted object is a `"brmsfit"` with
#' `brms::family(object)$family == "hurdle_poisson"` and has been
#' internally tagged with the class `"hurdle_poisson.brms"` by
#' [Zresidual()].
#'
#' In normal use, users should call [Zresidual()] directly on the
#' `brmsfit` object (for example `Zresidual(fit)`), rather than calling
#' `Zresidual.hurdle_poisson.brms()` explicitly.
#'
#' @param object A `brmsfit` object with hurdle Poisson family
#'   (`brms::family(object)$family == "hurdle_poisson"`).
#' @param type Character string specifying which part of the model to
#'   compute Z-residuals for:
#'   \itemize{
#'     \item `"zero"`   — the hurdle/zero part;
#'     \item `"count"`  — the truncated Poisson count part;
#'     \item `"hurdle"` — the full hurdle-Poisson model.
#'   }
#'   The default is `"hurdle"`.
#' @param method Character string specifying the residual calculation method:
#'   `"iscv"` for importance-sampled cross-validated randomized predictive
#'   p-values, `"rpost"` for randomized posterior predictive p-values,
#'   or `"mpost"` for middle-value posterior predictive p-values.
#'   Default is `"iscv"`.
#' @param nrep Integer; number of replicated Z-residual sets to generate.
#'   Default is `1`.
#' @param data Optional data frame used for prediction or residual computation.
#' If \code{NULL} (default), the data stored inside the \code{brmsfit} object are used.
#' @param ... Further arguments passed to the underlying implementation
#'   function [Zresidual_hurdle_poisson_brms()].
#' @return
#' A numeric matrix of Z-residuals (one column per replication) as returned
#' by [Zresidual_hurdle_poisson_brms()], with the class `"zresid"`
#' added to its class vector.
#'
#' @examples
#' \dontrun{
#'   library(brms)
#'   fit_hp <- brm(y ~ x1 + x2, data = df,
#'                 family = hurdle_poisson())
#'
#'   ## Counts part only
#'   z_count <- Zresidual(fit_hp, type = "count", method = "iscv")
#'
#'   ## Full hurdle model with 2 replicates
#'   z_hurdle <- Zresidual(fit_hp, type = "hurdle", nrep = 2)
#' }
#'
#' @method Zresidual hurdle_poisson.brms
#' @export
Zresidual.hurdle_poisson.brms <- function(object,nrep   = 1, data,
                                          type   = c("hurdle", "count", "zero"),
                                          method = "iscv",
                                           ...) {

  type <- match.arg(type)

  out <- Zresidual_hurdle_poisson_brms(
    fit    = object,
    type   = type,
    method = method,
    n.rep  = nrep, ...
  )

  class(out) <- c("zresid", class(out))
  out
}


#' Compute Z-residuals for hurdle or count Poisson brms models
#'
#' @description
#' Computes Z-residuals for fitted Bayesian hurdle or count models
#' with a Poisson distribution using a \pkg{brms} model with
#' `family = hurdle_poisson()`. Z-residuals can be calculated for
#' zeros, counts, or the overall hurdle distribution, and can be used
#' for model diagnostics.
#'
#' This is an internal workhorse for [Zresidual.hurdle_poisson.brms()]
#' and is not intended to be called directly by end users.
#'
#' @param fit A fitted \pkg{brms} model object for a hurdle or count
#'   Poisson outcome.
#' @param type Character string specifying which part of the model to
#'   calculate Z-residuals for:
#'   `"zero"` for the hurdle/zero portion,
#'   `"count"` for the truncated Poisson counts,
#'   `"hurdle"` for the full hurdle-Poisson model.
#' @param method Character string specifying the residual calculation method:
#'   `"iscv"` for importance-sampled cross-validated randomized predictive
#'   p-values, `"rpost"` for randomized posterior predictive p-values,
#'   or `"mpost"` for middle-value posterior predictive p-values.
#'   Default is `"iscv"`.
#' @param n.rep Integer; the number of replicated Z-residual sets to
#'   generate. Default is `1`.
#'
#' @details
#' A typical implementation:
#' \enumerate{
#'   \item Extracts the observed response vector from the model data.
#'   \item Computes the log-PMF and log-CDF for the specified part of the
#'         model using the corresponding \code{log_pred_dist_*} function,
#'         such as \code{\link{log_pred_dist_HP}} or
#'         \code{\link{log_pred_dist_TP}}.
#'   \item Generates posterior predictive p-values according to the
#'         specified \code{method}.
#'   \item Converts the p-values to Z-residuals via the negative quantile
#'         of the standard normal distribution.
#' }
#'
#' The output is a matrix of Z-residuals with one column per replication.
#'
#' @return A numeric matrix of Z-residuals with attributes such as:
#' \itemize{
#'   \item \code{type}: The requested model component.
#'   \item \code{zero_id}: Indices of zero outcomes.
#'   \item \code{log_pmf}: Log-probability mass function values.
#'   \item \code{log_cdf}: Log-cumulative distribution function values.
#'   \item \code{covariates}: Model covariates.
#'   \item \code{linear.pred}: Linear predictor values from the fitted model.
#' }
#' The S3 wrapper [Zresidual.hurdle_poisson.brms()] will additionally
#' attach the class \code{"zresid"} to the returned object.
#'
#' @examples
#' \dontrun{
#'   # Compute Z-residuals for counts
#'   zres_counts <- Zresidual_hurdle_poisson_brms(
#'     fit    = fit_hp,
#'     type   = "count",
#'     method = "iscv"
#'   )
#' }
#'
#' @seealso
#' \code{\link{log_pred_dist_HP}}, \code{\link{log_pred_dist_TP}},
#' \code{\link{post_logrpp}}, \code{\link{iscv_logrpp}},
#' and the S3 wrapper [Zresidual.hurdle_poisson.brms()].
#'
#' @keywords internal
Zresidual_hurdle_poisson_brms <- function(fit,  type , method = "iscv", n.rep = 1, ...){

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  sim.y <- as.matrix(model.data[, response])
  n <- length(sim.y)
  #if(type == "count") id <- which(sim.y > 0) else id <- 1:n
  zero_id <- which(sim.y == 0)

  # Argument should be one of the element in type_list
  type_list <- c("zero", "TP", "HP")
  names(type_list) <- c("zero", "count", "hurdle")

  ldist <- get(paste0("log_pred_dist_", type_list[type]))(fit)
  lpmf <- ldist$lpmf_hat
  lcdf <- ldist$lcdf_hat

  # Argument should be one of the element in rpp_list
  rpp_list <- c(iscv = "iscv_logrpp", post = "post_logrpp", "post_logmpp")
  names(rpp_list) <- c("iscv", "rpost", "mpost")

  #if(count_only) z_res<- matrix(NA, ncol = n.rep, nrow = dim(lpmf)[2])
  z_res<- matrix(NA, ncol = n.rep, nrow = n)

  for (i in 1:n.rep) {
    rpp <- get(rpp_list[[method]])(lcdf, lpmf)
    z_res[, i] <- -qnorm(rpp, log.p = T)
  }

  #if(type == "count") z_res <- z_res[-zero_id,]

  colnames(z_res) <- paste0("Z-residual ", 1:n.rep)

  attributes(z_res) <- c(attributes(z_res),list(
    type = type,
    #count_only = count_only,
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates = fit$data[, setdiff(names(fit$data), response), drop = FALSE],
    linear.pred = fitted(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
