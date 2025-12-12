#' Z-residuals for Bernoulli models fitted with brms
#'
#' @description
#' `Zresidual.bernoulli.brms()` is the S3 method for [Zresidual()] when
#' applied to Bernoulli (binary) regression models fitted with
#' [brms::brm()] and `family = bernoulli()`. Objects are dispatched here
#' when the fitted object is a `"brmsfit"` with family `"bernoulli"` and
#' has been internally tagged with the class `"bernoulli.brms"` by
#' [Zresidual()].
#'
#' In most cases users should call [Zresidual()] directly on the
#' `brmsfit` object, e.g. `Zresidual(fit)`, rather than calling
#' `Zresidual.bernoulli.brms()` explicitly.
#'
#' @param object A `brmsfit` object with `brms::family(object)$family == "bernoulli"`.
#' @param nrep Integer; number of replicated Z-residual sets to generate.
#'   Defaults to `1`.
#' @param data Optional data frame used for prediction. If `NULL`, the
#'   data stored inside the `brmsfit` object are used.
#' @param type Optional character string controlling the residual type,
#'   interpreted by the underlying implementation (if used).
#' @param method Character string specifying the residual calculation method:
#'   `"iscv"` for importance-sampled cross-validated randomized predictive
#'   p-values, `"rpost"` for randomized posterior predictive p-values, or
#'   `"mpost"` for middle-value posterior predictive p-values. Default is
#'   `"iscv"`.
#' @param ... Further arguments passed to the underlying implementation
#'
#' @return
#' A numeric matrix of Z-residuals with one column per replication, as
#' returned by [Zresidual_bernoulli_brms()], but with the class `"zresid"`
#' added to its class vector.
#'
#' @examples
#' \dontrun{
#'   library(brms)
#'   fit_bern <- brm(y ~ x1 + x2, data = df, family = bernoulli())
#'
#'   ## ISCV-based Z-residuals
#'   z1 <- Zresidual(fit_bern, method = "iscv", nrep = 2)
#'
#'   ## Posterior predictive Z-residuals
#'   z2 <- Zresidual(fit_bern, method = "rpost")
#' }
#'
#' @method Zresidual bernoulli.brms
#' @export
Zresidual.bernoulli.brms <- function(object,
                                     nrep   = 1,
                                     data   = NULL,
                                     type   = NULL,
                                     method = "iscv", ...) {

  out <- Zresidual_bernoulli_brms(
    fit    = object,
    method = method,
    n.rep  = nrep,
    data   = data,
    type   = type, ...
  )

  class(out) <- c("zresid", class(out))
  out
}


#' Compute Z-residuals for a Bernoulli/Logistic brms model
#'
#' @description
#' Computes Z-residuals for a fitted Bayesian Bernoulli/Logistic (binary)
#' model fitted with [brms::brm()] and `family = bernoulli()`. Z-residuals
#' are calculated using posterior predictive methods and can be used for
#' model diagnostics.
#'
#' This is an internal workhorse for [Zresidual.bernoulli.brms()] and is
#' not intended to be called directly by end users.
#'
#' @param fit A fitted `brmsfit` model object with Bernoulli family.
#' @param method Character string specifying the residual calculation method:
#'   `"iscv"` for importance-sampled cross-validated randomized predictive
#'   p-values, `"rpost"` for randomized posterior predictive p-values, or
#'   `"mpost"` for middle-value posterior predictive p-values. Default is
#'   `"iscv"`.
#' @param n.rep Integer; the number of replicated Z-residual sets to
#'   generate. Default is `1`.
#' @param data Optional data frame used to override the data stored inside
#'   `fit` for prediction and diagnostic calculation. If `NULL`, the data
#'   embedded in `fit` are used.
#' @param type Optional character string controlling the residual type; the
#'   meaning is determined by the underlying implementation (if used).
#' @param ... Further arguments passed to lower-level helpers.
#'
#' @details
#' The function typically performs the following steps:
#' \enumerate{
#'   \item Extracts the observed response vector from the model data.
#'   \item Computes the log-PMF and log-CDF for the Bernoulli model using
#'         \code{\link{log_pred_dist_bern}}.
#'   \item Generates posterior predictive p-values according to the
#'         specified \code{method}.
#'   \item Converts the p-values to Z-residuals via the negative quantile
#'         of the standard normal distribution.
#' }
#'
#' The output is a matrix of Z-residuals with one column per replication.
#'
#' @return A numeric matrix of Z-residuals with attributes:
#' \itemize{
#'   \item \code{type}: Type of outcome (Bernoulli).
#'   \item \code{zero_id}: Indices of zero outcomes.
#'   \item \code{log_pmf}: Log-probability mass function values.
#'   \item \code{log_cdf}: Log-cumulative distribution function values.
#'   \item \code{covariates}: Model covariates.
#'   \item \code{linear.pred}: Linear predictor values from the fitted model.
#' }
#'
#' @seealso
#' \code{\link{log_pred_dist_bern}}, \code{\link{post_logrpp}},
#' \code{\link{iscv_logrpp}},
#' and the S3 wrapper [Zresidual.bernoulli.brms()].
#'
#' @keywords internal
Zresidual_bernoulli_brms <- function(fit, method = "iscv", n.rep = 1,data = NULL,
                                     type = NULL, ...){

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
  # type <- "count"
  # type_list <- c("bern")
  # names(type_list) <- c("count")

  ldist <- log_pred_dist_bern(fit)
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
  #  type = type,
    #count_only = count_only,
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates = subset(fit$data, select = -get(response)),
    linear.pred = fitted(fit)[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
