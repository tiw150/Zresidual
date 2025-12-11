#' Compute Z-Residuals for Hurdle or Count Negative Binomial Models
#'
#' Computes Z-residuals for fitted Bayesian hurdle or count models
#' with a negative binomial distribution. Z-residuals can be calculated for
#' zeros, counts, or the overall hurdle distribution, and can be used
#' for model diagnostics.
#'
#' @param fit A fitted \pkg{brms} model object for a hurdle or count negative binomial outcome.
#' @param type Character string specifying which part of the model to calculate Z-residuals for:
#'   \code{"zero"} for the hurdle/zero portion,
#'   \code{"count"} for the truncated negative binomial counts,
#'   \code{"hurdle"} for the full hurdle-negative binomial model.
#' @param method Character string specifying the residual calculation method:
#'   \code{"iscv"} for importance-sampled cross-validated randomized predictive p-values,
#'   or \code{"rpost"} for randomized posterior predictive p-values or
#'   \code{"mpost"} for middle-value posterior predictive p-values. Default is \code{"iscv"}.
#' @param n.rep Integer; the number of replicated Z-residual sets to generate. Default is 1.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the observed response vector from the model data.
#'   \item Computes the log-PMF and log-CDF for the specified part of the model
#'         using the corresponding \code{log_pred_dist_*} function.
#'   \item Generates posterior predictive p-values according to the
#'         specified \code{method}.
#'   \item Converts the p-values to Z-residuals via the negative quantile of the
#'         standard normal distribution.
#' }
#'
#' The output is a matrix of Z-residuals with one column per replication.
#'
#' @return A numeric matrix of Z-residuals with attributes:
#' \itemize{
#'   \item \code{type}: The specified model component
#'   \item \code{zero_id}: Indices of zero outcomes
#'   \item \code{log_pmf}: Log-probability mass function values
#'   \item \code{log_cdf}: Log-cumulative distribution function values
#'   \item \code{covariates}: Model covariates
#'   \item \code{linear.pred}: Linear predictor values from the fitted model
#' }
#' The returned object has class \code{c("zresid", "matrix")}.
#'
#' @examples
#' \dontrun{
#' # Compute Z-residuals for counts
#' zres_counts <- Zresidual.hurdle.negbinomial(fit, type = "count", method = "iscv")
#'
#' # Compute Z-residuals for the full hurdle model with 2 replicates
#' zres_hurdle <- Zresidual.hurdle.negbinomial(fit, type = "hurdle", n.rep = 2)
#' }
#'
#' @seealso
#' \code{\link{log_pred_dist_HNB}}, \code{\link{log_pred_dist_TNB}}, \code{\link{post_logrpp}}, \code{\link{iscv_logrpp}}
#' @export
Zresidual.hurdle.negbinomial <- function(fit, type , method = "iscv", n.rep = 1){

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
  type_list <- c("zero", "TNB", "HNB")
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
    linear.pred = predict(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
