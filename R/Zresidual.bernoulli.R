#' Compute Z-Residuals for a Bernoulli/Logistic Model
#'
#' Computes Z-residuals for a fitted Bayesian Bernoulli/Logistic (binary) model.
#' Z-residuals are calculated using posterior predictive methods
#' and can be used for model diagnostics.
#'
#' @param fit A fitted \pkg{brms} model object.
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
#'   \item Computes the log-PMF and log-CDF for the Bernoulli model using
#'         \code{\link{log_pred_dist_bern}}.
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
#'   \item \code{type}: Type of outcome (Bernoulli)
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
#' # Compute Z-residuals using ISCV method
#' zres <- Zresidual.bernoulli(fit, method = "iscv", n.rep = 2)
#'
#' # Compute Z-residuals using posterior predictive method
#' zres_post <- Zresidual.bernoulli(fit, method = "post")
#' }
#'
#' @seealso
#' \code{\link{log_pred_dist_bern}}, \code{\link{post_logrpp}}, \code{\link{iscv_logrpp}}
#' @export
Zresidual.bernoulli <- function(fit, method = "iscv", n.rep = 1){

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
