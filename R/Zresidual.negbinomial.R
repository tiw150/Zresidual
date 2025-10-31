#' Compute Z-Residuals for Negative Binomial Models
#'
#' Computes Z-residuals for fitted Bayesian negative binomial models.
#' Z-residuals are useful for model diagnostics, including checking fit and
#' overdispersion, and can be calculated using posterior or cross-validated predictive p-values.
#'
#' @param fit A fitted \pkg{brms} model object for a negative binomial outcome.
#' @param method Character string specifying the residual calculation method:
#'   \code{"iscv"} for importance-sampled cross-validated randomized predictive p-values,
#'   \code{"rpost"} for posterior predictive p-values, or
#'   \code{"mpost"} for marginal posterior predictive p-values. Default is \code{"iscv"}.
#' @param n.rep Integer; the number of replicated Z-residual sets to generate. Default is 1.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the observed response vector from the model data.
#'   \item Computes the log-PMF and log-CDF for the negative binomial model
#'         using \code{\link{log.pred.dist.NB}}.
#'   \item Generates randomized or posterior predictive p-values according to the
#'         specified \code{method}.
#'   \item Converts the p-values to Z-residuals via the negative quantile of the
#'         standard normal distribution.
#' }
#'
#' The output is a matrix of Z-residuals with one column per replication.
#'
#' @return A numeric matrix of Z-residuals with attributes:
#' \itemize{
#'   \item \code{type}: The specified model component (currently NULL)
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
#' # Compute Z-residuals for a negative binomial model
#' zres_nb <- Zresidual.negbinomial(fit, method = "iscv")
#'
#' # Compute Z-residuals with 2 replicates using posterior predictive p-values
#' zres_nb_post <- Zresidual.negbinomial(fit, method = "rpost", n.rep = 2)
#' }
#'
#' @seealso
#' \code{\link{log.pred.dist.NB}}, \code{\link{post_logrpp}}, \code{\link{iscv_logrpp}}
#' @export
Zresidual.negbinomial <- function(fit, method = "iscv", n.rep = 1){

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
  #type_list <- c("NB")
  #names(type_list) <- c("nb")

  ldist <- log.pred.dist.NB(fit)
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
    type = NULL,
    #count_only = count_only,
    zero_id = zero_id,
    log_pmf = lpmf,
    log_cdf = lcdf,
    covariates = subset(fit$data, select = -get(response)),
    linear.pred = predict(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
