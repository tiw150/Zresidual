#' A function to calculate z-residuals of a 'brm' fit.
#' This function is to be used when the user needs to calculate the z-residuals of bernoulli
#'
#' @param fit A `brm` fit.
#' @param type The component which the z-residuals should be calculated.
#' @export
#'
zresidual_bernoulli <- function(fit, method = "iscv", nrep = 1){

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

  ldist <- log.pred.dist.bern(fit)
  lpmf <- ldist$lpmf_hat
  lcdf <- ldist$lcdf_hat

  # Argument should be one of the element in rpp_list
  rpp_list <- c(iscv = "iscv_logrpp", post = "post_logrpp")
  names(rpp_list) <- c("iscv", "post")

  #if(count_only) z_res<- matrix(NA, ncol = nrep, nrow = dim(lpmf)[2])
  z_res<- matrix(NA, ncol = nrep, nrow = n)

  for (i in 1:nrep) {
    rpp <- get(rpp_list[[method]])(lcdf, lpmf)
    z_res[, i] <- -qnorm(rpp, log.p = T)
  }

  #if(type == "count") z_res <- z_res[-zero_id,]

  colnames(z_res) <- paste0("zresidual", 1:nrep)

  attributes(z_res) <- c(attributes(z_res),list(
    type = type,
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
