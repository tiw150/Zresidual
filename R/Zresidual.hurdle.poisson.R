#' A function to calculate z-residuals of a 'brm' fit.
#' This function is to be used when the user needs to calculate the Z-residuals of TP/HP
#'
#' @param fit is the fit object are from 'brms' package.
#'
#' @export
#'
#' @return \itemize{
#'  \item{Zresid} {Z-residual}
#'  }
#'
Zresidual.hurdle.poisson <- function(fit,  type , method = "iscv", nrep = 1){

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

  ldist <- get(paste0("log.pred.dist.", type_list[type]))(fit)
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
    linear.pred = fitted(fit, type = "conditional")[,1]
  ))

  class(z_res) <- c("zresid", class(z_res))

  return(z_res)
}
