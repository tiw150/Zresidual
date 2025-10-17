#' A function to calculate posterior parameter estimates
#'
#' @param fit A `brm` fit.
#' @param dpar parameter.
#'

.LINK_FUNCTIONS <- list(
  log = function(x) exp(x),  # If link is "log", apply exp()
  identity = function(x) x,  # Identity link, no transformation
  softplus = function(x) log(exp(x) + 1),  # Softplus
  logit = function(x) exp(x) / (1 + exp(x)),  # Logit
  inverse = function(x) 1 / x,  # Inverse
  cloglog = function(x) 1 - exp(-exp(x)),  # Complementary log-log
  sqrt = function(x) sqrt(x), # square-root
  probit = function(x) qnorm(x), # probit
  probit_approx = function(x) x / sqrt(1 + x^2),  # Probit approximation
  cauchit = function(x) qcauchy(x), # quantile function of Cauchy distribution
  inverse_squared = function(x) 1/(x^2),
  tan_half = function(x) tan(x/2),
  squareplus = function (x) (x + sqrt(x^2 + 1))/2
)

posterior.pred <- function(fit, dpar, count.only = TRUE, mc_used){

  link_functions <- .LINK_FUNCTIONS

  n <- dim(fit$data)[1]
  chains <- fit$fit@sim$chains
  iter <- fit$fit@sim$iter
  warmup <- fit$fit@sim$warmup
  mc_used <- chains*(iter - warmup)

  data <- fit$data
  response <- fit$formula$resp
  model.data <- model.frame(fit$formula, data=data)
  model.var <- names(model.data)
  if(!response %in% model.var) response <- strsplit(as.character(fit$formula$formula), "~")[[2]]
  data$shape <- 1
  data$hu <- 1
  sim.y <- as.matrix(model.data[, response])

  if(dpar != "zero" & count.only) id <- which(sim.y > 0) else id <- 1:n
  n.id <- length(id)
  data.id <- data
  data.id[-id,] <- NA

  formula_list <- list(
    zero = fit$formula$pforms$hu,
    mu = fit$formula$formula,
    shape = fit$formula$pforms$shape)

  link_name_list <- list(
    zero = fit$family$link_hu,
    mu = fit$family$link,
    shape = fit$family$link_shape)

  # Define the column name prefixes in MCMC parameter estimates
  para_prefix_list <- list(
    zero = "^b_hu_", # For hurdle portion
    mu = "^b_", #For mu
    shape = "^b_shape|shape") # For shape

  # Calculating parameter
  formula <- formula_list[[dpar]]
  if(!is.null(formula)){
    mm <- model.matrix(formula,data=data.id)
  } else{
    mm <- matrix(1, nrow = length(sim.y), ncol = 1)
    mm[-id,] <- NA
  }

  to_brmsnames_mm <- brms:::rename(colnames(mm))
  para_col <- paste0(para_prefix_list[[dpar]], gsub("\\[|\\]|\\(|\\)", "", to_brmsnames_mm))
  para <- as.data.frame(fit, variable = para_col, regex = TRUE)
  if(ncol(mm) != ncol(para)) stop("Model matrix and estimate column names do not match.")
  log_parameter <- tcrossprod(as.matrix(para), mm)

  link_name <- link_name_list[[dpar]]
  activation_func <- link_functions[[link_name]]
  if (!is.null(activation_func)) {parameter <- activation_func(log_parameter)
  }else{ stop(paste("Unknown link function:", link_name))}


  return(parameter)
}
