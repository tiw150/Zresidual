log_pointpred_brmsfit_simulation <- function(fit, data, ...) {
  
  if (missing(data) || is.null(data)) {
    stop("log_pointpred_brmsfit_simulation: `data` must be provided.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data)
  model.var  <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  y_obs <- as.numeric(model.data[[response]])
  n <- length(y_obs)
  
  # Extract exact Log-Likelihood (S x n matrix)
  lpmf <- brms::log_lik(fit, newdata = data, ...)
  
  # Generate Simulated Responses (S x n matrix)
  y_sim <- brms::posterior_predict(fit, newdata = data, ...)
  
  # Calculate Point-wise Survival Indicator (S x n matrix)
  is_greater <- sweep(y_sim, 2, y_obs, FUN = ">")
  surv_indicator <- matrix(as.integer(is_greater), nrow = nrow(y_sim), ncol = n)
  
  # Apply log-scale: log(1) = 0, log(0) = -Inf
  lsf <- log(surv_indicator)
  
  # --- ADDED LINE: Extract the family from the fitted brms object ---
  fam <- family(fit)$family
  
  discrete_families <- c(
    "poisson", "negbinomial", "binomial", "bernoulli", "geometric", 
    "categorical", "multinomial", "beta_binomial", "discrete_weibull",
    "zero_inflated_poisson", "zero_inflated_negbinomial", "zero_inflated_binomial",
    "hurdle_poisson", "hurdle_negbinomial"
  )
  
  # Replicated to length n to ensure compatibility with downstream matrix checks
  is_discrete <- rep(as.integer(fam %in% discrete_families), n)
  
  list(
    log_like = lpmf,           # S x n matrix 
    log_surv = lsf,            # S x n matrix (0 and -Inf)
    is_discrete = is_discrete
  )
}
