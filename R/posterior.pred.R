## This is the function for making posterior prediction for brms
posterior.pred <- function(fit, dpar, data, count.only = TRUE, ...) {
  .LINK_FUNCTIONS <- list(
    log = function(x) exp(x),
    identity = function(x) x,
    softplus = function(x) ifelse(x > 20, x, log1p(exp(x))),
    logit = function(x) plogis(x),
    inverse = function(x) 1 / x,
    cloglog = function(x) 1 - exp(-exp(x)),
    sqrt = function(x) x^2,
    probit = function(x) pnorm(x),
    probit_approx = function(x) x / sqrt(1 + x^2),
    cauchit = function(x) pcauchy(x),
    inverse_squared = function(x) 1 / sqrt(x),
    tan_half = function(x) tan(x / 2),
    squareplus = function(x) (x + sqrt(x^2 + 1)) / 2
  )
  
  if (missing(data) || is.null(data)) {
    stop("posterior.pred: `data` must be provided.", call. = FALSE)
  }
  
  data <- as.data.frame(data)
  n <- nrow(data)
  
  if (n < 1L) {
    stop("posterior.pred: `data` has zero rows.", call. = FALSE)
  }
  
  model.data <- model.frame(fit$formula, data = data, na.action = stats::na.pass)
  model.var <- names(model.data)
  
  response <- fit$formula$resp
  if (is.null(response) || !response %in% model.var) {
    response <- all.vars(fit$formula$formula)[1]
  }
  
  if (is.null(response) || !response %in% model.var) {
    stop("posterior.pred: response variable not found in `data`.", call. = FALSE)
  }
  
  data$shape <- 1
  data$hu <- 1
  sim.y <- as.matrix(model.data[, response, drop = FALSE])[, 1]
  
  if (dpar != "zero" && isTRUE(count.only)) {
    id <- which(sim.y > 0)
  } else {
    id <- seq_len(n)
  }
  
  data.id <- data
  if (length(id) < n) {
    data.id[-id, ] <- NA
  }
  
  formula_list <- list(
    zero = fit$formula$pforms$hu,
    mu = fit$formula$formula,
    shape = fit$formula$pforms$shape
  )
  
  link_name_list <- list(
    zero = fit$family$link_hu,
    mu = fit$family$link,
    shape = fit$family$link_shape
  )
  
  formula <- formula_list[[dpar]]
  
  if (!is.null(formula)) {
    mm <- model.matrix(formula, data = data.id)
    to_brmsnames_mm <- sanitize_names(colnames(mm))
    
    if (dpar == "zero") {
      para_col <- paste0("^b_hu_", gsub("\\[|\\]|\\(|\\)", "", to_brmsnames_mm))
    } else if (dpar == "mu") {
      para_col <- paste0("^b_", gsub("\\[|\\]|\\(|\\)", "", to_brmsnames_mm))
    } else if (dpar == "shape") {
      para_col <- paste0("^b_shape_", gsub("\\[|\\]|\\(|\\)", "", to_brmsnames_mm))
    } else {
      stop("posterior.pred: unsupported dpar.", call. = FALSE)
    }
  } else {
    mm <- matrix(1, nrow = n, ncol = 1)
    if (length(id) < n) {
      mm[-id, ] <- NA
    }
    colnames(mm) <- "(Intercept)"
    
    if (dpar == "zero") {
      para_col <- "^b_hu_Intercept$"
    } else if (dpar == "mu") {
      para_col <- "^b_Intercept$"
    } else if (dpar == "shape") {
      para_col <- "^(b_shape_Intercept|shape)$"
    } else {
      stop("posterior.pred: unsupported dpar.", call. = FALSE)
    }
  }
  
  para <- as.data.frame(fit, variable = para_col, regex = TRUE)
  
  if (ncol(para) == 0) {
    stop(
      paste0(
        "posterior.pred: no posterior columns matched for dpar='", dpar,
        "' using regex ", paste(para_col, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  if (ncol(mm) != ncol(para)) {
    stop(
      paste0(
        "posterior.pred: model matrix and estimate columns do not match for dpar='",
        dpar, "'. ncol(mm)=", ncol(mm), ", ncol(para)=", ncol(para)
      ),
      call. = FALSE
    )
  }
  
  log_parameter <- tcrossprod(as.matrix(para), mm)
  
  link_name <- link_name_list[[dpar]]
  activation_func <- .LINK_FUNCTIONS[[link_name]]
  
  if (is.null(activation_func)) {
    stop(paste("Unknown link function:", link_name), call. = FALSE)
  }
  
  activation_func(log_parameter)
}