## This is the function for making posterior parameter predictions for brms
posterior.pred <- function(
    fit,
    dpar = "mu",
    data,
    count.only = TRUE,
    re_formula = NA,
    eps = .Machine$double.eps,
    ...
) {
  if (!inherits(fit, "brmsfit")) {
    stop("posterior.pred: `fit` must be a brmsfit object.", call. = FALSE)
  }
  
  if (missing(data) || is.null(data)) {
    stop("posterior.pred: `data` must be provided.", call. = FALSE)
  }
  
  data <- as.data.frame(data)
  n <- nrow(data)
  
  if (n < 1L) {
    stop("posterior.pred: `data` has zero rows.", call. = FALSE)
  }
  
  if (missing(dpar) || is.null(dpar) || length(dpar) != 1L) {
    stop("posterior.pred: `dpar` must be a single character value.", call. = FALSE)
  }
  
  ## Backward-compatible aliases used in the old package code.
  ## zero = hurdle probability in brms, usually called hu.
  dpar_alias <- list(
    zero = "hu",
    hurdle = "hu",
    zi_prob = "zi",
    zoi_prob = "zoi",
    coi_prob = "coi"
  )
  
  brms_dpar <- if (dpar %in% names(dpar_alias)) {
    dpar_alias[[dpar]]
  } else {
    dpar
  }
  
  ## brms main response parameter:
  ## For the main mean/location parameter, posterior_linpred usually works
  ## without specifying dpar. Some brms versions also allow dpar = "mu".
  parameter <- tryCatch(
    {
      if (brms_dpar %in% c("mu", "mean", "location")) {
        out <- try(
          brms::posterior_linpred(
            fit,
            newdata = data,
            dpar = "mu",
            transform = TRUE,
            re_formula = re_formula,
            ...
          ),
          silent = TRUE
        )

        if (inherits(out, "try-error")) {
          out <- brms::posterior_linpred(
            fit,
            newdata = data,
            transform = TRUE,
            re_formula = re_formula,
            ...
          )
        }

        out
      } else {
        brms::posterior_linpred(
          fit,
          newdata = data,
          dpar = brms_dpar,
          transform = TRUE,
          re_formula = re_formula,
          ...
        )
      }
    },
    error = function(e) {
      stop(
        paste0(
          "posterior.pred: failed to extract posterior predictions for dpar='",
          dpar,
          "' (mapped to brms dpar='",
          brms_dpar,
          "'). Original error: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  
  parameter <- as.matrix(parameter)
  
  if (ncol(parameter) != n) {
    stop(
      paste0(
        "posterior.pred: output dimension mismatch for dpar='",
        dpar,
        "'. ncol(parameter)=",
        ncol(parameter),
        ", nrow(data)=",
        n
      ),
      call. = FALSE
    )
  }
  
  if (any(!is.finite(parameter))) {
    stop(
      paste0(
        "posterior.pred: non-finite posterior predictions returned for dpar='",
        dpar,
        "' (mapped to brms dpar='",
        brms_dpar,
        "')."
      ),
      call. = FALSE
    )
  }
  
  ## Numerical guardrails by parameter type.
  ## Positive parameters.
  positive_dpars <- c(
    "mu", "shape", "sigma", "phi", "nu", "alpha",
    "disc", "bs", "ndt", "bias", "rate", "scale"
  )
  
  ## Probability parameters.
  probability_dpars <- c(
    "hu", "zi", "zoi", "coi", "theta", "prob"
  )
  
  if (brms_dpar %in% positive_dpars) {
    parameter <- pmax(parameter, eps)
  }

  if (brms_dpar %in% probability_dpars) {
    parameter <- pmin(pmax(parameter, eps), 1 - eps)
  }

  ## Keep old behavior:
  ## for hurdle count components, optionally return NA for zero observations.
  ## This should only be applied to count-side parameters, not to hu/zi/etc.
  if (isTRUE(count.only) && !(brms_dpar %in% c("hu", "zi", "zoi", "coi"))) {
    response <- fit$formula$resp

    if (is.null(response) || !response %in% names(data)) {
      response <- all.vars(fit$formula$formula)[1]
    }

    if (!is.null(response) && response %in% names(data)) {
      y <- data[[response]]

      ## Only apply this for count / hurdle style responses.
      ## Do not blindly mask Gaussian/Beta/binomial data.
      fam <- fit$family$family

      hurdle_or_count_family <- fam %in% c(
        "hurdle_poisson",
        "hurdle_negbinomial",
        "hurdle_gamma",
        "hurdle_lognormal"
      )

      if (hurdle_or_count_family) {
        parameter[, y <= 0] <- NA_real_
      }
    }
  }
  
  parameter
}
