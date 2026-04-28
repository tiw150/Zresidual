#' Compute Z-residuals from predictive tail probabilities
#'
#' Computes Z-residuals by transforming predictive upper-tail probabilities to
#' the standard normal scale. When point-wise predictive mass at the observed value
#' is available, randomized Z-residuals are constructed by adding a uniform
#' perturbation within that mass. This allows the same interface to handle both
#' continuous-type and discrete-type outcomes.
#' 
#' To avoid repeating expensive MCMC summaries across multiple randomizations,
#' users can precompute the predictive summaries using \code{\link{log_summary_pred}} 
#' and pass the result to the \code{log_sumpred} argument.
#'
#' @param fit A fitted model object. Optional if \code{log_sumpred} is provided.
#' @param data A data frame or list containing the data used to evaluate predictive 
#'   quantities. Optional if \code{log_sumpred} is provided.
#' @param log_pointpred Optional function (or function name as a character string) 
#'   returning predictive pointwise quantities. If \code{NULL}, the function name 
#'   is automatically constructed based on the model backend, family, and \code{pred_method}.
#' @param pred_method Character string indicating the prediction method. Must be one of 
#'   \code{"analytic"} or \code{"simulation"}.
#' @param mcmc_summarize Character string indicating how posterior predictive
#'   quantities are summarized. Must be one of \code{"post"} (posterior mean) or 
#'   \code{"iscv"} (importance sampling cross-validation).
#' @param type Optional character string used as a component selector or variant tag.
#' @param randomized Logical; if \code{TRUE}, generate randomized Z-residuals
#'   when point mass or discrete-type contribution is available.
#' @param nrep Integer giving the number of residual replicates to generate. If
#'   \code{randomized = FALSE}, this is forced to \code{1L}.
#' @param eps Small positive numeric constant used to keep probability values 
#'   strictly bounded away from exactly 0 and 1.
#' @param seed Optional integer seed for reproducible residual randomization.
#' @param log_sumpred Optional precomputed output list from \code{\link{log_summary_pred}} 
#'   to bypass MCMC summarization.
#' @param ... Additional arguments passed to the underlying prediction backend.
#'
#' @return A numeric matrix with one row per observation and \code{nrep} columns containing 
#'   the computed Z-residual replicates. The returned object has class \code{"zresid"}.
#'   Attribute \code{"rsp"} stores the probability-scale residuals used to form
#'   the Z-residuals. Other contextual attributes (e.g., covariates, linear predictors) 
#'   may be attached if \code{Zcov} successfully extracts them.
#' 
#' @seealso \code{\link{log_summary_pred}}
#' @export
#' 
Zresidual <- function(fit = NULL,
                      data = NULL,
                      log_pointpred = NULL,
                      pred_method = c("analytic", "simulation"),
                      mcmc_summarize = c("post", "iscv"),
                      type = NULL,
                      randomized = TRUE,
                      nrep = 30,
                      eps = 1e-12,
                      seed = NULL,
                      log_sumpred = NULL,
                      ...) {
  if (!is.null(seed)) set.seed(seed)
  
  .log_add_exp2 <- function(a, b) {
    both_ninf <- is.infinite(a) & a < 0 & is.infinite(b) & b < 0
    m <- pmax(a, b)
    res <- m + log(exp(a - m) + exp(b - m))
    res[both_ninf] <- -Inf
    
    a_ninf <- is.infinite(a) & a < 0 & !both_ninf
    res[a_ninf] <- b[a_ninf]
    
    b_ninf <- is.infinite(b) & b < 0 & !both_ninf
    res[b_ninf] <- a[b_ninf]
    
    res
  }
  
  .clip_logprob <- function(logp, eps = 1e-100) {
    lo <- log(eps)
    hi <- log1p(-eps)
    pmin(pmax(logp, lo), hi)
  }
  
  # Allow bypass of expensive MCMC summary extraction
  if (!is.null(log_sumpred)) {
    pre <- log_sumpred
  } else {
    if (missing(data) || is.null(data)) {
      stop("Zresidual: `data` must be provided unless `log_sumpred` is used.", call. = FALSE)
    }
    pred_method <- match.arg(tolower(as.character(pred_method)), c("analytic", "simulation"))
    mcmc_summarize <- match.arg(tolower(as.character(mcmc_summarize)), c("post", "iscv"))
    
    # --- ISCV Sample Size Check ---
    if (mcmc_summarize == "iscv" && !is.null(fit)) {
      n_fit <- tryCatch(nobs(fit), error = function(e) NULL)
      n_data <- NROW(data)
      
      if (!is.null(n_fit) && n_fit != n_data) {
        warning(
          sprintf("Zresidual: 'iscv' requested but fit sample size (%d) differs from data sample size (%d). Switching to 'post'.", n_fit, n_data),
          call. = FALSE
        )
        mcmc_summarize <- "post"
      }
    }
    # ------------------------------
    
    pre <- log_summary_pred(
      fit = fit,
      data = data,
      log_pointpred = log_pointpred,
      pred_method = pred_method,
      mcmc_summarize = mcmc_summarize,
      type = type,
      ...
    )
  }
  
  log_surv_hat <- as.numeric(pre$log_surv_hat)
  log_pmf_hat  <- as.numeric(pre$log_pmf_hat)
  is_discrete  <- as.integer(pre$is_discrete)
  
  n <- length(log_surv_hat)
  if (n < 1L) stop("Zresidual: empty prediction.", call. = FALSE)
  if (length(log_pmf_hat) != n) stop("Zresidual: `log_pmf_hat` length mismatch.", call. = FALSE)
  if (length(is_discrete) != n) stop("Zresidual: `is_discrete` length mismatch.", call. = FALSE)
  
  idx_disc <- which(is_discrete == 1L)
  if (length(idx_disc) > 0L) {
    bad <- idx_disc[!is.finite(log_pmf_hat[idx_disc])]
    if (length(bad) > 0L) {
      stop("Zresidual: discrete observations require finite `log_pmf_hat`.", call. = FALSE)
    }
  }
  
  if (!randomized) nrep <- 1L
  nrep <- as.integer(nrep)
  if (nrep < 1L) stop("Zresidual: `nrep` must be >= 1.", call. = FALSE)
  
  log_rsp <- matrix(NA_real_, nrow = n, ncol = nrep)
  
  for (r in seq_len(nrep)) {
    u <- if (randomized) runif(n) else rep(0.5, n)
    logu <- log(pmax(u, .Machine$double.xmin))
    
    lr <- log_surv_hat
    if (length(idx_disc) > 0L) {
      lr[idx_disc] <- .log_add_exp2(
        log_surv_hat[idx_disc],
        log_pmf_hat[idx_disc] + logu[idx_disc]
      )
    }
    
    log_rsp[, r] <- .clip_logprob(lr, eps = eps)
  }
  
  rsp <- exp(log_rsp)
  z <- -qnorm(log_rsp, log.p = TRUE) ## more robust than plugging rsp
  
  if (!is.matrix(z)) z <- matrix(z, nrow = n, ncol = nrep)
  colnames(z) <- paste0("rep", seq_len(nrep))
  
  class(z) <- c("zresid", class(z))
  attr(z, "rsp") <- rsp
  
  # Only attempt covariate extraction if fit/data were provided
  .type_from_zcov <- function(zcov) {
    if (is.null(zcov)) return(NULL)
    if (!is.null(zcov$type) && nzchar(as.character(zcov$type)[1])) {
      return(as.character(zcov$type)[1])
    }
    kind <- zcov$y_type_kind
    if (is.null(kind) || length(kind) < 1L) return(NULL)
    kind <- as.character(kind[1])
    if (identical(kind, "censor")) return("survival")
    if (identical(kind, "trunc"))  return("count")
    if (identical(kind, "hurdle")) return("hurdle")
    NULL
  }
  
  if (!is.null(fit) && !is.null(data)) {
    zcov <- tryCatch(
      Zcov(fit, data = data, type = type),
      error = function(e) NULL
    )
    
    if (!is.null(zcov)) {
      attr(z, "zcov") <- zcov
      attr(z, "linear.pred") <- zcov$linear_pred
      attr(z, "linear_pred") <- zcov$linear_pred
      attr(z, "covariates") <- zcov$covariates
      if (is.null(attr(z, "type"))) attr(z, "type") <- .type_from_zcov(zcov)
      
      if (identical(zcov$y_type_kind, "censor") && !is.null(zcov$y_type)) {
        attr(z, "censored.status") <- as.integer(zcov$y_type == 0L)
      }
      if (!is.null(zcov$extra) && !is.null(zcov$extra$zero_id)) {
        attr(z, "zero_id") <- zcov$extra$zero_id
      }
    }
  }
  
  z
}

#' Extract and summarize MCMC predictive tail probabilities
#'
#' Evaluates observation-level predictive log-survival probabilities and log-probability 
#' masses (for discrete data) over the MCMC draws, and aggregates them into a single 
#' summary vector per observation based on the chosen \code{mcmc_summarize} method. 
#' This function isolates the computationally expensive matrix operations 
#' (such as log-sum-exp) from the randomization step.
#'
#' @param fit A fitted model object.
#' @param data A data frame or list containing the data used to evaluate predictive 
#'   quantities. Must be provided.
#' @param log_pointpred Optional function (or function name as a character string) 
#'   returning predictive pointwise quantities. If \code{NULL}, the function is automatically 
#'   resolved based on the model backend, family, and \code{pred_method} using 
#'   \code{\link{required_log_pointpred}}. The resolved function must return a list 
#'   containing \code{log_surv}, \code{log_like}, and \code{is_discrete}.
#' @param pred_method Character string indicating the prediction method. Must be one of 
#'   \code{"analytic"} or \code{"simulation"}.
#' @param mcmc_summarize Character string indicating how posterior predictive
#'   quantities are summarized. Must be one of \code{"post"} (posterior mean) or 
#'   \code{"iscv"} (importance sampling cross-validation).
#' @param type Optional character string used as a component selector or variant tag.
#' @param ... Additional arguments passed to the underlying prediction backend.
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{\code{log_surv_hat}}{A numeric vector of the summarized log-survival probabilities.}
#'   \item{\code{log_pmf_hat}}{A numeric vector of the summarized log-probability masses (or \code{NA} for continuous observations).}
#'   \item{\code{is_discrete}}{An integer vector indicating whether each observation is treated as discrete (\code{1}) or continuous (\code{0}).}
#' }
#' 
#' @seealso \code{\link{Zresidual}}, \code{\link{required_log_pointpred}}
#' @export
#' 
log_summary_pred <- function(fit,
                             data,
                             log_pointpred = NULL,
                             pred_method = c("analytic", "simulation"),
                             mcmc_summarize = c("post", "iscv"),
                             type = NULL,
                             ...) {
  if (missing(data) || is.null(data)) {
    stop("log_summary_pred: `data` must be provided.", call. = FALSE)
  }
  
  pred_method <- match.arg(tolower(as.character(pred_method)), c("analytic", "simulation"))
  mcmc_summarize <- match.arg(tolower(as.character(mcmc_summarize)), c("post", "iscv"))
  
  # Resolve log_pointpred or fallback to automatic discovery
  if (!is.null(log_pointpred)) {
    if (is.function(log_pointpred)) {
      lp_fn <- log_pointpred
    } else if (is.character(log_pointpred) && length(log_pointpred) == 1L) {
      lp_fn <- match.fun(log_pointpred)
    } else {
      stop("log_summary_pred: `log_pointpred` must be a function or a character string.", call. = FALSE)
    }
  } else {
    lp_fn <- required_log_pointpred(
      fit = fit, 
      pred_method = pred_method, 
      type = type
    )
    
    # If required_log_pointpred returns a string, resolution failed. Carry the message into stop().
    if (is.character(lp_fn)) {
      stop(lp_fn, call. = FALSE)
    }
  }

  as_mat <- function(z, name) {
    if (is.null(z)) stop("log_summary_pred: missing `", name, "`.", call. = FALSE)
    if (is.vector(z)) return(matrix(as.vector(z), nrow = 1L))
    if (is.matrix(z)) return(z)
    stop("log_summary_pred: `", name, "` must be vector or matrix.", call. = FALSE)
  }
  
  safe_col_logmeanexp <- function(mat) {
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      stop("log_summary_pred: Bayesian summaries require package 'matrixStats'.", call. = FALSE)
    }
    Tdraw <- nrow(mat)
    out <- matrixStats::colLogSumExps(mat, na.rm = TRUE) - log(Tdraw)
    out[is.nan(out) | is.na(out)] <- -Inf
    out
  }
  
  post_hat <- function(log_surv_m, log_like_m, is_discrete) {
    n <- ncol(log_surv_m)
    log_pmf_hat <- rep(NA_real_, n)
    idx_disc <- which(is_discrete == 1L)
    if (length(idx_disc) > 0L) {
      log_pmf_hat[idx_disc] <- safe_col_logmeanexp(log_like_m[, idx_disc, drop = FALSE])
    }
    list(
      log_surv_hat = safe_col_logmeanexp(log_surv_m),
      log_pmf_hat = log_pmf_hat
    )
  }
  
  iscv_hat <- function(log_surv_m, log_like_m, is_discrete) {
    if (!requireNamespace("matrixStats", quietly = TRUE)) {
      stop("log_summary_pred: Bayesian summaries require package 'matrixStats'.", call. = FALSE)
    }
    if (anyNA(log_like_m)) {
      stop("log_summary_pred: `iscv` requires non-NA `log_like`.", call. = FALSE)
    }
    
    all_zero_lik <- apply(log_like_m, 2, function(v) all(is.infinite(v) & v < 0))
    if (any(all_zero_lik)) {
      stop("log_summary_pred: `iscv` is not compatible with columns whose likelihood is identically zero.", call. = FALSE)
    }
    
    log_sum_S_over_lik <- matrixStats::colLogSumExps(log_surv_m - log_like_m)
    log_sum_inv_lik <- matrixStats::colLogSumExps(-log_like_m)
    log_surv_hat <- log_sum_S_over_lik - log_sum_inv_lik
    
    n <- ncol(log_surv_m)
    log_pmf_hat <- rep(NA_real_, n)
    idx_disc <- which(is_discrete == 1L)
    if (length(idx_disc) > 0L) {
      log_pmf_hat[idx_disc] <- log(nrow(log_like_m)) - log_sum_inv_lik[idx_disc]
    }
    
    list(
      log_surv_hat = log_surv_hat,
      log_pmf_hat = log_pmf_hat
    )
  }
  
  call_log_pointpred <- function(fun, fit, data, type = NULL, ...) {
    fmls <- tryCatch(names(formals(fun)), error = function(e) character(0))
    accepts_dots <- "..." %in% fmls
    accepts_type <- "type" %in% fmls
    if (accepts_type || accepts_dots) {
      return(fun(fit, data = data, type = type, ...))
    }
    fun(fit, data = data, ...)
  }
  
  pp <- call_log_pointpred(lp_fn, fit = fit, data = data, type = type, ...)
  
  if (!is.list(pp) || is.null(pp$log_surv) || is.null(pp$log_like) || is.null(pp$is_discrete)) {
    stop("log_summary_pred: log_pointpred must return `log_surv`, `log_like`, and `is_discrete`.", call. = FALSE)
  }
  
  log_surv0 <- as_mat(pp$log_surv, "log_surv")
  log_like0 <- as_mat(pp$log_like, "log_like")
  isdisc0   <- as_mat(pp$is_discrete, "is_discrete")
  
  if (ncol(log_surv0) != ncol(log_like0)) {
    stop("log_summary_pred: `log_surv` and `log_like` dimension mismatch.", call. = FALSE)
  }
  if (nrow(log_surv0) != nrow(log_like0)) {
    stop("log_summary_pred: `log_surv` and `log_like` draw dimension mismatch.", call. = FALSE)
  }
  if (nrow(isdisc0) != 1L || ncol(isdisc0) != ncol(log_surv0)) {
    stop("log_summary_pred: `is_discrete` must be length n or a 1 x n matrix.", call. = FALSE)
  }
  
  is_discrete <- as.integer(isdisc0[1, ])
  is_bayes <- nrow(log_surv0) > 1L
  
  if (!is_bayes) {
    log_pmf_hat <- rep(NA_real_, ncol(log_surv0))
    idx_disc <- which(is_discrete == 1L)
    if (length(idx_disc) > 0L) {
      log_pmf_hat[idx_disc] <- as.numeric(log_like0[1L, idx_disc])
    }
    return(list(
      log_surv_hat = as.numeric(log_surv0[1L, ]),
      log_pmf_hat = log_pmf_hat,
      is_discrete = is_discrete
    ))
  }
  
  hat <- switch(
    mcmc_summarize,
    post = post_hat(log_surv0, log_like0, is_discrete),
    iscv = iscv_hat(log_surv0, log_like0, is_discrete)
  )
  
  list(
    log_surv_hat = as.numeric(hat$log_surv_hat),
    log_pmf_hat  = as.numeric(hat$log_pmf_hat),
    is_discrete  = is_discrete
  )
}

#' Resolve or Identify Required point-wise predictive Functions
#'
#' @description
#' Determines which native point-wise predictive function should be utilized 
#' for evaluating predictive tail probabilities. It dynamically searches the calling 
#' environment, global environment, and package namespace. 
#' 
#' This function evaluates multiple naming conventions by extracting the model's 
#' package namespace, S3 class, and statistical family. If no suitable function 
#' is found (or if \code{show_names} is \code{TRUE}), it compiles a descriptive 
#' message detailing the exact function names expected by the package.
#'
#' @param fit A fitted model object. Used to automatically extract the model 
#'   package, class, and family.
#' @param pred_method Character string indicating the prediction method. Must be 
#'   one of \code{"analytic"} or \code{"simulation"}.
#' @param type Optional character string used as a component selector or variant tag.
#' @param show_names Logical; if \code{TRUE}, bypasses the environment search, 
#'   prints the compiled string of expected function names to the console, and 
#'   returns it invisibly. Defaults to \code{FALSE}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The resolved function object. If no suitable function is found, it returns 
#'   a character string containing the compiled diagnostic message. If 
#'   \code{show_names = TRUE}, it prints the message and returns the string invisibly.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' # Example 1: Standard GLM
#' fit_glm <- glm(vs ~ mpg, data = mtcars, family = binomial())
#' required_log_pointpred(fit_glm, show_names = TRUE)
#' 
#' # Example 2: Survival Model (coxph)
#' library(survival)
#' fit_cox <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' required_log_pointpred(fit_cox, show_names = TRUE)
#' 
#' # Example 3: Bayesian Model (brms) using pre-compiled package template
#' template_path <- system.file("extdata", "brms_template.rds", package = "Zresidual")
#' if (file.exists(template_path)) {
#'   fit_brm <- readRDS(template_path)
#'   required_log_pointpred(fit_brm, show_names = TRUE)
#' }
#' }
required_log_pointpred <- function(fit = NULL, 
                                   pred_method = c("analytic", "simulation"), 
                                   type = NULL, 
                                   show_names = (sys.nframe() == 1),
                                   ...) {
  
  pred_method <- match.arg(pred_method)
  
  # --- 1. Safely Extract Model Attributes ---
  cls <- if (!is.null(fit)) class(fit)[1] else "UNKNOWN"
  
  pkg <- "UNKNOWN"
  if (cls != "UNKNOWN") {
    pred_fn <- tryCatch(utils::getS3method("predict", cls, optional = TRUE), error = function(e) NULL)
    if (!is.null(pred_fn)) {
      env_name <- environmentName(environment(pred_fn))
      if (nzchar(env_name) && !env_name %in% c("R_GlobalEnv", "base")) {
        pkg <- env_name
      }
    }
    if (pkg == "UNKNOWN" && !is.null(attr(class(fit), "package"))) {
      pkg <- attr(class(fit), "package")[1]
    }
  }
  
  fam <- if (!is.null(fit)) {
    res <- tryCatch(family(fit)$family, error = function(e) "UNKNOWN", warning = function(w) "UNKNOWN")
    if (is.null(res) || length(res) == 0 || !is.character(res)) "UNKNOWN" else res[1]
  } else {
    "UNKNOWN"
  }
  
  type_suffix <- if (!is.null(type)) paste0("_", type) else ""
  
  # --- 2. Construct Ordered Candidate Names ---
  # Build the exhaustive list of names we will look for, from most specific to least specific.
  cands_analytic <- character(0)
  cands_sim      <- character(0)
  
  if (cls != "UNKNOWN") {
    base_pkg_cls <- if (pkg != "UNKNOWN") paste0("log_pointpred_", pkg, "_", cls) else NULL
    base_cls     <- paste0("log_pointpred_", cls)
    
    # Analytic Candidates
    if (fam != "UNKNOWN") {
      if (!is.null(base_pkg_cls)) cands_analytic <- c(cands_analytic, paste0(base_pkg_cls, "_", fam, type_suffix))
      cands_analytic <- c(cands_analytic, paste0(base_cls, "_", fam, type_suffix))
    }
    if (!is.null(base_pkg_cls)) cands_analytic <- c(cands_analytic, paste0(base_pkg_cls, type_suffix))
    cands_analytic <- c(cands_analytic, paste0(base_cls, type_suffix))
    
    # Simulation Candidates
    if (!is.null(base_pkg_cls)) cands_sim <- c(cands_sim, paste0(base_pkg_cls, "_simulation", type_suffix))
    cands_sim <- c(cands_sim, paste0(base_cls, "_simulation", type_suffix))
  }
  
  # Select the active candidate list based on user request
  active_cands <- if (pred_method == "analytic") c(cands_analytic, cands_sim) else cands_sim
  
  # --- 3. Internal Helper: Ordered Environment Search ---
  find_func <- function(fn_name) {
    # 1. Calling Environment (User's script)
    out <- get0(fn_name, envir = parent.frame(n = 3), inherits = TRUE)
    if (is.function(out)) return(out)
    
    # 2. Global Environment
    out <- get0(fn_name, envir = .GlobalEnv, inherits = FALSE)
    if (is.function(out)) return(out)
    
    # 3. Zresidual Namespace (Native support)
    out <- get0(fn_name, envir = asNamespace("Zresidual"), inherits = FALSE)
    if (is.function(out)) return(out)
    
    # 4. Model's Package Namespace (If a 3rd party package exports Zresidual support)
    if (pkg != "UNKNOWN") {
      out <- tryCatch(get0(fn_name, envir = asNamespace(pkg), inherits = FALSE), error = function(e) NULL)
      if (is.function(out)) return(out)
    }
    
    return(NULL)
  }
  
  # --- 4. Execution Loop ---
  if (!show_names && length(active_cands) > 0) {
    for (cand in active_cands) {
      fn <- find_func(cand)
      if (!is.null(fn)) return(fn)
    }
  }
  
  # --- 5. Resolution Failed / show_names formatting ---
  msg <- if (show_names) {
    "Expected point-wise predictive function names:\n\n"
  } else {
    "log_summary_pred: The point-wise predictive function could not be resolved.\n\n"
  }
  
  if (cls == "UNKNOWN") {
    if (!show_names) msg <- paste0(msg, "--- ACTION REQUIRED ---\n")
    msg <- paste0(msg, "Because the model class could not be determined (or `fit` is NULL),\n")
    msg <- paste0(msg, "the expected function name cannot be automatically constructed.\n\n")
    msg <- paste0(msg, "You MUST pass your custom function directly to the `log_pointpred` argument.\n")
    msg <- paste0(msg, "Example: Zresidual(fit, data, log_pointpred = my_custom_function)\n\n")
  } else {
    if (!show_names) msg <- paste0(msg, "--- EXPECTED FUNCTION NAMES ---\n")
    msg <- paste0(msg, "To proceed, you must define a function in your environment using one of the exact\n")
    msg <- paste0(msg, "names below, OR pass your function directly to the `log_pointpred` argument.\n\n")
    
    msg <- paste0(msg, sprintf("Parsed Model Specs -> Package: '%s' | Class: '%s' | Family: '%s'\n\n", pkg, cls, fam))
    
    if (length(cands_analytic) > 0) {
      msg <- paste0(msg, "[Analytic Methods]\n")
      for (cand in cands_analytic) {
        msg <- paste0(msg, sprintf("  - %s\n", cand))
      }
      msg <- paste0(msg, "\n")
    }
    
    if (length(cands_sim) > 0) {
      msg <- paste0(msg, "[Simulation Methods]\n")
      for (cand in cands_sim) {
        msg <- paste0(msg, sprintf("  - %s\n", cand))
      }
      msg <- paste0(msg, "\n")
    }
  }
  
  # Detail the Function Structure
  msg <- paste0(msg, "--- REQUIRED FUNCTION STRUCTURE ---\n")
  msg <- paste0(msg, "Your defined function must have the following signature:\n\n")
  msg <- paste0(msg, "  Inputs:\n")
  msg <- paste0(msg, "    fit  : The fitted model object.\n")
  msg <- paste0(msg, "    data : The dataset used for evaluation.\n")
  msg <- paste0(msg, "    type : (Optional) Component selector or variant tag.\n")
  msg <- paste0(msg, "    ...  : Additional arguments passed from Zresidual.\n\n")
  msg <- paste0(msg, "  Output (A named list containing):\n")
  msg <- paste0(msg, "    log_surv    : Matrix (M draws x N obs) of log-survival probabilities.\n")
  msg <- paste0(msg, "    log_like    : Matrix (M draws x N obs) of log-likelihoods / log-probability masses.\n")
  msg <- paste0(msg, "    is_discrete : Integer vector (length N) where 1 = discrete, 0 = continuous.\n")
  
  if (show_names) {
    cat(msg)
    return(invisible(msg))
  }
  
  return(msg)
}
