get_logpre <- function(fit = NULL,
                       data = NULL,
                       log_pointpred = NULL,
                       logpred = NULL,
                       config = NULL,
                       ...) {
  
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  if (is.null(fit)) stop("get_logpre: `fit` is required.", call. = FALSE)
  if (!is.null(log_pointpred) && !is.null(logpred)) {
    stop("get_logpre: provide only one of `log_pointpred` or `logpred`.", call. = FALSE)
  }
  
  # resolve function input: function object or single name string
  resolve_fun <- function(f) {
    if (is.null(f)) return(NULL)
    if (is.function(f)) return(f)
    if (is.character(f) && length(f) == 1L) {
      out <- get0(f, envir = parent.frame(), inherits = TRUE)
      if (is.null(out)) out <- get0(f, envir = asNamespace(utils::packageName()), inherits = FALSE)
      if (!is.function(out)) stop("get_logpre: cannot resolve function: ", f, call. = FALSE)
      return(out)
    }
    stop("get_logpre: function must be a function or a single character name.", call. = FALSE)
  }
  
  config <- config %||% list()
  
  # Case A: user provides log_pointpred (already standardized)
  if (!is.null(log_pointpred)) {
    lp_fn <- resolve_fun(log_pointpred)
    return(log_summary_pred(fit, data = data, config = config, log_pointpred_fn = lp_fn, ...))
  }
  
  # Case B: user provides logpred (nonstandard). Adapt it to log_pointpred output.
  if (!is.null(logpred)) {
    user_fn <- resolve_fun(logpred)
    
    lp_fn <- function(fit, data = NULL, config = list(), ...) {
      out <- user_fn(fit = fit, data = data, config = config, ...)
      
      # if already standardized
      if (is.list(out) && !is.null(out$lsf_hat) && !is.null(out$lpmf_hat) && !is.null(out$y_type)) {
        return(out)
      }
      
      # accept log_sf/log_pmf/y_type (or log_cdf)
      if (!is.list(out) || is.null(out$y_type)) {
        stop("logpred adapter: output must be a list and contain y_type.", call. = FALSE)
      }
      
      log_sf <- out$log_sf
      if (is.null(log_sf) && !is.null(out$log_cdf)) {
        cdf <- exp(out$log_cdf)
        log_sf <- log1p(-cdf)
      }
      
      if (is.null(log_sf)) stop("logpred adapter: missing log_sf (or log_cdf).", call. = FALSE)
      
      list(
        lsf_hat  = log_sf,
        lpmf_hat = out$log_pmf %||% NA_real_,
        y_type   = out$y_type
      )
    }
    
    return(log_summary_pred(fit, data = data, config = config, log_pointpred_fn = lp_fn, ...))
  }
  
  # Case C: default path (use your package dispatcher log_pointpred)
  log_summary_pred(fit, data = data, config = config, log_pointpred_fn = NULL, ...)
}

#' @export 
Zresidual <- function(fit,
                      data = NULL,
                      log_pointpred = NULL,
                      logpred = NULL,
                      config = list(),
                      randomized = TRUE,
                      nrep = 30,
                      eps = 1e-12,
                      seed = NULL,
                      ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  .log_add_exp2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
  
  .clip_logprob <- function(logp, eps = 1e-12) {
    lo <- log(eps)
    hi <- log1p(-eps)
    pmin(pmax(logp, lo), hi)
  }
  
  # 1) get summarized predictions (vectors)
  pre <- get_logpre(
    fit = fit, data = data,
    log_pointpred = log_pointpred,
    logpred = logpred,
    config = config,
    ...
  )
  
  logS <- as.numeric(pre$log_surv_hat)
  logp <- as.numeric(pre$log_pmf_hat)
  y_type <- pre$y_type
  
  n <- length(logS)
  if (n < 1L) stop("Zresidual: empty prediction.", call. = FALSE)
  
  # 2) determine whether pmf is available (discrete) or not (continuous)
  pmf_available <- !(is.null(logp) || all(is.na(logp)))
  
  # 3) Monte Carlo replications for randomized residuals
  if (!randomized) nrep <- 1L
  nrep <- as.integer(nrep)
  if (nrep < 1L) stop("Zresidual: nrep must be >= 1.", call. = FALSE)
  
  logS <- .clip_logprob(logS, eps = eps)
  
  log_rsp <- matrix(NA_real_, nrow = n, ncol = nrep)
  
  for (r in seq_len(nrep)) {
    u <- if (randomized) runif(n) else rep(0.5, n)
    logu <- log(pmax(u, .Machine$double.xmin))
    
    if (pmf_available) {
      # discrete: rsp = S + u * p  (log-space)
      log_rsp[, r] <- .clip_logprob(.log_add_exp2(logS, logp + logu), eps = eps)
    } else {
      # continuous survival: event => rsp=S; censored => rsp=u*S
      lr <- logS
      if (!is.null(y_type) && length(y_type) == n && all(stats::na.omit(unique(y_type)) %in% c(0, 1))) {
        cens_id <- which(y_type == 0L)
        if (length(cens_id)) lr[cens_id] <- logS[cens_id] + logu[cens_id]
      }
      log_rsp[, r] <- .clip_logprob(lr, eps = eps)
    }
  }
  
  rsp <- exp(log_rsp)
  z <- -qnorm(rsp)
  
  if (!is.matrix(z)) z <- matrix(z, nrow = n, ncol = nrep)
  colnames(z) <- paste0("rep", seq_len(nrep))
  
  class(z) <- c("zresid", class(z))
 # attr(z, "log_rsp") <- log_rsp
  attr(z, "rsp") <- rsp
  attr(z, "y_type") <- y_type
  z
}
