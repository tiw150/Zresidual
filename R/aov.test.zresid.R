#' ANOVA test diagnostic test for Z-residuals
#'
#' @description
#' Computes an ANOVA-style diagnostic p-value for each column of a Z-residual
#' matrix against an index, linear predictor, covariate, or user-supplied
#' grouping variable.
#'
#' @param Zresidual A numeric vector or matrix of Z-residuals.
#' @param X X-axis specification or grouping variable.
#' @param k.anova Maximum number of bins for numeric \code{X}.
#' @param zcov Optional metadata returned by \code{\link{Zcov}}.
#' @param ... Reserved for forward compatibility.
#'
#' @return A numeric vector of p-values.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.2 * x))
#'   t_cens  <- rexp(n, rate = 0.5)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   z <- Zresidual(fit,  data = dat, nrep = 2, seed = 1)
#'   info <- Zcov(fit, data = dat)
#'   aov.test.zresid(z, X = "lp", zcov = info)
#' }
#'
#' @export
aov.test.zresid <- function(Zresidual,
                            X = c("lp", "covariate"),
                            k.anova = 10,
                            zcov = NULL,
                            ...) {
  
  # ---- helpers ----
  .as_matrix <- function(x) {
    if (is.null(dim(x))) matrix(x, ncol = 1)
    else x
  }
  
  .get_meta <- function(z, zcov) {
    if (!is.null(zcov)) return(zcov)
    
    # allow attaching meta as an attribute without forcing it (still "separated" in API)
    zcov_attr <- attr(z, "zcov")
    if (!is.null(zcov_attr)) return(zcov_attr)
    
    # fallback to legacy attributes stored on Zresidual itself
    list(
      linear_pred = attr(z, "linear_pred") %||% attr(z, "linear.pred"),
      covariates  = attr(z, "covariates")
    )
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  .coerce_user_x <- function(xv) {
    # keep numeric as numeric; logical -> integer; character -> numeric if fully numeric else factor
    if (is.logical(xv)) return(as.integer(xv))
    if (is.factor(xv))  return(droplevels(xv))
    if (is.character(xv)) {
      xv_num <- suppressWarnings(as.numeric(xv))
      if (all(!is.na(xv_num))) return(xv_num)
      return(droplevels(factor(xv)))
    }
    xv
  }
  
  .resolve_xv <- function(z, X, meta) {
    n <- NROW(z)
    
    keywords <- c("index", "lp", "covariate")
    
    # default behavior when user does not specify X
    if (missing(X)) X <- "lp"
    
    # If X is the default vector c("lp","covariate") or similar, take first element.
    if (is.character(X) && length(X) > 1L) X <- X[1L]
    
    # 1) user-supplied vector
    if (!is.null(X) && length(X) == n && !(is.character(X) && length(X) == 1L && X %in% keywords)) {
      return(.coerce_user_x(X))
    }
    
    # 2) keyword or covariate name
    if (!is.character(X) || length(X) != 1L) {
      stop("X must be: (1) a length-n vector, or (2) 'index'/'lp'/'covariate', or (3) a covariate name in zcov$covariates.",
           call. = FALSE)
    }
    
    if (X == "index") {
      return(seq_len(n))
    }
    
    if (X == "lp") {
      lp <- meta$linear_pred %||% meta$linear.pred %||% attr(z, "linear_pred") %||% attr(z, "linear.pred")
      if (is.null(lp)) stop("X='lp' requires `zcov$linear_pred` (Plan A) or legacy attr(Zresidual,'linear.pred').", call. = FALSE)
      if (length(lp) != n) {
        stop("Length mismatch: length(linear_pred) != nrow(Zresidual). Make sure Zresidual and Zcov were computed on the same (post-NA-drop) rows.",
             call. = FALSE)
      }
      return(lp)
    }
    
    # covariates
    covs <- meta$covariates %||% attr(z, "covariates")
    if (is.null(covs)) stop("Covariates not found: provide `zcov$covariates` (Plan A) or legacy attr(Zresidual,'covariates').", call. = FALSE)
    
    if (!is.data.frame(covs)) covs <- as.data.frame(covs)
    
    if (NROW(covs) != n) {
      stop("Length mismatch: nrow(zcov$covariates) != nrow(Zresidual). Make sure Zresidual and Zcov align.",
           call. = FALSE)
    }
    
    cn <- colnames(covs) %||% names(covs)
    
    if (X == "covariate") {
      # keep legacy behavior: print names, then use the first one
      cat("To test against other covariates, set X to a covariate name. Available covariates:\n",
          paste(cn, collapse = ", "), "\n")
      if (length(cn) < 1L) stop("No covariates available in zcov$covariates.", call. = FALSE)
      return(covs[[cn[1L]]])
    }
    
    if (!(X %in% cn)) {
      stop(sprintf("Unknown covariate name '%s'. Available covariates: %s", X, paste(cn, collapse = ", ")),
           call. = FALSE)
    }
    
    covs[[X]]
  }
  
  # ---- main ----
  Zresidual <- .as_matrix(Zresidual)
  n <- NROW(Zresidual)
  p <- NCOL(Zresidual)
  
  meta <- .get_meta(Zresidual, zcov)
  xv <- .resolve_xv(Zresidual, X, meta)
  
  # Replace +/-Inf inside Zresidual with large finite values (matrix-wise, consistent with old code)
  id_neginf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id_posinf <- which(is.infinite(Zresidual) & Zresidual > 0)
  if (length(id_neginf)) Zresidual[id_neginf] <- -1e10
  if (length(id_posinf)) Zresidual[id_posinf] <-  1e10
  
  aov.pv <- rep(NA_real_, p)
  
  for (j in seq_len(p)) {
    zj <- Zresidual[, j]
    
    ok <- !is.na(zj) & !is.na(xv)
    if (sum(ok) < 3L) {
      aov.pv[j] <- NA_real_
      next
    }
    
    zj_ok <- zj[ok]
    xv_ok <- xv[ok]
    if (is.factor(xv_ok)) xv_ok <- droplevels(xv_ok)
    
    aov.pv[j] <- tryCatch(
      test.nl.aov(zj_ok, xv_ok, k.anova = k.anova),
      error = function(e) NA_real_
    )
  }
  
  aov.pv
}

#' ANOVA Test Core (internal)
#'
#' Performs an ANOVA test for Z-residuals against a covariate.
#' Numeric covariates are binned; sparse bins (<=2 obs) are removed.
#'
#' @param Zresidual Numeric vector.
#' @param fitted.value Numeric or factor covariate.
#' @param k.anova Max bins for numeric covariates.
#' @keywords internal
test.nl.aov <- function(Zresidual, fitted.value, k.anova = 10) {
  
  # treat categorical / discrete
  if (is.factor(fitted.value) || length(unique(fitted.value)) <= k.anova) {
    fv <- if (is.factor(fitted.value)) droplevels(fitted.value) else factor(fitted.value)
    fit <- stats::lm(Zresidual ~ fv)
    return(stats::anova(fit)[["Pr(>F)"]][1])
  }
  
  # numeric -> bin
  fv <- as.numeric(fitted.value)
  
  make_bins <- function(v) droplevels(cut(v, k.anova))
  
  bin <- make_bins(fv)
  tab <- table(bin)
  bad <- names(tab)[tab <= 2]
  
  if (length(bad)) {
    keep <- !(bin %in% bad)
    bin2 <- droplevels(bin[keep])
    z2   <- Zresidual[keep]
  } else {
    bin2 <- bin
    z2   <- Zresidual
  }
  
  # if too few bins remain, log-transform and retry
  if (nlevels(bin2) < 2L) {
    fv2 <- log(pmax(fv, .Machine$double.eps))
    message("ANOVA binning: too few non-sparse bins; applying log() transform and re-binning.")
    bin <- make_bins(fv2)
    tab <- table(bin)
    bad <- names(tab)[tab <= 2]
    
    if (length(bad)) {
      keep <- !(bin %in% bad)
      bin2 <- droplevels(bin[keep])
      z2   <- Zresidual[keep]
    } else {
      bin2 <- bin
      z2   <- Zresidual
    }
  }
  
  if (nlevels(bin2) < 2L) return(NA_real_)
  
  fit <- stats::lm(z2 ~ bin2)
  stats::anova(fit)[["Pr(>F)"]][1]
}
