#' Bartlett test p-values for each column of a Z-residual matrix
#'
#' Plan A (Zresidual + metadata separated):
#' - Provide covariates / linear predictor via `zcov` returned by Zcov().
#' - Backward compatible with legacy attributes on Zresidual.
#'
#' @param x A numeric matrix of Z-residuals (n x p). Class "zresid" is allowed.
#' @param X X-axis specification. One of:
#'   - "index", "lp", "covariate"
#'   - a covariate name in `zcov$covariates`
#'   - a user-supplied vector of length nrow(x)
#' @param k.bl Integer. Number of bins for numeric covariates (default 10).
#' @param zcov Optional metadata returned by Zcov(). Should contain:
#'   - `linear_pred` (or `linear.pred`)
#'   - `covariates` (data.frame)
#' @param ... Reserved for forward compatibility.
#'
#' @return A numeric vector of Bartlett-test p-values, length = ncol(x).
#' @method bartlett.test zresid
#' @export bartlett.test.zresid
bartlett.test.zresid <- function(x,
                                 X = c("lp", "covariate"),
                                 k.bl = 10,
                                 zcov = NULL,
                                 ...) {
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  .as_matrix <- function(z) {
    if (is.null(dim(z))) matrix(z, ncol = 1)
    else z
  }
  
  .get_meta <- function(z, zcov) {
    if (!is.null(zcov)) return(zcov)
    
    zcov_attr <- attr(z, "zcov")
    if (!is.null(zcov_attr)) return(zcov_attr)
    
    # legacy fallback
    list(
      linear_pred = attr(z, "linear_pred") %||% attr(z, "linear.pred"),
      covariates  = attr(z, "covariates")
    )
  }
  
  .coerce_user_x <- function(xv) {
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
    
    if (missing(X)) X <- "lp"
    if (is.character(X) && length(X) > 1L) X <- X[1L]
    
    # user-supplied vector
    if (!is.null(X) && length(X) == n && !(is.character(X) && length(X) == 1L && X %in% keywords)) {
      return(.coerce_user_x(X))
    }
    
    if (!is.character(X) || length(X) != 1L) {
      stop("X must be: (1) a length-n vector, or (2) 'index'/'lp'/'covariate', or (3) a covariate name in zcov$covariates.",
           call. = FALSE)
    }
    
    if (X == "index") return(seq_len(n))
    
    if (X == "lp") {
      lp <- meta$linear_pred %||% meta$linear.pred %||% attr(z, "linear_pred") %||% attr(z, "linear.pred")
      if (is.null(lp)) stop("X='lp' requires `zcov$linear_pred` (Plan A) or legacy attr(Zresidual,'linear.pred').", call. = FALSE)
      if (length(lp) != n) stop("Length mismatch: length(linear_pred) != nrow(Zresidual). Make sure Zresidual and Zcov align.",
                                call. = FALSE)
      return(lp)
    }
    
    covs <- meta$covariates %||% attr(z, "covariates")
    if (is.null(covs)) stop("Covariates not found: provide `zcov$covariates` (Plan A) or legacy attr(Zresidual,'covariates').",
                            call. = FALSE)
    if (!is.data.frame(covs)) covs <- as.data.frame(covs)
    if (NROW(covs) != n) stop("Length mismatch: nrow(zcov$covariates) != nrow(Zresidual). Make sure Zresidual and Zcov align.",
                              call. = FALSE)
    
    cn <- colnames(covs) %||% names(covs)
    
    if (X == "covariate") {
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
  Zresidual <- .as_matrix(x)
  n <- NROW(Zresidual)
  p <- NCOL(Zresidual)
  
  meta <- .get_meta(Zresidual, zcov)
  xv <- .resolve_xv(Zresidual, X, meta)
  
  # Replace +/-Inf in residuals (matrix-wise)
  id_neginf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id_posinf <- which(is.infinite(Zresidual) & Zresidual > 0)
  if (length(id_neginf)) Zresidual[id_neginf] <- -1e10
  if (length(id_posinf)) Zresidual[id_posinf] <-  1e10
  
  bl.pv <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    zj <- Zresidual[, j]
    
    ok <- !is.na(zj) & !is.na(xv)
    if (sum(ok) < 3L) {
      bl.pv[j] <- NA_real_
      next
    }
    
    zj_ok <- zj[ok]
    xv_ok <- xv[ok]
    if (is.factor(xv_ok)) xv_ok <- droplevels(xv_ok)
    
    bl.pv[j] <- tryCatch(
      test.var.bartl(zj_ok, xv_ok, k.bl = k.bl),
      error = function(e) NA_real_
    )
  }
  
  bl.pv
}


#' Bartlett test core (internal): variance homogeneity across bins/levels
#'
#' @param Zresidual numeric vector
#' @param fitted.value numeric/factor vector (same length as Zresidual)
#' @param k.bl number of bins for numeric covariates
#' @return p-value (numeric) or NA_real_
#' @keywords internal
test.var.bartl <- function(Zresidual, fitted.value, k.bl = 10) {
  
  # pairwise complete
  ok <- !is.na(Zresidual) & !is.na(fitted.value)
  z  <- Zresidual[ok]
  x  <- fitted.value[ok]
  
  if (length(z) < 3L) return(NA_real_)
  
  # categorical/discrete
  if (is.factor(x) || length(unique(x)) <= k.bl) {
    f <- if (is.factor(x)) droplevels(x) else factor(x)
    
    # need >=2 non-empty groups
    tab <- table(f)
    f <- droplevels(f[tab[f] > 2L])
    z2 <- z[tab[factor(f, levels = levels(f))] > 2L]  # safe fallback; but easiest: recompute with keep
    keep <- tab[f] > 2L
    z2 <- z[keep]
    f2 <- droplevels(f[keep])
    
    if (nlevels(f2) < 2L) return(NA_real_)
    return(stats::bartlett.test(z2, f2)$p.value)
  }
  
  # numeric -> binning
  x_num <- as.numeric(x)
  
  .bin_equal_width <- function(v, k) droplevels(cut(v, k))
  .bin_quantile <- function(v, k) {
    br <- unique(stats::quantile(v, probs = seq(0, 1, length.out = k + 1), na.rm = TRUE, type = 7))
    if (length(br) < 3L) return(NULL) # <2 bins
    droplevels(cut(v, breaks = br, include.lowest = TRUE))
  }
  
  make_bins_and_filter <- function(v, k) {
    b <- .bin_equal_width(v, k)
    tab <- table(b)
    keep <- tab[b] > 2L
    list(bin = droplevels(b[keep]), z = z[keep])
  }
  
  # try equal-width bins
  tmp <- make_bins_and_filter(x_num, k.bl)
  bin <- tmp$bin
  zf  <- tmp$z
  
  # if insufficient bins, try log + quantile bins
  if (nlevels(bin) < 2L) {
    v2 <- log(pmax(x_num, .Machine$double.eps))
    bq <- .bin_quantile(v2, k.bl)
    if (!is.null(bq)) {
      tab <- table(bq)
      keep <- tab[bq] > 2L
      bin <- droplevels(bq[keep])
      zf  <- z[keep]
    } else {
      return(NA_real_)
    }
  }
  
  if (nlevels(bin) < 2L) return(NA_real_)
  stats::bartlett.test(zf, bin)$p.value
}
