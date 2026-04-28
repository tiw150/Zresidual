# File: R/upper_bound_pvalue.R

# Internal cache for H(J,r)
.zres_H_cache <- new.env(parent = emptyenv())


#' Compute H_{J,r} for replicated p-values (internal)
#'
#' H_{J,r} = sup_{p in (0,1]} P{Bin(J,p) >= r} / p
#' @noRd
compute_HJr <- function(J,
                        r = ceiling(J / 2),
                        eps = 1e-12,
                        grid_size = 6000,
                        top_k = 30,
                        local_span_logit = 1.2) {
  if (!is.numeric(J) || length(J) != 1 || is.na(J) || J < 1) {
    stop("J must be a positive integer.")
  }
  J <- as.integer(J)
  
  if (!is.numeric(r) || length(r) != 1 || is.na(r) || r < 1 || r > J) {
    stop("r must be an integer in [1, J].")
  }
  r <- as.integer(r)
  
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0 || eps >= 1e-2) {
    stop("eps must be in (0, 1e-2) (small positive).")
  }
  
  logit <- function(x) log(x / (1 - x))
  inv_logit <- function(z) 1 / (1 + exp(-z))
  
  # Tail prob: P(Bin(J,p) >= r) = P(Bin(J,p) > r-1)
  tail_prob <- function(p) {
    stats::pbinom(q = r - 1, size = J, prob = p, lower.tail = FALSE)
  }
  
  obj <- function(p) {
    p <- pmax(p, eps)
    tail_prob(p) / p
  }
  
  # Coarse grid on logit-scale for robust candidate selection
  z_grid <- seq(logit(eps), logit(1 - eps), length.out = grid_size)
  p_grid <- inv_logit(z_grid)
  v_grid <- obj(p_grid)
  
  ord <- order(v_grid, decreasing = TRUE)
  k <- min(top_k, length(ord))
  cand <- p_grid[ord[seq_len(k)]]
  
  best_val <- -Inf
  best_p <- NA_real_
  
  # Local refine around top candidates
  for (p0 in cand) {
    z0 <- logit(p0)
    lo <- inv_logit(z0 - local_span_logit)
    hi <- inv_logit(z0 + local_span_logit)
    
    lo <- max(lo, eps)
    hi <- min(hi, 1 - eps)
    if (lo >= hi) next
    
    opt <- stats::optimize(f = function(p) -obj(p), lower = lo, upper = hi)
    val <- -opt$objective
    if (is.finite(val) && val > best_val) {
      best_val <- val
      best_p <- opt$minimum
    }
  }
  
  list(H = best_val, p_star = best_p, J = J, r = r)
}


#' Get cached H_{J,r} (internal)
#' @noRd
get_H <- function(J, r = ceiling(J / 2), ...) {
  if (!is.numeric(J) || length(J) != 1 || is.na(J) || J < 1) {
    stop("J must be a positive integer.")
  }
  J <- as.integer(J)
  
  if (!is.numeric(r) || length(r) != 1 || is.na(r) || r < 1 || r > J) {
    stop("r must be an integer in [1, J].")
  }
  r <- as.integer(r)
  
  key <- paste0("J=", J, "|r=", r)
  if (exists(key, envir = .zres_H_cache, inherits = FALSE)) {
    return(get(key, envir = .zres_H_cache, inherits = FALSE))
  }
  
  H <- compute_HJr(J = J, r = r, ...)$H
  assign(key, H, envir = .zres_H_cache)
  H
}

#' Distribution-free upper-bound p-value from replicated p-values
#'
#' @description
#' Combines replicated p-values into a single valid upper-bound p-value using
#' an order-statistic correction.
#'
#' @param p_rep Numeric vector of replicated p-values in \code{[0, 1]}.
#' @param r Integer order-statistic rank. The default is
#'   \code{ceiling(length(p_rep) / 2)}.
#' @param H Optional precomputed correction constant.
#' @param ... Additional arguments passed to the internal computation of
#'   \eqn{H_{J,r}} when \code{H = NULL}.
#'
#' @return A single numeric p-value in \code{[0, 1]}.
#'
#' @examples
#' p <- c(0.03, 0.08, 0.12, 0.05, 0.09)
#' upper_bound_pvalue(p)
#' upper_bound_pvalue(p, r = 2)
#'
#' @export
upper_bound_pvalue <- function(p_rep,
                               r = ceiling(length(p_rep) / 2),
                               H = NULL,
                               ...) {
  if (!is.numeric(p_rep) || length(p_rep) < 1) {
    stop("p_rep must be a non-empty numeric vector.")
  }
  if (anyNA(p_rep)) {
    stop("p_rep contains NA. Please remove or impute before calling.")
  }
  if (any(p_rep < 0 | p_rep > 1)) {
    stop("p_rep must be within [0,1].")
  }
  
  J <- length(p_rep)
  
  if (!is.numeric(r) || length(r) != 1 || is.na(r) || r < 1 || r > J) {
    stop("r must be an integer in [1, length(p_rep)].")
  }
  r <- as.integer(r)
  
  # Use the r-th ORDER STATISTIC (not median() when J is even)
  p_r <- sort(p_rep, decreasing = FALSE)[r]
  
  if (is.null(H)) {
    H <- get_H(J = J, r = r, ...)
  } else {
    if (!is.numeric(H) || length(H) != 1 || is.na(H) || H <= 0) {
      stop("H must be a single positive numeric value.")
    }
  }
  
  min(1, H * p_r)
}

