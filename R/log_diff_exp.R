.log1mexp <- function(x) {
  out <- rep(NaN, length(x))
  
  na_idx <- is.na(x)
  out[na_idx] <- NA_real_
  
  ok_idx <- !na_idx & x >= 0
  small_idx <- ok_idx & x <= log(2)
  large_idx <- ok_idx & x > log(2)
  
  out[small_idx] <- log(-expm1(-x[small_idx]))
  out[large_idx] <- log1p(-exp(-x[large_idx]))
  
  out
}

#' A function to compute the logarithm of the difference between the exponentials of two log values.
#'
#' @param la A log value.
#' @param lb A log value.
#' @return The logarithm of exp(la) - exp(lb) when la >= lb; otherwise NaN.
log_diff_exp <- function(la, lb) {
  d <- la - lb
  a <- la + 0 * d
  
  out <- rep(NaN, length(d))
  
  na_idx <- is.na(a) | is.na(d)
  out[na_idx] <- NA_real_
  
  ok_idx <- !na_idx & d >= 0
  out[ok_idx] <- a[ok_idx] + .log1mexp(d[ok_idx])
  
  out
}

log_sum_exp <- function(lx) {
  if (length(lx) == 0L) {
    return(-Inf)
  }
  
  mlx <- max(lx)
  if (is.infinite(mlx) && mlx < 0) {
    return(-Inf)
  }
  
  log(sum(exp(lx - mlx))) + mlx
}