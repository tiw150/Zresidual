#' A function to compute the logarithm of the difference between the exponentials of two log values.
#'
#' @param la A log value.
#' @param lb A log value.
#' @import DPQ

log_diff_exp <- function (la, lb)
{
  ifelse(la >= lb, la + log1mexp(la-lb), NaN)
}

# compute the logarithm of the sum of exponential values.
log_sum_exp <- function (lx)
{
  mlx <- max (lx)
  log (sum (exp (lx - mlx))) + mlx
}
