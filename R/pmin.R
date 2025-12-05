#' A function to calculate the min p-value
#'
#' @param pv Numeric vector of p-values.
#'
#' @return A single numeric value: the minimum adjusted p-value.
#' @export
pvalue.min <- function (pv)
{
  pv <- pv[is.finite(pv)]
  n <- length (pv)
  if(n<=0) stop("There is no value in 'pv'")
  pv.sorted <- sort(pv)
  corrected.pv <- pmin(1, pv.sorted*n/1:n)
  pv.min <- min(corrected.pv)
  pv.min
}
