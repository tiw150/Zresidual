#' A function to calculate the min p-value
#' @export
#' @param pv
#'
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
