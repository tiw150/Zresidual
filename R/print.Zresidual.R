#' Zresidual: Diagnostics with Z-residuals
#'
#'
#' @importFrom graphics boxplot plot
#' @keywords internal
"_PACKAGE"

#' @export
print.zresid <- function(x, ...){
  print(as.vector(x))
}

#' @export
print.cvzresid <- function(x, ...){
  print(as.vector(x))
}
