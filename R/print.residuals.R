#' Print a residual vector
#'
#' @description
#' Simple print method for residual-like vectors used in this package.
#'
#' @param x A residual object or numeric vector.
#' @param ... Further arguments passed to \code{print()}.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' x <- structure(rnorm(5), class = "resid")
#' print(x)
#'
#' @export
print.resid <- function(x, ...){
  print(as.vector(x))
}
