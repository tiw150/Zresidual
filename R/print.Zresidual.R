#' @export
print.zresid <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 2L) {
    cat(sprintf("<zresid> n=%d, reps=%d\n", d[1], d[2]))
    print(utils::head(unclass(x)), ...)
  } else {
    cat("<zresid>\n")
    print(unclass(x), ...)
  }
  invisible(x)
}

#' @export
print.cvzresid <- function(x, ...) {
  d <- dim(x)
  if (length(d) == 2L) {
    cat(sprintf("<cvzresid> n=%d, reps=%d\n", d[1], d[2]))
    print(utils::head(unclass(x)), ...)
  } else {
    cat("<cvzresid>\n")
    print(unclass(x), ...)
  }
  invisible(x)
}
