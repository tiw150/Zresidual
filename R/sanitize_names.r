#' Sanitize column names (internal)
#'
#' Replace problematic characters in column names
#' to get something closer to brms' naming scheme.
#'
#' @param x Character vector of names.
#' @return Character vector with sanitized names.
#' @keywords internal
sanitize_names <- function(x) {
  out <- x
  out <- gsub("\\[|\\]|\\(|\\)", "", out)
  out <- gsub(" ", "_", out, fixed = TRUE)
  out
}
