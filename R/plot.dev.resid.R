#' Plot deviance residuals for survival models
#'
#' Draws a ggplot2 diagnostic plot against observation index, linear predictor,
#' or a selected covariate. The original LOWESS and external `is.outlier`
#' behavior are retained.
#'
#' @param x A deviance residual object with `censored.status`, `linear.pred`,
#'   and `covariates` attributes as required by the selected x-axis.
#' @param ylab Y-axis label.
#' @param x_axis_var `"index"`, `"lp"`, `"covariate"`, or a covariate name.
#' @param main.title Plot title.
#' @param outlier.return Whether to use the calling environment's logical
#'   `is.outlier` vector and return its indices.
#' @param point.args Named arguments for `ggplot2::geom_point()`.
#' @param smooth.args Named arguments for the LOWESS `geom_line()`.
#' @param theme A ggplot2 theme.
#' @param ... Additional plotting arguments.
#'
#' @return Invisibly returns `NULL`, or `list(outliers = ...)` when
#'   `outlier.return = TRUE`, matching the original method.
#'
#' @method plot dev.resid
#' @export
plot.dev.resid <- function(x,
                           ylab = "Deviance Residual",
                           x_axis_var = c("index", "lp", "covariate"),
                           main.title = "Deviance Residual Plot",
                           outlier.return = TRUE,
                           point.args = list(),
                           smooth.args = list(),
                           theme = ggplot2::theme_bw(),
                           ...) {
  X <- if (missing(x_axis_var)) "lp" else x_axis_var
  is_outlier <- NULL
  if (isTRUE(outlier.return) &&
      exists("is.outlier", envir = parent.frame(), inherits = TRUE)) {
    is_outlier <- get("is.outlier", envir = parent.frame(), inherits = TRUE)
  }

  plot <- .plot_survival_residual_ggplot(
    x = x,
    residual_name = "Deviance",
    ylab = ylab,
    x_axis_var = X,
    main_title = main.title,
    outlier_return = outlier.return,
    is_outlier = is_outlier,
    point_args = point.args,
    smooth_args = smooth.args,
    theme = theme,
    dots = list(...)
  )
  print(plot)

  if (isTRUE(outlier.return)) {
    indices <- attr(plot, "survival_residual_outliers")
    message("Outlier Indices: ", paste(indices, collapse = " "))
    return(invisible(list(outliers = indices)))
  }
  invisible(NULL)
}
