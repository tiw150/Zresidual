#' Plot martingale residuals for survival models
#'
#' Produce diagnostic plots for martingale residuals from survival models,
#' with the option to plot against the observation index, the linear predictor,
#' or a selected covariate. The function expects a vector (or one-column
#' matrix) of martingale residuals with attributes attached by the residual
#' computation functions in this package.
#'
#' The input \code{x} is typically obtained from
#' \code{residual.coxph()}, \code{residual.coxph.frailty()}, or
#' \code{residual.survreg()} with \code{residual.type = "martingale"}.
#'
#' The \code{x_axis_var} argument controls the x-axis:
#' \itemize{
#'  \item \code{"index"}: plot martingale residuals against observation index.
#'  \item \code{"lp"}: plot martingale residuals against the linear predictor
#'  (attribute \code{"linear.pred"}).
#'  \item \code{"covariate"}: prompt the user and print the available
#'  covariate names to the console.
#'  \item a character string matching one of the covariate names in
#'  \code{attr(x, "covariates")}: plot martingale residuals against that covariate.
#' }
#'
#' In the \code{"lp"} and covariate cases, a LOWESS smooth is added to the
#' plot to highlight systematic patterns in the residuals.
#'
#' @param x Numeric vector (or one-column matrix) of
#'  martingale residuals, typically returned by one of the residual
#'  functions in this package with \code{residual.type = "martingale"}.
#'  It must carry the attributes \code{"censored.status"},
#'  \code{"linear.pred"}, and \code{"covariates"} as described above.
#' @param ylab Character string for the y-axis label. Default is
#'  \code{"Martingale Residual"}.
#' @param x_axis_var Character string controlling the x-axis. Must be one of
#'  \code{"index"}, \code{"lp"}, \code{"covariate"}, or the name of a
#'  covariate contained in \code{attr(x, "covariates")}.
#'  The default is effectively \code{"lp"} if \code{x_axis_var} is not supplied.
#' @param main.title Character string for the main plot title. Default is
#'  \code{"Martingale Residual Plot"}.
#' @param outlier.return Logical; if \code{TRUE}, attempted outliers (as
#'  indicated by an external logical vector \code{is.outlier} in the calling
#'  environment) are highlighted in the plot and their indices are returned
#'  invisibly. If \code{FALSE} (default), no outlier indices are returned.
#'  Note that this function does not compute outliers internally: it assumes
#'  that a logical vector \code{is.outlier} of the same length as
#'  \code{x} is available if outlier highlighting is desired.
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @details
#' Non-finite martingale residuals are detected and truncated to lie slightly
#' beyond the largest finite residual, with a warning message printed to alert
#' the user that there may be problems with the model fit. Censored and
#' uncensored observations are distinguished by color and plotting symbol in
#' all display modes.
#'
#' @return
#' The function is primarily called for its side-effect of producing a plot.
#' If \code{outlier.return = TRUE}, it prints the indices of outlying points
#' to the console and invisibly returns a list with component
#' \code{outliers}, containing the indices where \code{is.outlier} is
#' \code{TRUE}. Otherwise, it returns \code{NULL} invisibly.
#'
#' @seealso
#' \code{residual.coxph}, \code{residual.coxph.frailty},
#' \code{residual.survreg}, \code{\link[survival]{Surv}},
#' \code{\link[survival]{survfit}}
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' data(lung)
#' fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#' r_m <- residual.coxph(fit, newdata = lung,
#'                       residual.type = "martingale")
#'
#' ## Basic plot vs. index
#' plot(r_m, x_axis_var = "index")
#'
#' ## Plot vs. linear predictor
#' plot(r_m, x_axis_var = "lp")
#'
#' ## Plot vs. a specific covariate, e.g. "age"
#' plot(r_m, x_axis_var = "age")
#' }
#' @method plot martg.resid
#' @export plot.martg.resid
plot.martg.resid <- function(x,ylab = "Martingale Residual",
                             x_axis_var = c("index", "lp", "covariate"),
                             main.title = "Martingale Residual Plot",
                             outlier.return = FALSE,
                             ...)
{
  Martingale.residual <- x
  X <- x_axis_var
  sign.na <- function(x)
  {
    sign.x <- sign(x)
    sign.x[is.na(x)] <- 1
    sign.x
  }

  as.character.na <- function(x)
  {
    label.x <- as.character(x)
    label.x[is.na(x)] <- "NA"
    label.x
  }
  if (missing(X)) X = "lp"

  id.infinity <- which (!is.finite(Martingale.residual))
  if (length(id.infinity) > 0L) {
    value.notfinite <- as.character.na(Martingale.residual[id.infinity])
    max.non.infinity <- max(abs(Martingale.residual[-id.infinity]))
    Martingale.residual[id.infinity] <-
      sign.na(Martingale.residual[id.infinity]) * (max.non.infinity + 0.1)
    message("Non-finite martingale residual exist! The model or the fitting process has a problem!")
  }
  ylim0<- max(Martingale.residual)
  censored <- attr(Martingale.residual, "censored.status")

  if (X == "index") {
    plot.default (
      Martingale.residual,
      ylab = ylab,
      ylim = c(min(Martingale.residual), max(Martingale.residual) + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      xlab = "Index",
      main = main.title
    )
    legend(
      x = "topleft",
      legend = c("Uncensored", "Censored"),
      col = c("darkolivegreen4", "blue"),
      pch = c(2, 3),
      cex = 0.8,
      xpd = TRUE,
      bty = "L",
      horiz = TRUE
    )

    if (isTRUE(outlier.return)) {
      if (identical(which(is.outlier), integer(0))) {
        return(invisible(NULL))
      } else {
        symbols(
          which(is.outlier),
          Martingale.residual[which(is.outlier)],
          circles = rep(5, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          which(is.outlier),
          Martingale.residual[which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }

    }
  }
  if (X == "lp") {
    fitted.value <- attr(Martingale.residual, "linear.pred")
    plot(
      fitted.value,
      Martingale.residual,
      ylab = ylab,
      ylim = c(min(Martingale.residual), max(Martingale.residual) + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      main = main.title,
      xlab = "Linear Predictor"
    )
    lines(lowess(Martingale.residual ~ fitted.value),
          col = "red",
          lwd = 3)
    legend(
      x = "topleft",
      legend = c("Uncensored", "Censored"),
      col = c("darkolivegreen4", "blue"),
      pch = c(2, 3),
      cex = 0.8,
      xpd = TRUE,
      bty = "L",
      horiz = TRUE
    )

    if (isTRUE(outlier.return)) {
      if (identical(which(is.outlier), integer(0))) {
        return(invisible(NULL))
      } else {
        symbols(
          fitted.value[which(is.outlier)],
          Martingale.residual[which(is.outlier)],
          circles = rep(0.03, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[which(is.outlier)],
          Martingale.residual[which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }
    }
  }

  if (X != "index" && X != "lp") {

    fitted.value <- attr(Martingale.residual, "covariates")

    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop("X must be the one of covariate name.") }


    plot(
      fitted.value[,i],
      Martingale.residual,
      ylab = ylab,
      ylim = c(min(Martingale.residual), max(Martingale.residual) + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      xlab = colnames(fitted.value)[i],
      main = main.title,
      ...
    )
    lines(lowess(Martingale.residual ~ fitted.value[, i]),
          col = "red",
          lwd = 3)
    legend(
      x = "topleft",
      legend = c("Uncensored", "Censored"),
      col = c("darkolivegreen4", "blue"),
      pch = c(2, 3),
      cex = 0.5,
      xpd = TRUE,
      bty = "L",
      horiz = TRUE
    )

    if (isTRUE(outlier.return)) {
      if (identical(which(is.outlier), integer(0))) {
        return(invisible(NULL))
      } else {
        symbols(
          fitted.value[, i][which(is.outlier)],
          Martingale.residual[which(is.outlier)],
          circles = rep(5, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[,i][which(is.outlier)],
          Martingale.residual[which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }
    }
  }


  if (outlier.return) {
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers = which(is.outlier)))
  }

}
