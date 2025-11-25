#' Cox–Snell residual plot for survival models
#'
#' Produce a Cox–Snell residual diagnostic plot for survival models, based on
#' the cumulative hazard of the Cox–Snell residuals. Under a correctly
#' specified model, the Cox–Snell residuals should follow an exponential
#' distribution with mean 1, so the nonparametric estimate of the cumulative
#' hazard should lie close to the 45-degree line.
#'
#' The input \code{cs.residual} is typically obtained from the residual
#' functions in this package (e.g., \code{residual.coxph()},
#' \code{residual.coxph.frailty()}, or \code{residual.survreg()}) with
#' \code{residual.type = "Cox-Snell"}.
#'
#' @importFrom survival Surv survfit
#'
#' @param cs.residual Numeric vector (or one-column matrix) of Cox–Snell
#'   residuals, typically returned by one of the residual functions in this
#'   package with \code{residual.type = "Cox-Snell"}. It must carry the
#'   attribute \code{"censored.status"}, giving the event indicator
#'   (1 = event, 0 = censored).
#' @param ylab Character string for the y-axis label. Default is
#'   \code{"Cumulative Hazard Function"}.
#' @param main.title Character string for the main plot title. Default is
#'   \code{"Cox-Snell Residuals Scatterplot"}.
#' @param outlier.return Logical; if \code{TRUE}, potential outliers are
#'   identified using a simple cutoff \code{|cs.residual| > 3.5}. Their
#'   indices are printed to the console and returned invisibly. If
#'   \code{FALSE} (default), no outlier indices are returned. Note that the
#'   current implementation attempts to highlight outliers using additional
#'   plotting calls and assumes access to objects named \code{Zresidual} and
#'   \code{j} in the calling environment; users may wish to adapt this part
#'   of the code for their own workflows.
#' @param ... Additional arguments passed to \code{\link[survival]{plot.survfit}}
#'   or to the underlying base graphics functions.
#'
#' @details
#' Non-finite Cox–Snell residuals are detected and truncated to lie slightly
#' beyond the largest finite residual, with a message printed to alert the
#' user that there may be problems with the model fit. The cumulative hazard
#' is drawn using \code{fun = "cumhaz"} and compared visually to the
#' exponential(1) reference line \eqn{H(t) = t}.
#'
#' @return
#' The function is primarily called for its side-effect of producing a plot.
#' If \code{outlier.return = TRUE}, it prints the indices of points flagged
#' as outliers (\code{|cs.residual| > 3.5}) and invisibly returns a list with
#' component \code{outliers}, containing these indices. Otherwise, it returns
#' \code{NULL} invisibly.
#'
#' @seealso
#' \code{residual.coxph}, \code{residual.coxph.frailty},
#' \code{residual.survreg}, \code{\link[survival]{Surv}},
#' \code{\link[survival]{survfit}}
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   data(lung)
#'   fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   cs_resid <- residual.coxph(fit, newdata = lung,
#'                              residual.type = "Cox-Snell")
#'
#'   ## Cox–Snell residual plot
#'   plot.cs.residual(cs_resid)
#'
#'   ## Return indices of large residuals
#'   out <- plot.cs.residual(cs_resid, outlier.return = TRUE)
#'   out$outliers
#' }
#' @export plot.cs.residual
#'
plot.cs.residual <- function(cs.residual,ylab = "Cumulative Hazard Function",
                             main.title = "Cox-Snell Residuals Scatterplot",
                             outlier.return = FALSE,
                             ...)
{
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
  censored <- attr(cs.residual, "censored.status")
  km <- survfit(Surv(cs.residual, censored)~1,type='fleming')
  id <- order(cs.residual)

  id.infinity <- which (!is.finite(cs.residual))
  if (length(id.infinity) > 0L) {
    value.notfinite <- as.character.na(cs.residual[id.infinity])
    max.non.infinity <- max(abs(cs.residual[-id.infinity]))
    cs.residual[id.infinity] <-
      sign.na(cs.residual[id.infinity]) * (max.non.infinity + 0.1)
    message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
  }

  is.outlier <- (abs(cs.residual) > 3.5)

  plot(
    km,
    fun  = "cumhaz",
    xlab = "Cox-Snell Residuals",
    ylab = ylab,
    main = main.title,
    ...
  )
  abline(0, 1, col="red", lty=2)
  points(km$time, -log(km$surv),
         col=c("blue","darkolivegreen4")[censored[id]+1],
         pch=c(3,2)[censored[id]+1] )
  legend(x = "topleft",
         legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
         pch=c(2,3),cex=1,xpd = TRUE,bty="L")


  if (isTRUE(outlier.return)) {
    if (identical(which(is.outlier), integer(0))) {
      return(invisible(NULL))
    } else {
      symbols(
        which(is.outlier),
        Zresidual[,j][which(is.outlier)],
        circles = rep(5, length(which(is.outlier))),
        fg = rep('red', length(which(is.outlier))),
        add = T,
        inches = F
      )
      text(
        which(is.outlier),
        Zresidual[,j][which(is.outlier)],
        pos = 1,
        label = which(is.outlier),
        cex = 0.8,
        col = "red"
      )
    }

  }

  if (outlier.return) {
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers = which(is.outlier)))
  }

}



