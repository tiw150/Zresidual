#' A function to draw scatter plot of a Cox-snell Residual
#'
#' @param cs.residual A numeric vector of Cox-Snell residuals.
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

  plot(km, fun="cumhaz", xlab=("Cox-Snell Residuals"),
       ylab=ylab,main=main.title)
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



