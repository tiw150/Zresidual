#' @export plot.dev.resid
#' @param Deviance.residual
#'
plot.dev.resid <- function(Deviance.residual,ylab = "Deviance Residual",
                           X = c("index", "lp", "covariate"),
                           main.title = "Deviance Residual Plot",
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
  if (missing(X)) X = "lp"

  id.infinity <- which (!is.finite(Deviance.residual))
  if (length(id.infinity) > 0L) {
    value.notfinite <- as.character.na(Deviance.residual[id.infinity])
    max.non.infinity <- max(abs(Deviance.residual[-id.infinity]))
    Deviance.residual[id.infinity] <-
      sign.na(Deviance.residual[id.infinity]) * (max.non.infinity + 0.1)
    message("Non-finite deviance residual exist! The model or the fitting process has a problem!")
  }
  censored <- attr(Deviance.residual, "censored.status")

  if (X == "index") {
    plot.default (
      Deviance.residual,
      ylab = ylab,
      ylim = c(min(Deviance.residual), max(Deviance.residual) + 1),
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
          Deviance.residual[which(is.outlier)],
          circles = rep(5, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          which(is.outlier),
          Deviance.residual[which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }

    }
  }
  if (X == "lp") {
    fitted.value <- attr(Deviance.residual, "linear.pred")
    plot(
      fitted.value,
      Deviance.residual,
      ylab = ylab,
      ylim = c(min(Deviance.residual), max(Deviance.residual) + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      main = main.title,
      xlab = "Linear Predictor"
    )
    lines(lowess(Deviance.residual ~ fitted.value),
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
          Deviance.residual[which(is.outlier)],
          circles = rep(0.03, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[which(is.outlier)],
          Deviance.residual[which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }
    }
  }

  if (X != "index" && X != "lp") {

    fitted.value <- attr(Deviance.residual, "covariates")

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
      Deviance.residual,
      ylab = ylab,
      ylim = c(min(Deviance.residual), max(Deviance.residual) + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      xlab = colnames(fitted.value)[i],
      main = main.title
    )
    lines(lowess(Deviance.residual ~ fitted.value[, i]),
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
          fitted.value[, i][which(is.outlier)],
          Deviance.residual[which(is.outlier)],
          circles = rep(5, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[,i][which(is.outlier)],
          Deviance.residual[which(is.outlier)],
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
