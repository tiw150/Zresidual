#' @export plot.martg.resid
#' @param Martingale.residual
#'
#'
plot.martg.resid <- function(Martingale.residual,ylab = "Martingale Residual",
                             X = c("index", "lp", "covariate"),
                             main.title = "Martingale Residual Plot",
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
      main = main.title
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
