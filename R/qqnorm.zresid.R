#' Normal Q-Q Plot for Z-Residuals with Outlier Detection and Normality Diagnostics
#'
#' Produces a normal Q-Q plot for Z-residuals, with optional
#' Shapiro–Wilk normality testing, automatic handling of infinite and extreme values,
#' axis breaks for very large residuals, and visual annotation of detected outliers.
#' This diagnostic is designed for checking the normality assumption of Z-residuals
#' obtained from Bayesian predictive model checks (posterior, LOOCV, ISCV, etc.).
#'
#' @param y A numeric matrix of Z-residuals where each column corresponds
#'   to an iteration or predictive draw. The function also uses the attribute
#'   \code{"type"} (optional model name) to construct default titles.
#' @param irep Integer or vector of integers indicating which column(s) of
#'   \code{Zresidual} to plot. Defaults to \code{1}.
#' @param diagnosis.test Character string indicating the normality test to perform.
#'   Currently only \code{"SW"} (Shapiro–Wilk test; using \code{sw.test.zresid()})
#'   is supported.
#' @param main.title Main title of the plot. If missing, it is automatically
#'   constructed using the \code{"type"} attribute of \code{Zresidual}.
#' @param xlab,ylab Axis labels for the Q-Q plot.
#' @param outlier.return Logical; if \code{TRUE}, the function prints and returns
#'   the indices of detected outliers.
#' @param outlier.value Numeric threshold used to classify an observation as
#'   an outlier based on \code{|Zresidual| > outlier.value}. Default is \code{3.5}.
#' @param outlier.set Optional named list of graphical parameters passed to
#'   \code{\link[graphics]{symbols}} and \code{\link[graphics]{text}} for customizing
#'   the outlier annotation.
#' @param my.mar Numeric vector giving margin sizes, passed internally to
#'   \code{par(mar = ...)}. Default is \code{c(5, 4, 4, 6) + 0.1}.
#' @param legend.settings Optional named list of parameters to override default
#'   legend appearance settings.
#' @param ... Additional graphical arguments passed to \code{\link[stats]{qqnorm}}
#'   and \code{\link[graphics]{plot}}.
#'
#' @details
#' This function extends the base R Q-Q plot to better handle typical behavior
#' of Z-residuals in Bayesian predictive checking:
#'
#' \itemize{
#'   \item Infinite values (\code{Inf}/\code{-Inf}) are replaced with large finite
#'         values and trigger a warning.
#'   \item Very large Z-residuals (\code{|Z| > 6}) are shown using axis breaks to
#'         avoid plot distortion.
#'   \item Outliers (\code{|Z| > outlier.value}) are highlighted and labeled.
#'   \item Column-wise Shapiro–Wilk tests assess normality.
#'   \item Legends summarize model type, selected Q-Q lines, and diagnostic results.
#' }
#'
#' This diagnostic is suitable for Z-residuals such as randomized quantile
#' residuals, posterior predictive Z-residuals, LOOCV/ISCV Z-residuals, and
#' residuals from hurdle or zero-inflated Bayesian models.
#'
#' @return
#' Invisibly returns a list with:
#' \describe{
#'   \item{outliers}{An integer vector containing the indices of detected outliers.}
#' }
#' If no outliers are detected, an empty integer vector is returned.
#'
#' A Q-Q plot is produced as a side effect.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics}, 5(3), 236–244.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., Dunson, D. B.,
#' Vehtari, A., & Rubin, D. B. (2013).
#' \emph{Bayesian Data Analysis}. CRC Press.
#'
#' @seealso
#'   \code{\link{boxplot.zresid}},
#'   \code{\link[stats]{qqnorm}},
#'   \code{\link[stats]{qqline}},
#'   \code{\link[graphics]{plot}}.
#'
#' @examples
#' library(Zresidual)
#' set.seed(1)
#' Z <- matrix(rnorm(200), ncol = 2)
#' attr(Z, "type") <- "Example Model"
#'
#' # Basic Q-Q plot
#' qqnorm.zresid(Z)
#'
#' # Use the second column with custom outlier threshold
#' qqnorm.zresid(Z, irep = 2, outlier.value = 2.5)
#'
#' # Modify legend settings
#' qqnorm.zresid(Z, legend.settings = list(cex = 0.8))
#'
#' @method qqnorm zresid
#' @export qqnorm.zresid
qqnorm.zresid <- function (y, irep=1, diagnosis.test = "SW",
                           #X.anova = c("lp", "covariate"),k.anova=10,
                           main.title = ifelse(is.null(attr(y, "type")),
                                               "Normal Q-Q Plot",
                                               paste("Normal Q-Q Plot -", attr(y, "type"))),
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                           outlier.return=TRUE, outlier.value = 3.5, outlier.set = list(),
                           my.mar = c(5, 4, 4, 6) + 0.1,legend.settings = list(), ...)
{
  #i<-index

  Zresidual <- y
  if(missing(diagnosis.test)) diagnosis.test <- "SW"
#  if(missing(X.anova) && diagnosis.test == "ANOVA") X.anova <- "lp"
  # test <- list("SW" = shapiro.test,
  #              "ANOVA" = aov.test.zresid)
  if(diagnosis.test == "SW") test <- sw.test.zresid(Zresidual)

  for (i in irep) {
    par(mar = my.mar)
    # Get the range of the QQ plot data to set the plot window
    qq_x <- qqnorm(Zresidual[,i], plot.it = FALSE)$x
    qq_y <- qqnorm(Zresidual[,i], plot.it = FALSE)$y
    xlim0_qq <- c(min(qq_x, na.rm = TRUE) - 0.5, max(qq_x, na.rm = TRUE) + 0.6)
    ylim0_qq <- range(qq_y, na.rm = TRUE) # Use the actual range of QQ values

    plot.window(xlim = xlim0_qq, ylim = ylim0_qq) # Explicitly set the plot window

   # type <- attr(Zresidual, "type")
    id.negtv.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] < 0)
    id.pos.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] > 0)
    Zresidual[,i][id.negtv.inf]<- -1e10
    Zresidual[,i][id.pos.inf]<- 1e10

    id.nan <- which(is.nan(Zresidual[,i]))
    id.infinity <- which (is.infinite(Zresidual[,i]))
    id.outlier <- which(abs(Zresidual[,i]) > outlier.value)

    if (length(id.infinity) > 0L) message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    if(length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")

    test.pv <- test[i]
    ### Diagnosis Test
    #default.test.args <- list(Zresidual, X.anova, k.anova)
    #test.pv <- do.call(test[[diagnosis.test]], default.test.args[!sapply(default.test.args, is.null)])[i]
    #shapiro.test(Zresidual[,i])$p.value

    ### Legend
    default.legend1 <- list(
      x = grconvertX(1.02, "npc", "user"),
      y = grconvertY(0.65, "npc", "user"),
      legend = c("qqline", "45\u00B0 Line") ,
      col = c("black", "green"),
      lty = c(1, 1),
      lwd = 1.5,
      cex = 0.6,
      bty = "n",
      xpd = NA,
      seg.len = 1.5,
      adj = 0
    )

    default.legend2 <- list(
      x = grconvertX(0.99, "npc", "user"),
      y = grconvertY(0.4, "npc", "user"),
      legend = c(expression(bold("P-value:")),
                 paste0("Z-", diagnosis.test, " = ", sprintf("%.2f", test.pv))),
      cex = 0.6,
      bty = "n",
      xpd = NA,
      adj = 0
    )
    legend.args1 <- modifyList(default.legend1, legend.settings)
    legend.args2 <- modifyList(default.legend2, legend.settings)

    ### outlier
    default.outlier <- list(
      pos = 4,
      labels = id.outlier,
      cex = 0.8,
      col = "red",
      add = T,
      inches = F,
      circles = rep((par("usr")[2]-par("usr")[1])*0.027, length(id.outlier)),
      fg = "red"
    )

    # Add user-provided `...` arguments, overwriting defaults if necessary
    outlier.args <- modifyList(default.outlier, outlier.set)

    # Extract necessary arguments seperately for text() and symbols()
    text.args <- outlier.args[names(outlier.args) %in% names(formals(text.default))]
    symbols.args <- outlier.args[names(outlier.args) %in% names(formals(symbols))]

    if(max(abs(Zresidual[,i]), na.rm = T)<6){
      xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)-0.5,
               max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)+0.6)
      x_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]
      y_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]]
      outlier_points <- (y_values > 3.5) | (y_values < -3.5)
      qqnorm(Zresidual[,i],main=main.title,xlab=xlab, ylab=ylab,xlim=xlim0, ...)
      qqline(Zresidual[,i],col=1)
      abline(a=0,b=1,col=3)
      do.call(legend, legend.args1)
      do.call(legend, legend.args2)

      if(isTRUE(outlier.return)){
        if(identical(id.outlier, integer(0))){
          #return(invisible(NULL))
          next
        } else {
          points(x_values[outlier_points], y_values[outlier_points], pch = 16)
          x.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier]
          y.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier]

          do.call(symbols, c(list(x=x.outlier,y=y.outlier), symbols.args))
          do.call(text, c(list(x=x.outlier,y=y.outlier), text.args))
        }
        }
    }

    if(max(abs(Zresidual[,i]), na.rm = T)>=6){
      max.values<-which(Zresidual[,i]>=6)
      min.values<-which(Zresidual[,i] < -6)

      if(identical(min.values, integer(0))){
        gap=c(6,max(abs(Zresidual[,i]), na.rm = T))
        gapsize1 <- gap[2] - gap[1]
        ylim0<-c(min(Zresidual[,i][-max.values], na.rm = T),7)
        xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)-0.5,
                 max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)+0.6)
        x_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]
        y_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]]
        outlier_points <- (y_values > 3.5 & y_values < 6) | (y_values < -3.5 & y_values > -6 )
        plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-max.values],
             qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-max.values],
             xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab, ...)
        axis(2, at=7 ,labels=round(max(Zresidual[,i], na.rm = T),1))
        qqline(Zresidual[,i],col=1)
        abline(a=0,b=1,col=3)
        points(x_values[outlier_points], y_values[outlier_points])
        axis.break(2,6,style="gap")
        axis.break(2, 6, breakcol="snow", style="gap")
        axis.break(2, 6,breakcol="black", style="slash")
        axis.break(4, 6,breakcol="black", style="slash")
        points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
               pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize1+1,6.5))
        do.call(legend, legend.args1)
        do.call(legend, legend.args2)

        if(isTRUE(outlier.return)){
          if(identical(id.outlier, integer(0))){
            #return(invisible(NULL))
            next
          } else {
            points(x_values[outlier_points], y_values[outlier_points], pch = 16)
            points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
                   pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize1+1,6.5),
                   col="black",pch=16)
            x.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier]
            y.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier]

            if(any(y.outlier > 6)) {
              idx_max <- which(y.outlier > 6)
              y.outlier[idx_max] <- pmax(y.outlier[idx_max] - gapsize1 + 1, 6.5)
            }

            do.call(symbols, c(list(x=x.outlier,y=y.outlier), symbols.args))
            do.call(text, c(list(x=x.outlier,y=y.outlier), text.args))
          }
        }

      }

      if(identical(max.values, integer(0))){
        gap=c(min(Zresidual[,i], na.rm = T),-6)
        gapsize2 <- gap[2] + gap[1]
        ylim0<-c(-7 ,max(Zresidual[,i][-min.values], na.rm = T))
        xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)-0.5,
                 max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T)+0.6)
        x_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]
        y_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]]
        outlier_points <- (y_values > 3.5 & y_values < 6) | (y_values < -3.5 & y_values > -6 )
        plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-min.values],
             qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-min.values],
             xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab, ...)
        qqline(Zresidual[,i],col=1)
        abline(a=0,b=1,col=3)
        points(x_values[outlier_points], y_values[outlier_points])
        axis(2, at= -7 ,labels=round(min(Zresidual[,i], na.rm = T),1))
        axis.break(2,-6.1,style="gap")
        axis.break(2, -6.1, breakcol="snow", style="gap")
        axis.break(2, -6.1,breakcol="black", style="slash")
        axis.break(4, -6.1, breakcol="black", style="slash")
        points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][min.values],
               pmin(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5))
        do.call(legend, legend.args1)
        do.call(legend, legend.args2)

        if(isTRUE(outlier.return)){
          if(identical(id.outlier, integer(0))){
            #return(invisible(NULL))
            next
          } else {
            points(x_values[outlier_points], y_values[outlier_points], pch = 16)
            points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][min.values],
                   pmin(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5),
                   col="black",pch=16)

            x.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier]
            y.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier]

            if (any(y.outlier < -6)) {
              idx_min <- which(y.outlier< -6)
              y.outlier[idx_min] <- pmin(y.outlier[idx_min] - gapsize2 - 1, -6.5)
            }

            do.call(symbols, c(list(x=x.outlier,y=y.outlier), symbols.args))
            do.call(text, c(list(x=x.outlier,y=y.outlier), text.args))
          }
        }

      }

      if(!identical(max.values, integer(0)) && !identical(min.values, integer(0))){
        gap1=c(6,max(Zresidual[,i],na.rm=T))
        gapsize3 <- gap1[2] - gap1[1]
        gap2=c(min(Zresidual[,i],na.rm=T),6)
        gapsize4 <- gap2[2] + gap2[1]
        x_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]
        y_values<-qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]]
        ylim0<-c(-7.5,7.5)
        xlim0<-c(min(x_values, na.rm = T)-0.5,
                 max(x_values, na.rm = T)+0.6)
        outlier_points <- (y_values > 3.5 & y_values < 6) | (y_values < -3.5 & y_values > -6 )
        plot(x_values[-c(max.values,min.values)],y_values[-c(max.values,min.values)],
             xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab,
             pch = 1, ...)
        points(x_values[outlier_points], y_values[outlier_points])
        qqline(Zresidual[,i],col=1)
        abline(a=0,b=1,col=3)
        axis(2, at=7 ,labels=round(max(Zresidual[,i], na.rm = T),1))
        axis.break(2,6,style="gap")
        axis.break(2, 6, breakcol="snow", style="gap")
        axis.break(2, 6,breakcol="black", style="slash")
        axis.break(4, 6,breakcol="black", style="slash")
        points(x_values[max.values],
               pmax(y_values[max.values]-gapsize3+1,6.5))

        axis(2, at= -7 ,labels=round(min(Zresidual[,i], na.rm = T),1))
        axis.break(2,-6.1,style="gap")
        axis.break(2, -6.1, breakcol="snow", style="gap")
        axis.break(2, -6.1,breakcol="black", style="slash")
        axis.break(4, -6.1, breakcol="black", style="slash")
        points(x_values[min.values],
               pmin(y_values[min.values]-gapsize4-1,-6.5))

        ### Legend
        default.legend1 <- list(
          x = grconvertX(1.02, "npc", "user"),
          y = grconvertY(0.65, "npc", "user"),
          legend = c("qqline", "Diagonal Line"),
          col = c("black", "green"),
          lty = c(1, 1),
          lwd = 1.5,
          cex = 0.6,
          bty = "n",
          xpd = NA,
          seg.len = 1.5,
          adj = 0
        )

        default.legend2 <- list(
          x = grconvertX(0.99, "npc", "user"),
          y = grconvertY(0.4, "npc", "user"),
          legend = c(expression(bold("P-value:")),
                     paste0("Z-", diagnosis.test, " = ", sprintf("%.2f", test.pv))),
          cex = 0.6,
          bty = "n",
          xpd = NA,
          adj = 0
        )
        legend.args1 <- modifyList(default.legend1, legend.settings)
        legend.args2 <- modifyList(default.legend2, legend.settings)

        do.call(legend, legend.args1)
        do.call(legend, legend.args2)

        if(isTRUE(outlier.return)){
          if(identical(id.outlier, integer(0))){
            #return(invisible(NULL))
            next
          } else {
            points(x_values[outlier_points], y_values[outlier_points], pch = 16)
            points(x_values[max.values],
                   pmax(y_values[max.values]-gapsize3+1,6.5),
                   col="black",pch=16)
            points(x_values[min.values],
                   pmin(y_values[min.values]-gapsize4-1,-6.5),
                   col="black",pch=16)

            x.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier]
            y.outlier<-(qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier]

            if(any(y.outlier > 6)) {
              idx_max <- which(y.outlier > 6)
              y.outlier[idx_max] <- pmax(y.outlier[idx_max] - gapsize3 + 1, 6.5)
            }

            if (any(y.outlier < -6)) {
              idx_min <- which(y.outlier< -6)
              y.outlier[idx_min] <- pmin(y.outlier[idx_min] - gapsize4 - 1, -6.5)
            }

            do.call(symbols, c(list(x=x.outlier,y=y.outlier), symbols.args))
            do.call(text, c(list(x=x.outlier,y=y.outlier), text.args))
          }
        }

      }

    }


    if(outlier.return){
      cat("Outlier Indices :", id.outlier, "\n")
      invisible(list(outliers=id.outlier))
    }

  }

}

