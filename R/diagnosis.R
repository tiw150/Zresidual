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

# sel.finite <- function(x) {
#     id.inf <- !is.finite(x);
#    # if(sum(id.inf)>1L) message("Infinity Exists:",which(id.inf));
#     x[!id.inf]
# }

#' @importFrom nortest sf.test
#' @export
sf.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  sf.pv<-rep(0,ncol(Zresidual))
  for(i in 1:ncol(Zresidual)){
    sf.pv[i]<-sf.test(Zresidual[,i])$p.value
  }
  sf.pv
}

#' @export
sw.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  sw.pv<-rep(0,ncol(Zresidual))
  for(i in 1:ncol(Zresidual)){
    sw.pv[i]<-shapiro.test(Zresidual[,i])$p.value
  }
  sw.pv
}

#' @export qqnorm.zresid
qqnorm.zresid <- function (Zresidual, index=1,
                           main.title = "Normal Q-Q Plot",
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                           outlier.return=FALSE,...)
{
  i<-index
  id.negtv.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] < 0)
  id.pos.inf <- which(is.infinite(Zresidual[,i]) & Zresidual[,i] > 0)
  Zresidual[,i][id.negtv.inf]<- -1e10
  Zresidual[,i][id.pos.inf]<- 1e10
  is.outlier <- (abs(Zresidual[,i]) >3.5)
  sw.pv<-shapiro.test(Zresidual[,i])$p.value

  if(max(abs(Zresidual[,i]))<6){
    qqnorm(Zresidual[,i],main=main.title,xlab=xlab, ylab=ylab)
    qqline(Zresidual[,i],col=1)
    abline(a=0,b=1,col=3)
    legend(x = "topleft",cex = 0.6,
           legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
  }
  if(max(abs(Zresidual[,i]))>=6){
    max.values<-which(Zresidual[,i]>=6)
    min.values<-which(Zresidual[,i]< -6)
    if(identical(min.values, integer(0))){
      gap=c(6,max(abs(Zresidual[,i])))
      gapsize <- gap[2] - gap[1]
      ylim0<-c(min(Zresidual[,i][-max.values]),7)
      xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-max.values],
           qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-max.values],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      axis(2, at=7 ,labels=round(max(Zresidual[,i]),1))
      qqline(Zresidual[,i],col=1)
      abline(a=0,b=1,col=3)
      axis.break(2,6,style="gap")
      axis.break(2, 6, breakcol="snow", style="gap")
      axis.break(2, 6,breakcol="black", style="slash")
      axis.break(4, 6,breakcol="black", style="slash")
      points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
             pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize+1,6.5),col="red")
      legend(x = "topleft",cex = 0.6,
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
    }

    if(identical(max.values, integer(0))){
      gap=c(min(Zresidual[,i]),6)
      gapsize2 <- gap[2] + gap[1]
      ylim0<-c(-7 ,max(Zresidual[,i][-min.values]))
      xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-min.values],
           qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-min.values],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      qqline(Zresidual[,i],col=1)
      abline(a=0,b=1,col=3)
      axis(2, at= -7 ,labels=round(min(Zresidual[,i]),1))
      axis.break(2,-6.1,style="gap")
      axis.break(2, -6.1, breakcol="snow", style="gap")
      axis.break(2, -6.1,breakcol="black", style="slash")
      axis.break(4, -6.1, breakcol="black", style="slash")
      points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][min.values],
             pmin(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5),col="red")

      legend(x = "topleft",cex = 0.6,
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))

    }

    if(!identical(max.values, integer(0)) && !identical(min.values, integer(0))){
      gap1=c(6,max(Zresidual[,i]))
      gapsize1 <- gap1[2] - gap1[1]
      gap2=c(min(Zresidual[,i]),6)
      gapsize2 <- gap2[2] + gap2[1]
      ylim0<-c(-7,7)
      xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-c(max.values,min.values)],
           qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-c(max.values,min.values)],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      qqline(Zresidual[,i],col=1)
      abline(a=0,b=1,col=3)
      axis(2, at=7 ,labels=round(max(Zresidual[,i]),1))
      axis.break(2,6,style="gap")
      axis.break(2, 6, breakcol="snow", style="gap")
      axis.break(2, 6,breakcol="black", style="slash")
      axis.break(4, 6,breakcol="black", style="slash")
      points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
             pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize1+1,6.5),col="red")

      axis(2, at= -7 ,labels=round(min(Zresidual[,i]),1))
      axis.break(2,-6.1,style="gap")
      axis.break(2, -6.1, breakcol="snow", style="gap")
      axis.break(2, -6.1,breakcol="black", style="slash")
      axis.break(4, -6.1, breakcol="black", style="slash")
      points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][min.values],
             pmin(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5),col="red")
      legend(x = "topleft",cex = 0.6,
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
    }

  }

  if(isTRUE(outlier.return)){

    if(identical(which(is.outlier), integer(0))){
      return(invisible(NULL))
    } else {
      symbols((qqnorm(Zresidual[,i],plot.it=FALSE)$x)[which(is.outlier)],
              (qqnorm(Zresidual[,i],plot.it=FALSE)$y)[which(is.outlier)],
              circles=rep(0.1,length(which(is.outlier))),
              fg=rep('red',length(which(is.outlier))),
              add=T, inches=F)
      text((qqnorm(Zresidual[,i],plot.it=FALSE)$x)[which(is.outlier)],
           (qqnorm(Zresidual[,i],plot.it=FALSE)$y)[which(is.outlier)],
           pos=1,label = which(is.outlier),
           cex = 0.8,col="red")
    }
  }
  if(outlier.return){
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers=which(is.outlier)))
  }
}



#' @export
gof.censore.zresid <- function (censored.Zresidual)
{
  id.negtv.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual < 0)
  id.pos.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual > 0)
  censored.Zresidual[id.negtv.inf]<- -1e10
  censored.Zresidual[id.pos.inf]<- 1e10
  censored.status<-attr(censored.Zresidual, "censored.status")
  censored.Zresidual<- as.vector(censored.Zresidual)
  gofTestCensored(censored.Zresidual,censored=censored.status, test = "sf",
                  censoring.side = "right",
                  distribution = "norm")$p.value
}

#' @export
anov.test.zresid <- function (Zresidual,X = c("lp", "covariate"), k.anova=10)
{
  if (missing(X))
    X = "lp"
  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")
    id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
    id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
    Zresidual[id.negtv.inf]<- -1e10
    Zresidual[id.pos.inf]<- 1e10
    aov.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      aov.pv[j]<- test.nl.aov(Zresidual[,j], fitted.value, k.anova)
    }
    aov.pv
    }
  if (X != "lp") {
    fitted.value <- attr(Zresidual, "covariates")
    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". "))}

    id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
    id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
    Zresidual[id.negtv.inf]<- -1e10
    Zresidual[id.pos.inf]<- 1e10
    aov.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      aov.pv[j]<- test.nl.aov(Zresidual[,j], fitted.value[,i], k.anova)
    }
    aov.pv

  }
}

#' @export bartlett.test.zresid
bartlett.test.zresid <- function (Zresidual, X = c("lp", "covariate"), k.bl=10)
{
  if (missing(X))
    X = "lp"

  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10

  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")
    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      bl.pv[j]<-test.var.bartl(Zresidual[,j], fitted.value, k.bl)
    }
    bl.pv
  }
  if (X != "lp") {
    fitted.value <- attr(Zresidual, "covariates")
    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". ")) }

    bl.pv<-rep(0,ncol(Zresidual))
    for(j in 1:ncol(Zresidual)){
      bl.pv[j]<-test.var.bartl(Zresidual[,j], fitted.value[,i], k.bl)
    }
    bl.pv
  }

}

#' @export boxplot.zresid
boxplot.zresid <- function(Zresidual,index=1,
                           X = c("lp", "covariate"),
                           num.bin = 10,
                           main.title="Z-residual Boxplot",
                           outlier.return = FALSE,
                           ...)
{
  j<-index

  if (missing(X))
    X = "lp"

  id.infinity <- which (!is.finite(Zresidual[,j]))
  if (length(id.infinity) > 0L) {
    value.notfinite <- as.character.na(Zresidual[,j][id.infinity])
    max.non.infinity <- max(abs(Zresidual[,j][-id.infinity]))
    Zresidual[,j][id.infinity] <-
      sign.na(Zresidual[,j][id.infinity]) * (max.non.infinity + 0.1)
    message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
  }
  ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[,j])))
  is.outlier <- (abs(Zresidual[,j]) > 3.5)


  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")

    if (is.factor(fitted.value)) {
      bin <- fitted.value
    } else{
      bin <- cut(fitted.value, num.bin)
    }
    plot(
      bin,
      Zresidual[,j],
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
      main = main.title,
      xlab = "Linear Predictor"
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("Z-AOV p-value = ",
               sprintf(
                 "%3.2f",
                 anov.test.zresid(Zresidual[,j], fitted.value, k.anova = num.bin)
               )),
        paste0("Z-BL p-value = ",
               sprintf(
                 "%3.2f",
                 bartlett.test.zresid(Zresidual[,j], fitted.value, k.bl = num.bin)
               ))
      ),
      cex = 0.6,
      horiz = TRUE
    )
  }

  if (X != "lp") {
    fitted.value <- attr(Zresidual, "covariates")
    if(X == "covariate"){
      i<-1
      cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
          variable.names(fitted.value))

    } else if(X %in% variable.names(fitted.value)){
      cov.name<-variable.names(fitted.value)
      i<- which(cov.name==X)
    } else{stop( paste0("X must be the one of covariate name: ", variable.names(fitted.value),". "))}

    if (is.factor(fitted.value[,i])) {
      bin <- fitted.value[,i]
    } else{
      bin <- cut(fitted.value[,i], num.bin)
    }
    plot(
      bin,
      Zresidual[,j],
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
      main = main.title,
      xlab = colnames(fitted.value)[i]
    )
    legend(
      x = "topleft",
      legend = c(
        paste0("Z-AOV p-value = ",
               sprintf(
                 "%3.2f",
                 anov.test.zresid(Zresidual[,j], fitted.value[,i], k.anova = num.bin)
               )),
        paste0("Z-BL p-value = ",
               sprintf(
                 "%3.2f",
                 bartlett.test.zresid(Zresidual[,j], fitted.value[,i], k.bl = num.bin)
               ))
      ),
      cex = 0.6,
      horiz = TRUE
    )
  }

  if (outlier.return)
  {
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers = which(is.outlier)))
  }
}


#' @export plot.zresid
plot.zresid <- function(Zresidual,index=1,
                        X = c("index", "lp", "covariate"),
                        main.title = "Z-residual Scatterplot",
                        outlier.return = FALSE,
                        ...)
{
  j<-index

  if (missing(X))
    X = "lp"
  id.infinity <- which (!is.finite(Zresidual[,j]))
  if (length(id.infinity) > 0L) {
    value.notfinite <- as.character.na(Zresidual[,j][id.infinity])
    max.non.infinity <- max(abs(Zresidual[,j][-id.infinity]))
    Zresidual[,j][id.infinity] <-
      sign.na(Zresidual[,j][id.infinity]) * (max.non.infinity + 0.1)
    message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
  }
  ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[,j])))
  is.outlier <- (abs(Zresidual[,j]) > 3.5)
  censored <- attr(Zresidual, "censored.status")

  if (X == "index") {
    plot.default (
      Zresidual[,j],
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
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
  }
  if (X == "lp") {
    fitted.value <- attr(Zresidual, "linear.pred")
    plot(
      fitted.value,
      Zresidual[,j],
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      main = main.title,
      xlab = "Linear Predictor"
    )
    lines(lowess(Zresidual[,j] ~ fitted.value),
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
          fitted.value[which(is.outlier)],
          Zresidual[,j][which(is.outlier)],
          circles = rep(0.03, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[which(is.outlier)],
          Zresidual[,j][which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }
    }
  }

  if (X != "index" && X != "lp") {

    fitted.value <- attr(Zresidual, "covariates")

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
      Zresidual[,j],
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
      col = c("blue", "darkolivegreen4")[censored + 1],
      #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
      pch = c(3, 2)[censored + 1],
      xlab = colnames(fitted.value)[i],
      main = main.title
    )
    lines(lowess(Zresidual[,j] ~ fitted.value[, i]),
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
          Zresidual[,j][which(is.outlier)],
          circles = rep(5, length(which(is.outlier))),
          fg = rep('red', length(which(is.outlier))),
          add = T,
          inches = F
        )
        text(
          fitted.value[,i][which(is.outlier)],
          Zresidual[,j][which(is.outlier)],
          pos = 1,
          label = which(is.outlier),
          cex = 0.8,
          col = "red"
        )
      }
    }
  }

  if (length(id.infinity) > 0L) {
    text(id.infinity + 5,
         Zresidual[,j][id.infinity],
         labels = value.notfinite,
         #adj = c(0,0),
         col = 2)
  }
  hlines <- c(1.96, 3)
  hlines2 <- -hlines
  abline(h = c(hlines, hlines2),
         lty = 3,
         col = "grey")
  if (outlier.return) {
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers = which(is.outlier)))
  }

}
