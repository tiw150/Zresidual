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
  sf.test(Zresidual)$p.value
}

#' @export
sw.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  shapiro.test(Zresidual)$p.value
}

#' @export qqnorm.zresid
qqnorm.zresid <- function (Zresidual,main.title = "Normal Q-Q Plot",
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                           outlier.return=FALSE,...)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  is.outlier <- (abs(Zresidual) >3.5)
  sw.pv<-shapiro.test(Zresidual)$p.value

  if(max(abs(Zresidual))<6){
    qqnorm(Zresidual,main=main.title,xlab=xlab, ylab=ylab)
    qqline(Zresidual,col=1)
    abline(a=0,b=1,col=3)
    legend(x = "topleft",
           legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
  }
  if(max(abs(Zresidual))>=6){
    max.values<-which(Zresidual>=6)
    min.values<-which(Zresidual< -6)
    if(identical(min.values, integer(0))){
      gap=c(6,max(abs(Zresidual)))
      gapsize <- gap[2] - gap[1]
      ylim0<-c(min(Zresidual[-max.values]),7)
      xlim0<-c(min(qqnorm(Zresidual,plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual,plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual,plot.it = FALSE)[["x"]][-max.values],
           qqnorm(Zresidual,plot.it = FALSE)[["y"]][-max.values],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      axis(2, at=7 ,labels=round(max(Zresidual),1))
      qqline(Zresidual,col=1)
      abline(a=0,b=1,col=3)
      axis.break(2,6,style="gap")
      axis.break(2, 6, breakcol="snow", style="gap")
      axis.break(2, 6,breakcol="black", style="slash")
      axis.break(4, 6,breakcol="black", style="slash")
      points(qqnorm(Zresidual,plot.it = FALSE)[["x"]][max.values],
             pmax(qqnorm(Zresidual,plot.it = FALSE)[["y"]][max.values]-gapsize+1,6.5),col="red")
      legend(x = "bottomright",
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
    }

    if(identical(max.values, integer(0))){
      gap=c(min(Zresidual),6)
      gapsize2 <- gap[2] + gap[1]
      ylim0<-c(-7 ,max(Zresidual[-min.values]))
      xlim0<-c(min(qqnorm(Zresidual,plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual,plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual,plot.it = FALSE)[["x"]][-min.values],
           qqnorm(Zresidual,plot.it = FALSE)[["y"]][-min.values],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      qqline(Zresidual,col=1)
      abline(a=0,b=1,col=3)
      axis(2, at= -7 ,labels=round(min(Zresidual),1))
      axis.break(2,-6.1,style="gap")
      axis.break(2, -6.1, breakcol="snow", style="gap")
      axis.break(2, -6.1,breakcol="black", style="slash")
      axis.break(4, -6.1, breakcol="black", style="slash")
      points(qqnorm(Zresidual,plot.it = FALSE)[["x"]][min.values],
             pmin(qqnorm(Zresidual,plot.it = FALSE)[["y"]][min.values]-gapsize-1,-6.5),col="red")

      legend(x = "bottomright",
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))

    }

    if(!identical(max.values, integer(0)) && !identical(min.values, integer(0))){
      gap1=c(6,max(Zresidual))
      gapsize1 <- gap1[2] - gap1[1]
      gap2=c(min(Zresidual),6)
      gapsize2 <- gap2[2] + gap2[1]
      ylim0<-c(-7,7)
      xlim0<-c(min(qqnorm(Zresidual,plot.it = FALSE)[["x"]]),
               max(qqnorm(Zresidual,plot.it = FALSE)[["x"]]))
      plot(qqnorm(Zresidual,plot.it = FALSE)[["x"]][-c(max.values,min.values)],
           qqnorm(Zresidual,plot.it = FALSE)[["y"]][-c(max.values,min.values)],
           xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
      qqline(Zresidual,col=1)
      abline(a=0,b=1,col=3)
      axis(2, at=7 ,labels=round(max(Zresidual),1))
      axis.break(2,6,style="gap")
      axis.break(2, 6, breakcol="snow", style="gap")
      axis.break(2, 6,breakcol="black", style="slash")
      axis.break(4, 6,breakcol="black", style="slash")
      points(qqnorm(Zresidual,plot.it = FALSE)[["x"]][max.values],
             pmax(qqnorm(Zresidual,plot.it = FALSE)[["y"]][max.values]-gapsize1+1,6.5),col="red")

      axis(2, at= -7 ,labels=round(min(Zresidual),1))
      axis.break(2,-6.1,style="gap")
      axis.break(2, -6.1, breakcol="snow", style="gap")
      axis.break(2, -6.1,breakcol="black", style="slash")
      axis.break(4, -6.1, breakcol="black", style="slash")
      points(qqnorm(Zresidual,plot.it = FALSE)[["x"]][min.values],
             pmin(qqnorm(Zresidual,plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5),col="red")
      legend(x = "bottomright",
             legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
    }

  }

  if(isTRUE(outlier.return)){

    if(identical(which(is.outlier), integer(0))){
      return(invisible(NULL))
    } else {
      symbols((qqnorm(Zresidual,plot.it=FALSE)$x)[which(is.outlier)],
              (qqnorm(Zresidual,plot.it=FALSE)$y)[which(is.outlier)],
              circles=rep(0.1,length(which(is.outlier))),
              fg=rep('red',length(which(is.outlier))),
              add=T, inches=F)
      text((qqnorm(Zresidual,plot.it=FALSE)$x)[which(is.outlier)],
           (qqnorm(Zresidual,plot.it=FALSE)$y)[which(is.outlier)],
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
anov.test.zresid <- function (Zresidual,fitted.values, k.anova=10)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  test.nl.aov(Zresidual, fitted.values, k.anova)
}

#' @export bartlett.test.zresid
bartlett.test.zresid <- function (Zresidual,fitted.values, k.bl=10)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  test.var.bartl(Zresidual, fitted.values, k.bl)
}

#' @export boxplot.zresid
boxplot.zresid <- function(Zresidual,fitted.values, num.bin=10,
                           main.title="Z-residual Boxplot",
                           xlab=NULL,
                           outlier.return=FALSE,...)
{
  id.infinity <- which (!is.finite(Zresidual))
  if(length(id.infinity)>0L){
    value.notfinite <- as.character.na(Zresidual[id.infinity])
    max.non.infinity <- max(abs(Zresidual[-id.infinity]))
    Zresidual[id.infinity] <-
      sign.na(Zresidual[id.infinity])* (max.non.infinity + 0.1)
    message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
  }
  ylim0 <- max( qnorm(c(0.9999)), max(abs(Zresidual)))
  is.outlier <- (abs(Zresidual) >3)

  if(is.factor(fitted.values)){
    bin<-fitted.values
  }else{
    bin<-cut(fitted.values,num.bin)
  }
  plot(bin, Zresidual, ylab = "Z-Residual",
       ylim = c(-ylim0,ylim0+1),
       main=main.title,xlab=xlab
  )
  legend(x = "topleft",
         legend = c(paste0("Z-AOV p-value = ",
                         sprintf("%3.2f",anov.test.zresid(Zresidual,fitted.values, k.anova=num.bin))),
                    paste0("Z-BL p-value = ",
                           sprintf("%3.2f",bartlett.test.zresid(Zresidual,fitted.values, k.bl=num.bin)))),
         cex=0.6,horiz=TRUE)


  if(length(id.infinity)>0L){
    text(id.infinity+5, Zresidual[id.infinity],
         labels = value.notfinite,
         #adj = c(0,0),
         col=2)
  }
  hlines <- c(1.96,3); hlines2 <- -hlines
  abline(h= c(hlines,hlines2), lty=3, col="grey")
  if(outlier.return)
  {
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers=which(is.outlier)))
  }
}

#' @export plot.zresid
plot.zresid<- function(Zresidual,X=NULL,
                       main.title="Z-residual Scatterplot",xlab=NULL,
                       outlier.return=TRUE,...)
{
  id.infinity <- which (!is.finite(Zresidual))
  if(length(id.infinity)>0L){
    value.notfinite <- as.character.na(Zresidual[id.infinity])
    max.non.infinity <- max(abs(Zresidual[-id.infinity]))
    Zresidual[id.infinity] <-
      sign.na(Zresidual[id.infinity])* (max.non.infinity + 0.1)
    message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
  }
  ylim0 <- max( qnorm(c(0.9999)), max(abs(Zresidual)))
  is.outlier <- (abs(Zresidual) >3)
  censored<-attr(Zresidual, "censored.status")

  if(is.null(X)){
    plot.default (Zresidual, ylab = "Z-Residual",
                  ylim = c(-ylim0,ylim0+1),
                  col=c("blue","darkolivegreen4")[censored+1],
                  #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
                  pch=c(3,2)[censored+1],
                  main =main.title)
    legend(x = "topleft",
           legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
           pch=c(2,3),cex=0.5,xpd = TRUE,bty="L",horiz=TRUE)

    if(isTRUE(outlier.return)){

      if(identical(which(is.outlier), integer(0))){
        return(invisible(NULL))
      } else {
      symbols(which(is.outlier),
              Zresidual[which(is.outlier)],
              circles=rep(5,length(which(is.outlier))),
              fg=rep('red',length(which(is.outlier))),
              add=T, inches=F)
      text(which(is.outlier),
           Zresidual[which(is.outlier)],
           pos=1,label = which(is.outlier),
           cex = 0.8,col="red")
    }

  }
}
  if(!is.null(X)){
    plot(X, Zresidual, ylab = "Z-Residual",
         ylim = c(-ylim0,ylim0+1),
         col=c("blue","darkolivegreen4")[censored+1],
         #col = ifelse(is.outlier, "darkgoldenrod2", ifelse(censored,"darkolivegreen4","blue")),
         pch=c(3,2)[censored+1],
         main =main.title,xlab=xlab
    )
    lines(lowess(Zresidual ~ X),col = "red",lwd = 3)
    legend(x = "topleft",
           legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
           pch=c(2,3),cex=0.5,xpd = TRUE,bty="L",horiz=TRUE)

    if(isTRUE(outlier.return)){
      if(identical(which(is.outlier), integer(0))){
        return(invisible(NULL))
      } else {
        symbols(X[which(is.outlier)],
                Zresidual[which(is.outlier)],
                circles=rep(5,length(which(is.outlier))),
                fg=rep('red',length(which(is.outlier))),
                add=T, inches=F)
        text(X[which(is.outlier)],
             Zresidual[which(is.outlier)],
             pos=1,label = which(is.outlier),
             cex = 0.8,col="red")

      }

    }

  }


  if(length(id.infinity)>0L){
    text(id.infinity+5, Zresidual[id.infinity],
         labels = value.notfinite,
         #adj = c(0,0),
         col=2)
  }
  hlines <- c(1.96,3); hlines2 <- -hlines
  abline(h= c(hlines,hlines2), lty=3, col="grey")
  if(outlier.return){
    cat("Outlier Indices:", which(is.outlier), "\n")
    invisible(list(outliers=which(is.outlier)))
  }

}

