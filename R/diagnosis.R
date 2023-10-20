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

plot.zresid <- function(Zresidual,main.title="Z-Residual Scatterplot",
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
    plot.default (Zresidual, ylab = "Z-Residual",
          ylim = c(-ylim0,ylim0),
          col = c(1,2)[is.outlier + 1],
          main =main.title
    )
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


# sel.finite <- function(x) {
#     id.inf <- !is.finite(x);
#    # if(sum(id.inf)>1L) message("Infinity Exists:",which(id.inf));
#     x[!id.inf]
# }

sf.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  sf.test(Zresidual)$p.value
}

sw.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  shapiro.test(Zresidual)$p.value
}

qqnorm.zresid <- function (Zresidual,main.title = "Normal Q-Q Plot",
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
{
    id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
    id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
    Zresidual[id.negtv.inf]<- -1e10
    Zresidual[id.pos.inf]<- 1e10
    sw.pv<-shapiro.test(Zresidual)$p.value
    qqnorm(Zresidual,main=main.title,xlab=xlab, ylab=ylab)
    qqline(Zresidual,col=1)
    abline(a=0,b=1,col=2)
    legend(x = "topleft",
           legend = paste0("Z-SW p-value = ",sprintf("%3.2f",sw.pv)))
}


gof.censore.zresid <- function (censored.Zresidual,censored)
{
  id.negtv.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual < 0)
  id.pos.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual > 0)
  censored.Zresidual[id.negtv.inf]<- -1e10
  censored.Zresidual[id.pos.inf]<- 1e10

  gofTestCensored(censored.Zresidual,censored, test = "sf",
                  censoring.side = "right",
                  distribution = "norm")$p.value
}


anov.test.zresid <- function (Zresidual,fitted.values, k.anova=10)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  test.nl.aov(Zresidual, fitted.values, k.anova)
}

bartlett.test.zresid <- function (Zresidual,fitted.values, k.bl=10)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  test.var.bartl(Zresidual, fitted.values, k.bl)
}


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
       ylim = c(-ylim0,ylim0),
       main=main.title,xlab=xlab
  )
  legend(x = "topleft",
         legend = c(paste0("Z-AOV p-value = ",
                         sprintf("%3.2f",anov.test.zresid(Zresidual,fitted.values, k.anova=num.bin))),
                    paste0("Z-BL p-value = ",
                           sprintf("%3.2f",bartlett.test.zresid(Zresidual,fitted.values, k.bl=num.bin)))),
         cex=0.5,horiz=TRUE)



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


plot.zresid.fitted.value <- function(Zresidual,fitted.values,cenered,
                                     main.title="Z-residual Scatterplot",
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

  plot(fitted.values, Zresidual, ylab = "Z-Residual",
       ylim = c(-ylim0,ylim0),
       col=c("blue","darkolivegreen4")[cenered+1],
       pch=c(3,2)[cenered+1],
#       col = c(1,2)[is.outlier + 1],
       main =main.title,xlab=xlab
  )
  lines(lowess(Zresidual ~ fitted.values),col = "red",lwd = 3)
  legend(x = "topleft",
         legend = c("Uncensored", "Censored"), col=c("darkolivegreen4","blue"),
         pch=c(2,3),cex=0.5,xpd = TRUE,bty="L",horiz=TRUE)

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

