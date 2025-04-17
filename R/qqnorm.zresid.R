#' A function to draw qq plot of a z-residual
#'
#' @param Zresidual A z-residual.
#' @param point.clr A vector of colors of points.
#' @param point.type A vector of type of points.
#' @param legend.settings A list of parameters available in legend().
#' @param outlier.settings A list of parameters available in symbols() and text().
#' @param diagnosis.test diagnosis test
#' @param X.anova X variable to run if ANOVA is selected in diagnosis.test. Possible c("lp", "covariate")
#' @param k.anova Number of bins if X.anova is 'covariates'.
#' @export qqnorm.zresid
#'
qqnorm.zresid <- function (Zresidual, irep=1, diagnosis.test = c("SW", "ANOVA"), X.anova = c("lp", "covariate"),
                           k.anova=10, main.title = paste("Normal Q-Q Plot -", attr(Zresidual, "type")),
                           xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                           outlier.return=FALSE, outlier.value = 3.5, legend.settings = list(), ...)
{
  #i<-index
  
  if(missing(diagnosis.test)) diagnosis.test <- "SW"
  if(missing(X.anova) && diagnosis.test == "ANOVA") X.anova <- "lp"
  # test <- list("SW" = shapiro.test,
  #              "ANOVA" = aov.test.zresid)
  if(diagnosis.test == "ANOVA") test <- aov.test.zresid(Zresidual, X.anova, k.anova) else test <- sw.test.zresid(Zresidual)
  
  for (i in irep) {
    
    type <- attr(Zresidual, "type")
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
    default.legend <- list(
      x = "topleft",
      cex = 0.6,
      bty = "n",
      horiz = T,
      #text.width = strwidth(paste0("Z-", diagnosis.test, " p-value = ",sprintf("%3.2f",test.pv)), units = "inches"),
      legend = paste0("Z-", diagnosis.test, " p-value = ",sprintf("%3.2f",test.pv))
    )
    
    # Add user-provided `...` arguments, overwriting defaults if necessary
    legend.args <- modifyList(default.legend, legend.settings)
    
    if(max(abs(Zresidual[,i]), na.rm = T)<6){
      qqnorm(Zresidual[,i],main=main.title,xlab=xlab, ylab=ylab)
      qqline(Zresidual[,i],col=1)
      abline(a=0,b=1,col=3)
      legend.box(legend.args)
    }
    if(max(abs(Zresidual[,i]), na.rm = T)>=6){
      max.values<-which(Zresidual[,i]>=6)
      min.values<-which(Zresidual[,i] < -6)
      if(identical(min.values, integer(0))){
        gap=c(6,max(abs(Zresidual[,i]), na.rm = T))
        gapsize <- gap[2] - gap[1]
        ylim0<-c(min(Zresidual[,i][-max.values], na.rm = T),7)
        xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T),
                 max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T))
        plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-max.values],
             qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-max.values],
             xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
        axis(2, at=7 ,labels=round(max(Zresidual[,i], na.rm = T),1))
        qqline(Zresidual[,i],col=1)
        abline(a=0,b=1,col=3)
        axis.break(2,6,style="gap")
        axis.break(2, 6, breakcol="snow", style="gap")
        axis.break(2, 6,breakcol="black", style="slash")
        axis.break(4, 6,breakcol="black", style="slash")
        points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
               pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize+1,6.5),col="red")
        legend.box(legend.args)
      }
      
      if(identical(max.values, integer(0))){
        gap=c(min(Zresidual[,i], na.rm = T),6)
        gapsize2 <- gap[2] + gap[1]
        ylim0<-c(-7 ,max(Zresidual[,i][-min.values], na.rm = T))
        xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T),
                 max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T))
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
        
        legend.box(legend.args)
      }
      
      if(!identical(max.values, integer(0)) && !identical(min.values, integer(0))){
        gap1=c(6,max(Zresidual[,i]))
        gapsize1 <- gap1[2] - gap1[1]
        gap2=c(min(Zresidual[,i]),6)
        gapsize2 <- gap2[2] + gap2[1]
        ylim0<-c(-7,7)
        xlim0<-c(min(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T),
                 max(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]], na.rm = T))
        plot(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][-c(max.values,min.values)],
             qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][-c(max.values,min.values)],
             xlim=xlim0,ylim=ylim0,main=main.title,xlab=xlab, ylab=ylab)
        qqline(Zresidual[,i],col=1)
        abline(a=0,b=1,col=3)
        axis(2, at=7 ,labels=round(max(Zresidual[,i], na.rm = T),1))
        axis.break(2,6,style="gap")
        axis.break(2, 6, breakcol="snow", style="gap")
        axis.break(2, 6,breakcol="black", style="slash")
        axis.break(4, 6,breakcol="black", style="slash")
        points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][max.values],
               pmax(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][max.values]-gapsize1+1,6.5),col="red")
        
        axis(2, at= -7 ,labels=round(min(Zresidual[,i], na.rm = T),1))
        axis.break(2,-6.1,style="gap")
        axis.break(2, -6.1, breakcol="snow", style="gap")
        axis.break(2, -6.1,breakcol="black", style="slash")
        axis.break(4, -6.1, breakcol="black", style="slash")
        points(qqnorm(Zresidual[,i],plot.it = FALSE)[["x"]][min.values],
               pmin(qqnorm(Zresidual[,i],plot.it = FALSE)[["y"]][min.values]-gapsize2-1,-6.5),col="red")
        legend.box(legend.args)
      }
      
    }
    
    if(isTRUE(outlier.return)){
      
      if(identical(id.outlier, integer(0))){
        #return(invisible(NULL))
        next
      } else {
        symbols((qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier],
                (qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier],
                circles=rep(0.1,length(id.outlier)),
                fg=rep('red',length(id.outlier)),
                add=T, inches=F)
        text((qqnorm(Zresidual[,i],plot.it=FALSE)$x)[id.outlier],
             (qqnorm(Zresidual[,i],plot.it=FALSE)$y)[id.outlier],
             pos=1,label = id.outlier,
             cex = 0.8,col="red")
      }
    }
    if(outlier.return){
      cat("Outlier Indices (", type, ") :", id.outlier, "\n")
      invisible(list(outliers=id.outlier))
    }
    
  }
}