#' A function to draw box plot of a z-residual
#'
#' @param Zresidual A Z-residual.
#' @export boxplot.zresid

boxplot.zresid <- function(Zresidual,irep=1,
                           X = c("lp", "covariate"),
                           num.bin = 10,
                           main.title=paste("Z-residual Boxplot -", attr(Zresidual, "type")),
                           outlier.return = FALSE, outlier.value = 3.5,
                           ...)
{
  
  sign.na <- function(x)
  {
    sign.x <- sign(x)
    sign.x[is.infinite(x)] <- 1
    sign.x
  }
  
  as.character.na <- function(x)
  {
    label.x <- as.character(x)
    label.x[is.infinite(x)] <- "Inf"
    label.x
  }
  
  if (missing(X))
    X = "lp"
  
  for (j in irep) {
    
    
    type <- attr(Zresidual, "type")
    zero.id <- attr(Zresidual, "zero_id")
    id.nan <- which(is.nan(Zresidual[,j]))
    id.infinity <- which (is.infinite(Zresidual[,j]))
    id.outlier <- which(abs(Zresidual[,j]) > outlier.value | is.infinite(Zresidual[,j]))
    
    # Converting Inf/-Inf to maximum finite value
    if (length(id.infinity) > 0L) {
      value.notfinite <- as.character.na(Zresidual[,j][id.infinity])
      max.non.infinity <- max(abs(Zresidual[,j]), na.rm = T)
      Zresidual[,j][id.infinity] <- sign.na(Zresidual[,j][id.infinity]) * (max.non.infinity + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    
    if(length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
    
    ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[,j]), na.rm = T))
    
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
                   aov.test.zresid(Zresidual, X=X, k.anova = num.bin)[j]
                 )),
          paste0("Z-BL p-value = ",
                 sprintf(
                   "%3.2f",
                   bartlett.test.zresid(Zresidual, X=X, k.bl = num.bin)[j]
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
                   aov.test.zresid(Zresidual, X= X, k.anova = num.bin)[j]
                 )),
          paste0("Z-BL p-value = ",
                 sprintf(
                   "%3.2f",
                   bartlett.test.zresid(Zresidual, X=X, k.bl = num.bin)[j]
                 ))
        ),
        cex = 0.6,
        horiz = TRUE
      )
    }
    
    if (outlier.return)
    {
      cat("Outlier Indices:", id.outlier, "\n")
      invisible(list(outliers = id.outlier))
    }
    
  }
}