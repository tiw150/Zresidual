#' A function to draw scatter plot of a z-residual
#'
#' @param Zresidual A Z-residual.
#' @param outlier.set A list of parameters available in symbols() and text().
#' @export plot.zresid
#' @import stringr
#'

plot.zresid <- function(Zresidual,irep=1:ncol(Zresidual),ylab = "Z-Residual", normality.test = c("SW", "AOV", "BL"),
                        k.test=10, X = c("index", "covariate", "lp"),
                        main.title = paste("Z-residual Scatterplot -", attr(Zresidual, "type")),
                        outlier.return = FALSE, outlier.value = 3.5, category=NULL, outlier.set = list(), ...)
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
  
  args <- list(...)
  var.call <- match.call()
  
  unique.cats <- NULL
  default.legend.title <- NULL
  
  type <- attr(Zresidual, "type")
  zero.id <- attr(Zresidual, "zero_id")
  
  # Override defaults if category is provided
  if (!is.null(category)) {
    unique.cats <- unique(category)
    default.legend.title <- deparse(substitute(category))
    if (!is.null(args[["col"]])) {
      col <- args[["col"]]
    } else {
      col <- if (length(unique.cats) == 1) "red" else rainbow(100)[1:length(unique.cats)]
    }
    if (!is.null(args[["pch"]])) {
      pch <- args[["pch"]]
    } else {
      pch <- if (length(unique.cats) == 1) 1 else c(1:25)[1:length(unique.cats)]
    }
  } else{
    
    if(type == 'hurdle'){
      unique.cats <- c("zero", "count")
      col <- c("red", "blue")[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
      pch <- c(1,3)[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
    } else if(type %in% c('count')){
      unique.cats <- type
      col <- "blue"
      pch <- 3
    } else if(type %in% c('zero')){
      unique.cats <- type
      col <- "red"
      pch <- 1
    } else {
      col <- "red"
      pch <- 1
    }
    
    if(!is.null(args[["col"]])){
      col <-  args[["col"]]
      #default.legend.title <- deparse(var.call[["col"]])
      unique.cats <- if(is.symbol(var.call["col"])) deparse(var.call[["col"]]) else unique(args[["col"]])
    }
    if(!is.null(args[["pch"]])) {
      pch <- args[["pch"]]
      #default.legend.title <- deparse(var.call[["pch"]])
      unique.cats <- if(is.symbol(var.call["pch"])) deparse(var.call[["pch"]]) else unique(args[["pch"]])
    }
  }
  
  
  test.legend <- NULL
  for (j in irep) {
    
    par(mar = c(5, 4, 4, 6) + 0.1)
    
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
    
    # Test Legend
    test.pv <- NULL
    if (X != "index") {
      test.list <- c("SW"="sw", "AOV" = "aov", "BL" = "bartlett")
      #if(!missing(normality.test)){
      for (a in normality.test) {
        test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
        test.pv <- c(test.pv, paste(a, "-", sprintf("%3.2f",test[j])))
      }
      
      #legend.text.width <- max(strwidth(c(test.pv, unique.cats)))
      
      test.legend <- list(
        cex = 0.6,
        bty = "n",
        text.width= max(convertWidth(stringWidth(c(test.pv, unique.cats)), "inches", valueOnly = TRUE)),
        xpd=TRUE,
        title = "P-value",
        title.font = 2,
        legend = test.pv
      )
    }
    
    default.legend <- list(
      legend = unique.cats,
      col = unique(col),
      pch = unique(pch),
      cex = 0.6,
      text.width = max(convertWidth(stringWidth(c(test.pv, unique.cats)), "inches", valueOnly = TRUE)),
      xpd = TRUE,
      bty = "n",
      title = if(!hasArg("title")) default.legend.title else title,
      horiz = FALSE
    )
    
    
    # Add user-provided `...` arguments, overwriting defaults if necessary
    legend.args <- modifyList(default.legend, args[!names(args) %in% c("col", "pch")])
    legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
    
    default.outlier <- list(
      pos = 4,
      labels = id.outlier,
      cex = 0.8,
      col = "blue",
      add = T,
      inches = F,
      circles = rep((par("usr")[2]-par("usr")[1])*0.05, length(id.outlier)),
      fg = "blue"
    )
    
    # Add user-provided `...` arguments, overwriting defaults if necessary
    outlier.args <- modifyList(default.outlier, outlier.set)
    
    # Extract necessary arguments seperately for text() and symbols()
    text.args <- outlier.args[names(outlier.args) %in% formalArgs(text.default)]
    symbols.args <- outlier.args[names(outlier.args) %in% formalArgs(symbols)]
    
    
    default.plot <- modifyList(list(col = col,
                                    pch = pch,
                                    ylab = ylab,
                                    ylim = c(-ylim0, ylim0 + 1),
                                    main = main.title),
                               args)
    if (X == "index") {
      do.call(plot.default, c(list(x=Zresidual[,j], xlab = "Index"), default.plot))
      plot_limits <- par("usr")
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          do.call(symbols, c(list(id.outlier, Zresidual[,j][id.outlier]), symbols.args))
          do.call(text, c(list(id.outlier, Zresidual[,j][id.outlier]), text.args))
        }
      }
      if (length(legend.args$legend) > 0) do.call(legend, c(list(x=plot_limits[2], y = plot_limits[4]), legend.args))
      
    }
    
    if (X == "lp") {
      
      fitted.value <- attr(Zresidual, "linear.pred")
      
      do.call(plot, c(list(x=fitted.value, y=Zresidual[,j], xlab = "Linear Predictor"), default.plot))
      plot_limits <- par("usr")
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          do.call(symbols, c(list(fitted.value[id.outlier], Zresidual[,j][id.outlier]), symbols.args))
          do.call(text, c(list(fitted.value[id.outlier], Zresidual[,j][id.outlier]), text.args))
        }
      }
      if (length(legend.args$legend) > 0) do.call(legend, c(list(x=plot_limits[2], y = plot_limits[4]), legend.args))
      if(!is.null(test.legend)) do.call(legend, c(list(x=plot_limits[2], y = mean(plot_limits[3:4])), test.legend))
      
    }
    
    if (X != "index" & X != "lp") {
      
      fitted.value <- attr(Zresidual, "covariates")
      
      if(X == "covariate"){
        i<-1
        cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:", variable.names(fitted.value))
      } else if(X %in% variable.names(fitted.value)){
        cov.name<-variable.names(fitted.value)
        i<- which(cov.name==X)
      } else{stop("X must be the one of covariate name.") }
      
      do.call(plot, c(list(x=fitted.value[,i], y=Zresidual[,j], xlab = colnames(fitted.value)[i]), default.plot))
      plot_limits <- par("usr")
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          do.call(symbols, c(list(fitted.value[, i][id.outlier], Zresidual[,j][id.outlier]), symbols.args))
          do.call(text, c(list(fitted.value[,i][id.outlier], Zresidual[,j][id.outlier]), text.args))
        }
      }
      if (length(legend.args$legend) > 0) do.call(legend, c(list(x=plot_limits[2], y = plot_limits[4]), legend.args))
      if(!is.null(test.legend)) do.call(legend, c(list(x=plot_limits[2], y = mean(plot_limits[3:4])), test.legend))
    }
    
    if (length(id.infinity) > 0L) {
      text(id.infinity + 5,
           Zresidual[,j][id.infinity],
           labels = value.notfinite,
           col = 2)
    }
    hlines <- c(1.96, 3)
    hlines2 <- -hlines
    abline(h = c(hlines, hlines2),
           lty = 3,
           col = "grey")
    if (outlier.return & !identical(id.outlier, integer(0))) {
      cat("Outlier Indices(", attr(Zresidual, "type"), "):", id.outlier, "\n", sep = " ")
      invisible(list(outliers = id.outlier))
    }
    
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}