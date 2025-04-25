#' A function to calculate ANOVA of Zresidual
#'
#' @param Zresidual A z-residual.
#' @param X Linear predictor or covariate
#' @param k.anova Number of bins if applicable
#' @export

aov.test.zresid <- function (Zresidual,X = c("lp", "covariate"), k.anova=10)
{
  if (missing(X)) X = "lp"
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
  return(aov.pv)
}

#' A function to draw scatter plot of a z-residual
#'
#' @param Zresidual A Z-residual.
#' @param outlier.set A list of parameters available in symbols() and text().
#' @param xlab Optional label for the x-axis. If NULL (default), the function
#'   determines the label based on 'X'. Accepts LaTeX code enclosed in 'tex()'
#'   for enhanced formatting using the latex2exp package.
#' @export plot.zresid
#' @import stringr

plot.zresid <- function(Zresidual, irep = 1:ncol(Zresidual), ylab = "Z-Residual",
                        normality.test = c("SW", "AOV", "BL"), k.test = 10,
                        X = c("index", "covariate", "lp"),
                        main.title = paste("Z-residual Scatterplot -",
                                           attr(Zresidual, "type")),
                        outlier.return = FALSE, outlier.value = 3.5,
                        category = NULL, outlier.set = list(), xlab = NULL,
                        ...) {

  sign.na <- function(x) {
    sign.x <- sign(x)
    sign.x[is.infinite(x)] <- 1
    sign.x
  }
  as.character.na <- function(x) {
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
  test.pv <- NULL # Initialize test.pv outside the loop

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
  } else {
    if (type == "hurdle") {
      unique.cats <- c("zero", "count")
      col <- c("red", "blue")[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
      pch <- c(1, 3)[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
    } else if (type %in% c("count")) {
      unique.cats <- type
      col <- "blue"
      pch <- 3
    } else if (type %in% c("zero")) {
      unique.cats <- type
      col <- "red"
      pch <- 1
    } else {
      col <- "red"
      pch <- 1
    }
    if (!is.null(args[["col"]])) {
      col <- args[["col"]]
      unique.cats <- if (is.symbol(var.call["col"]))
        deparse(var.call[["col"]]) else unique(args[["col"]])
    }
    if (!is.null(args[["pch"]])) {
      pch <- args[["pch"]]
      unique.cats <- if (is.symbol(var.call["pch"]))
        deparse(var.call[["col"]]) else unique(args[["pch"]])
    }
  }

  for (j in irep) {
    par(mar = c(5, 4, 4, 6) + 0.1)
    id.nan <- which(is.nan(Zresidual[, j]))
    id.infinity <- which(is.infinite(Zresidual[, j]))
    id.outlier <- which(abs(Zresidual[, j]) > outlier.value |
                          is.infinite(Zresidual[, j]))
    if (length(id.infinity) > 0L) {
      value.notfinite <- as.character.na(Zresidual[, j][id.infinity])
      max.non.infinity <- max(abs(Zresidual[, j]), na.rm = TRUE)
      Zresidual[, j][id.infinity] <- sign.na(Zresidual[, j][id.infinity]) *
        (max.non.infinity + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
    ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[, j]), na.rm = TRUE))
    if (X != "index") {
      test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
      current_test_pv <- NULL # Initialize for each iteration
      for (a in normality.test) {
        test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
        current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
      }
      if (length(irep) == 1) { # Only create test.legend if plotting one column
        test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
                            cex = 0.6, bty = "n", xpd = TRUE, adj = c(0, 0.5))
      }
    }
    default.legend <- list(legend = unique.cats, col = unique(col),
                           pch = unique(pch), cex = 0.6, xpd = TRUE, bty = "n",
                           title = if (!hasArg("title")) default.legend.title else title,
                           horiz = FALSE, y.intersp = 1)
    legend.args <- modifyList(default.legend, args[!names(args) %in% c("col", "pch")])
    legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
    default.outlier <- list(pos = 4, labels = id.outlier, cex = 0.8,
                            col = "blue", add = TRUE, inches = FALSE,
                            circles = rep((par("usr")[2] - par("usr")[1]) * 0.05,
                                          length(id.outlier)),
                            fg = "blue")
    outlier.args <- modifyList(default.outlier, outlier.set)
    text.args <- outlier.args[names(outlier.args) %in% formalArgs(text.default)]
    symbols.args <- outlier.args[names(outlier.args) %in% formalArgs(symbols)]

    current_xlab <- if (!is.null(xlab)) {
      if (startsWith(as.character(xlab), "tex(")) {
        if (requireNamespace("latex2exp", quietly = TRUE)) {
          latex_string <- sub("tex\\((.*)\\)", "\\1", as.character(xlab))
          latex2exp::TeX(latex_string)
        } else {
          warning("The 'latex2exp' package is not installed. LaTeX formatting for xlab will not be used.")
          sub("tex\\((.*)\\)", "\\1", as.character(xlab))
        }
      } else {
        xlab
      }
    } else if (X == "index") {
      "Index"
    } else if (X == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      if (!is.null(fitted.value)) "Linear Predictor" else "Fitted Value"
    } else if (X != "index" & X != "lp") {
      fitted.value <- attr(Zresidual, "covariates")
      if (!is.null(fitted.value)) {
        cov.name <- variable.names(fitted.value)
        i <- if (X == "covariate") 1 else which(cov.name == X)
        colnames(fitted.value)[i]
      } else {
        if (X == "covariate") "Covariate" else X
      }
    } else {
      "x"
    }

    default.plot <- modifyList(list(col = col, pch = pch, ylab = ylab,
                                    ylim = c(-ylim0, ylim0 + 1),
                                    main = main.title, xlab = current_xlab,
                                    font.lab = 2), args)
    if (X == "index") {
      do.call(plot.default, c(list(x = Zresidual[, j]), default.plot))
      plot_limits <- par("usr")
      do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          do.call(symbols, c(list(id.outlier, Zresidual[, j][id.outlier]), symbols.args))
          do.call(text, c(list(id.outlier, Zresidual[, j][id.outlier]), text.args))
        }
      }
    }
    if (X == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      do.call(plot, c(list(x = fitted.value, y = Zresidual[, j]), default.plot))
      plot_limits <- par("usr")
      do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          do.call(symbols, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), symbols.args))
          do.call(text, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), text.args))
        }
      }
    }
    if (X != "index" & X != "lp") {
      fitted.value <- attr(Zresidual, "covariates")
      if (!is.null(fitted.value)) {
        cov.name <- variable.names(fitted.value)
        i <- if (X == "covariate") 1 else which(cov.name == X)
        do.call(plot, c(list(x = fitted.value[, i], y = Zresidual[, j]), default.plot))
        plot_limits <- par("usr")
        do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
        if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
        if (isTRUE(outlier.return)) {
          if (!identical(id.outlier, integer(0))) {
            do.call(symbols, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), symbols.args))
            do.call(text, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), text.args))
          }
        }
      }
    }

    if (length(id.infinity) > 0L) {
      text(id.infinity + 5, Zresidual[, j][id.infinity],
           labels = value.notfinite, col = 2)
    }
    hlines <- c(1.96, 3)
    hlines2 <- -hlines
    abline(h = c(hlines, hlines2), lty = 3, col = "grey")
    if (outlier.return & !identical(id.outlier, integer(0))) {
      cat("Outlier Indices(", attr(Zresidual, "type"), "):", id.outlier, "\n",
          sep = " ")
      invisible(list(outliers = id.outlier))
    }
  }
  par(mar = c(5, 4, 4, 2) + 0.1)
}
