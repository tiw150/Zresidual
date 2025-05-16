#' A function to draw scatter plot of a Z-residual
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
                        main.title = ifelse(is.null(attr(Zresidual, "type")),
                                            "Z-residual Scatterplot",
                                            paste("Z-residual Scatterplot -",
                                                  attr(Zresidual, "type"))),
                        outlier.return = TRUE, outlier.value = 3.5,
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
  } else if (is.null(type)){
    unique.cats <- c("Uncensored", "Censored")
    censored <- attr(Zresidual, "censored.status")
    col<- c("blue","red")[censored+1]
    pch<- c(3,2)[censored+1]
  } else {
    if (type == "hurdle") {
      unique.cats <- c("count","zero")
      col <- c("blue", "red")[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
      pch <- c(3,2)[seq_along(Zresidual) %in% attr(Zresidual, "zero_id") + 1]
    } else if (type %in% c("count")) {
      unique.cats <- type
      col <- "blue"
      pch <- 3
    } else if (type %in% c("zero")) {
      unique.cats <- type
      col <- "red"
      pch <- 2
    } else {
      col <- "red"
      pch <- 2
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
 #   plot.new()
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
    test.legend <- NULL
    if (X != "index") {
      test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
      current_test_pv <- NULL # Initialize for each iteration
      for (a in normality.test) {
        test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
        current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
      }

      test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
                            cex = 0.6, bty = "n", xpd = TRUE, adj = c(0, 0.5))
    }
    default.legend <- list(legend = unique.cats, col = unique(col),
                           pch = unique(pch), cex = 0.6, xpd = TRUE, bty = "n",
                           title = if (!hasArg("title")) default.legend.title else title,
                           horiz = FALSE, y.intersp = 1)
    legend.args <- modifyList(default.legend, args[!names(args) %in% c("col", "pch")])
    legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
    default.outlier <- list(pos = 4, labels = id.outlier, cex = 0.8,
                            col = "darkolivegreen4", add = TRUE, inches = FALSE,
                            circles = rep((par("usr")[2] - par("usr")[1]) * 0.03,
                                          length(id.outlier)),
                            fg = "darkolivegreen4", font = 2)
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
      plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
      if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

      do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          # Recalculate circles based on current par("usr")
          outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
          # Update symbols.args with the recalculated circles
          symbols.args$circles <- outlier.circles
          symbols.args$add <- NULL
          effective_symbols_args <- c(list(x = id.outlier, y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
          do.call(symbols, effective_symbols_args)
          # Adjust text position for the last outlier
          text_x <- id.outlier
          text_y <- Zresidual[, j][id.outlier]
          text_pos <- outlier.args$pos # Default position

          y_median <- median(Zresidual[, j], na.rm = TRUE)

          for (i in seq_along(id.outlier)) {
            if (!is.na(text_y[i])) {
              if (text_y[i] < y_median) {
                text_pos[i] <- 3 # Position above for lower outliers
              } else {
                text_pos[i] <- 1 # Position below for upper outliers
              }
            }
          }
          text.args_no_pos <- text.args[names(text.args) != "pos"]
          do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
          #do.call(text, c(list(id.outlier, Zresidual[, j][id.outlier]), text.args))
        }
      }
    }
    if (X == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      do.call(plot, c(list(x = fitted.value, y = Zresidual[, j]), default.plot))

      plot_limits <- par("usr")
      plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
      if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

      do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {
          outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
          symbols.args$circles <- outlier.circles
          symbols.args$add <- NULL
          effective_symbols_args <- c(list(x = fitted.value[id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
          do.call(symbols, effective_symbols_args)
          #do.call(symbols, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), symbols.args))

          text_x <- fitted.value[id.outlier]
          text_y <- Zresidual[, j][id.outlier]
          text_pos <- outlier.args$pos # Default position
          y_median <- median(Zresidual[, j], na.rm = TRUE)
          for (i in seq_along(id.outlier)) {
            if (!is.na(text_y[i])) {
              if (text_y[i] < y_median) {
                text_pos[i] <- 3 # Position above for lower outliers
              } else {
                text_pos[i] <- 1 # Position below for upper outliers
              }
            }
          }
          text.args_no_pos <- text.args[names(text.args) != "pos"]
          do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
        #  do.call(text, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), text.args))
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
        plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
        if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

        do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
        if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2] - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
        if (isTRUE(outlier.return)) {
          if (!identical(id.outlier, integer(0))) {
            # Recalculate circles based on current par("usr")
            outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
            # Update symbols.args with the recalculated circles
            symbols.args$circles <- outlier.circles
            symbols.args$add <- NULL
            effective_symbols_args <- c(list(x = fitted.value[, i][id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
            do.call(symbols, effective_symbols_args)
            #do.call(symbols, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), symbols.args))
            text_x <- fitted.value[, i][id.outlier]
            text_y <- Zresidual[, j][id.outlier]
            text_pos <- outlier.args$pos
            y_median <- median(Zresidual[, j], na.rm = TRUE)
            for (i in seq_along(id.outlier)) {
              if (!is.na(text_y[i])) {
                if (text_y[i] < y_median) {
                  text_pos[i] <- 3 # Position above for lower outliers
                } else {
                  text_pos[i] <- 1 # Position below for upper outliers
                }
              }
            }
            text.args_no_pos <- text.args[names(text.args) != "pos"]
            do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
           # do.call(text, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), text.args))
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

