#' Plot Z-Residuals for Bayesian and Frequentist Count / Hurdle / Zero / Survival Models
#'
#' @description
#' Produces diagnostic scatterplots of Z-residuals for a wide range of model
#' types, including hurdle models, zero models, count models, and survival models.
#' This function is designed to be compatible with Z-residual matrices generated
#' from Bayesian models (e.g., **brms**, Stan-based models) as well as classical models.
#'
#' The function supports:
#' - Plotting residuals against index, covariates, linear predictors, or user-specified vectors.
#' - Visual outlier detection with customizable coloring, emphasis, and labels.
#' - Automatic handling of censored/un-censored and hurdle count/zero classifications.
#' - Multiple normality tests (Shapiro-Wilk, ANOVA, Bartlett) with per-iteration reporting.
#' - Extensive customization through graphical parameters.
#'
#' @usage
#' \method{plot}{zresid}(
#'   x,
#'   irep = 1:ncol(x),
#'   ylab = "Z-Residual",
#'   normality.test = c("SW", "AOV", "BL"),
#'   k.test = 10,
#'   x_axis_var = c("index", "covariate", "lp"),
#'   main.title = ifelse(is.null(attr(x, "type")),
#'                       "Z-residual Scatterplot",
#'                       paste("Z-residual Scatterplot -", attr(x, "type"))),
#'   outlier.return = TRUE,
#'   outlier.value = 3.5,
#'   category = NULL,
#'   outlier.set = list(),
#'   xlab = NULL,
#'   my.mar = c(5, 4, 4, 6) + 0.1,
#'   add_lowess = FALSE,
#'   ...
#' )
#'
#' @param x
#' A numeric matrix of Z-residuals (dimensions: *n* × *m*) of class `"zresid"`.
#' Attributes used by the function include: `"type"`, `"zero_id"`, `"censored.status"`,
#' `"covariates"`, and `"linear.pred"`.
#'
#' @param irep
#' A vector specifying the residual column(s) (iterations) to plot. Defaults to
#' all columns of **`x`**.
#'
#' @param ylab Label for the y-axis.
#'
#' @param normality.test
#' A character vector specifying normality tests applied per iteration: `"SW"`,
#' `"AOV"`, or `"BL"`. Helper functions (e.g., `sw.test.zresid`) must exist and
#' accept inputs (`x`, `x_axis_var`, `k.test`).
#'
#' @param k.test
#' Bin size for normality tests (used for grouping continuous predictors).
#'
#' @param x_axis_var
#' Specifies the x-axis values. Options include: `"index"`, `"covariate"`, `"lp"`,
#' a covariate name present in `attr(x, "covariates")`, or a numeric vector of length *n*.
#'
#' @param main.title
#' Title of the plot. Defaults to a type-based informative title.
#'
#' @param outlier.return
#' If `TRUE`, outliers are printed to console and returned invisibly.
#'
#' @param outlier.value
#' Threshold above which a residual is flagged as an outlier. Defaults to `3.5`.
#'
#' @param category
#' Optional vector categorizing observations (length *n*). Used for coloring and
#' shaping points in scatterplots.
#'
#' @param outlier.set
#' A named list of arguments passed to `symbols()` and `text()` for marking and
#' labeling outliers. Overrides defaults.
#'
#' @param xlab
#' Label for the x-axis. May include LaTeX syntax using the form `tex("...")`,
#' which will be interpreted via **latex2exp** (if installed).
#'
#' @param my.mar
#' A numeric vector passed to `par(mar=...)` to adjust plot margins.
#'
#' @param ...
#' Additional graphical arguments passed to `plot()`, `legend()`, `symbols()`,
#' and `text()`.
#'
#' @details
#' This function employs S3 method dispatch. The input `x` must be an object of
#' class `"zresid"`.
#'
#' ## **Model-Type-Specific Behavior**
#'
#' ### **Hurdle models**
#' - Zeros and counts are colored differently (`red` vs `blue`).
#'
#' ### **Survival models**
#' - Censored and uncensored observations are visually separated.
#'
#' ## **Outlier Detection**
#' Outliers are defined as:
#'
#' \deqn{|Z| > \mbox{outlier.value} \hspace{1em} \mbox{or non-finite values}}
#'
#' They are marked, labeled, and returned to the user if `outlier.return = TRUE`.
#'
#' @return
#' Invisibly returns (when `outlier.return = TRUE`) a list:
#' \item{outliers}{Vector of outlier indices}
#' Otherwise returns `NULL`. Always produces a scatterplot as its primary output.
#'
#' @note
#' This function modifies graphical parameters (`par(mar=...)`) during execution
#' and resets them at the end.
#'
#' @exportS3Method graphics::plot zresid
#' @export plot.zresid
plot.zresid <- function(x, irep = 1:ncol(x), ylab = "Z-Residual",
                        normality.test = c("SW", "AOV", "BL"), k.test = 10,
                        x_axis_var = c("index", "covariate", "lp"),
                        main.title = ifelse(is.null(attr(x, "type")),
                                            "Z-residual Scatterplot",
                                            paste("Z-residual Scatterplot -",
                                                  attr(x, "type"))),
                        outlier.return = TRUE, outlier.value = 3.5,
                        category = NULL, outlier.set = list(), xlab = NULL,
                        my.mar=c(5,4,4,6)+0.1, add_lowess = FALSE, ...) {
  Zresidual <- x
  X <- x_axis_var

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
  
  add_lowess_line <- function(xv, yv, args) {
    if (!isTRUE(add_lowess)) return(invisible(NULL))
    
    log_opt <- args[["log"]]
    log_x <- !is.null(log_opt) && grepl("x", log_opt)
    
    ok <- is.finite(xv) & is.finite(yv)
    if (log_x) ok <- ok & (xv > 0)
    
    if (sum(ok) < 3L) return(invisible(NULL))
    
    if (log_x) {
      lx <- log10(xv[ok])
      lw <- stats::lowess(x = lx, y = yv[ok])
      graphics::lines(x = 10^lw$x, y = lw$y, col = "red", lwd = 3)
    } else {
      lw <- stats::lowess(x = xv[ok], y = yv[ok])
      graphics::lines(lw$x, lw$y, col = "red", lwd = 3)
    }
    
    invisible(NULL)
  }

  var.call <- match.call()
  unique.cats <- NULL
  default.legend.title <- NULL
  type <- attr(Zresidual, "type")
  zero.id <- attr(Zresidual, "zero_id")
  test.pv <- NULL # Initialize test.pv outside the loop

  ## --- choose colors/pchs & legend defaults (fixed) ---
  legend_colors <- NULL
  legend_pchs   <- NULL

  if (!is.null(category)) {
    # map per-observation by category
    unique.cats <- unique(category)
    default.legend.title <- deparse(substitute(category))

    pal  <- if (!is.null(args[["col"]])) args[["col"]] else
      if (length(unique.cats) == 1) "red" else rainbow(max(length(unique.cats), 2))[seq_along(unique.cats)]
    pchs <- if (!is.null(args[["pch"]])) args[["pch"]] else
      if (length(unique.cats) == 1) 1 else (1:25)[seq_along(unique.cats)]

    col <- pal[match(category, unique.cats)]
    pch <- pchs[match(category, unique.cats)]
    legend_colors <- pal
    legend_pchs   <- pchs

  } else if (is.null(type)) {
    unique.cats <- NULL
    col <- "black"
    pch <- 1
    legend_colors <- NULL
    legend_pchs   <- NULL

  } else if (type == "survival") {
    unique.cats <- c("Uncensored", "Censored")
    censored <- attr(Zresidual, "censored.status")
    col <- c("blue","red")[censored + 1]
    pch <- c(3,2)[censored + 1]
    legend_colors <- c("blue", "red")
    legend_pchs   <- c(3, 2)

  } else {
    if (type == "hurdle") {
      legend_labels <- c("count", "zero")
      legend_colors <- c("blue", "red")
      legend_pchs   <- c(3, 2)
      n <- NROW(Zresidual)
      is_zero <- seq_len(n) %in% attr(Zresidual, "zero_id")
      col <- c("blue", "red")[is_zero + 1]
      pch <- c(3, 2)[is_zero + 1]
      unique.cats <- legend_labels

    } else if (type %in% c("count")) {
      unique.cats <- type
      col <- "blue"; pch <- 3
      legend_colors <- "blue"; legend_pchs <- 3

    } else if (type %in% c("zero")) {
      unique.cats <- type
      col <- "red"; pch <- 2
      legend_colors <- "red"; legend_pchs <- 2

    } else {
      col <- "red"; pch <- 2
      legend_colors <- "red"; legend_pchs <- 2
    }

    # If user passed col/pch in ..., honor them without altering legend labels
    if (!is.null(args[["col"]])) col <- args[["col"]]
    if (!is.null(args[["pch"]])) pch <- args[["pch"]]
  }

  ## --- build legend args only when needed (fixed) ---
  legend.args <- NULL
  if (!is.null(unique.cats) && length(unique.cats) > 0 &&
      !is.null(legend_colors) && length(legend_colors) > 0 &&
      !is.null(legend_pchs)   && length(legend_pchs)   > 0) {
    default.legend <- list(
      legend = unique.cats,
      col = legend_colors,
      pch = legend_pchs,
      cex = 1, xpd = TRUE, bty = "n",
      title = if (!hasArg("title")) default.legend.title else title,
      horiz = FALSE, y.intersp = 1
    )
    # Only allow args that legend() knows, but keep title if present
    user_legend_overrides <- args[!names(args) %in% c("col", "pch")]
    legend.args <- modifyList(default.legend, user_legend_overrides)
    legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
  }
  choices <- c("index", "covariate", "lp")
  n_obs   <- NROW(Zresidual)
  is_uservec <- (length(X) == n_obs)
  is_keyword <- is.character(X) && length(X) == 1 && X %in% choices
  is_covname <- is.character(X) && length(X) == 1 && !is_keyword

  if (!is_uservec && !is_keyword && !is_covname) {
    stop("X must be: (1) a length-n vector, or (2) one of 'index','lp','covariate', or (3) a covariate name present in attr(Zresidual,'covariates').")
  }


  for (j in irep) {
    #   plot.new()
    par(mar = my.mar)
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

    if (is_uservec) {
      xvec <- X

      {
        test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
        current_test_pv <- NULL
        for (a in normality.test) {
          testfun <- get(paste0(test.list[a], ".test.zresid"))
          tt <- try(testfun(Zresidual, xvec, k.test), silent = TRUE)
          pv_str <- if (inherits(tt, "try-error")) "NA" else sprintf("%3.2f", tt[j])
          current_test_pv <- c(current_test_pv, paste(a, "-", pv_str))
        }
        test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
                            cex = 1, bty = "n", xpd = TRUE, adj = c(0, 0.5))
      }

      current_xlab <- if (!is.null(xlab)) {
        if (startsWith(as.character(xlab), "tex(")) {
          if (requireNamespace("latex2exp", quietly = TRUE)) {
            latex_string <- sub("tex\\((.*)\\)", "\\1", as.character(xlab))
            latex2exp::TeX(latex_string)
          } else {
            warning("The 'latex2exp' package is not installed. LaTeX formatting for xlab will not be used.")
            sub("tex\\((.*)\\)", "\\1", as.character(xlab))
          }
        } else xlab
      } else "X"

      default.plot <- modifyList(list(col = col, pch = pch, ylab = ylab,
                                      ylim = c(-ylim0, ylim0 + 1),
                                      main = main.title,
                                      xlab = current_xlab, font.lab = 2), args)
      do.call(plot, c(list(x = xvec, y = Zresidual[, j]), default.plot))
      add_lowess_line(xvec, Zresidual[, j], args)

      plot_limits <- par("usr")
      plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
      if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

      if (!is.null(legend.args)) {
        do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      }
      if (!is.null(test.legend)) {
        do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05,
                               y = plot_limits[4] * 0.7), test.legend))
      }

      if (isTRUE(outlier.return)) {
        if (!identical(id.outlier, integer(0))) {

          default.outlier <- list(pos = 4, labels = id.outlier, cex = 0.8,
                                  col = "darkolivegreen4", add = TRUE, inches = FALSE,
                                  circles = rep((par("usr")[2] - par("usr")[1]) * 0.03,
                                                length(id.outlier)),
                                  fg = "darkolivegreen4", font = 2)
          outlier.args  <- modifyList(default.outlier, outlier.set)
          text.args     <- outlier.args[names(outlier.args) %in% formalArgs(text.default)]
          symbols.args  <- outlier.args[names(outlier.args) %in% formalArgs(symbols)]

          symbols.args$add <- NULL
          effective_symbols_args <- c(list(x = xvec[id.outlier], y = Zresidual[, j][id.outlier], add = TRUE),
                                      symbols.args)
          do.call(symbols, effective_symbols_args)

          text_x <- xvec[id.outlier]
          text_y <- Zresidual[, j][id.outlier]
          text_pos <- outlier.args$pos
          y_median <- median(Zresidual[, j], na.rm = TRUE)
          for (ii in seq_along(id.outlier)) {
            if (!is.na(text_y[ii])) {
              if (text_y[ii] < y_median) text_pos[ii] <- 3 else text_pos[ii] <- 1
            }
          }
          text.args_no_pos <- text.args[names(text.args) != "pos"]
          do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
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

      next
    }

    if (X != "index") {
      test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
      current_test_pv <- NULL # Initialize for each iteration
      for (a in normality.test) {
        test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
        current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
      }

      test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
                          cex = 1, bty = "n", xpd = TRUE, adj = c(0, 0.5))
    }
    # default.legend <- list(legend = unique.cats, col = legend_colors,pch = legend_pchs,
    #                        #  col = unique(col),pch = unique(pch),
    #                        cex = 1, xpd = TRUE, bty = "n",
    #                        title = if (!hasArg("title")) default.legend.title else title,
    #                        horiz = FALSE, y.intersp = 1)
    # legend.args <- modifyList(default.legend, args[!names(args) %in% c("col", "pch")])
    # legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
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
                                    main = main.title,
                                    xlab = current_xlab,font.lab = 2), args)
    if (X == "index") {
      do.call(plot.default, c(list(x = Zresidual[, j]), default.plot))

      plot_limits <- par("usr")
      plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
      if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

      if (!is.null(unique.cats)) {
        do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))

      }

      #do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
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
      add_lowess_line(fitted.value, Zresidual[, j], args)

      plot_limits <- par("usr")
      plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
      if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

      if (!is.null(legend.args)) {
        do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
      }

      if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
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
        add_lowess_line(fitted.value[, i], Zresidual[, j], args)

        plot_limits <- par("usr")
        plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
        if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])

        if (!is.null(legend.args)) {
          do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
        }
        if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
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








# plot.zresid <- function(Zresidual, irep = 1:ncol(Zresidual), ylab = "Z-Residual",
#                         normality.test = c("SW", "AOV", "BL"), k.test = 10,
#                         X = c("index", "covariate", "lp"),
#                         main.title = ifelse(is.null(attr(Zresidual, "type")),
#                                             "Z-residual Scatterplot",
#                                             paste("Z-residual Scatterplot -",
#                                                   attr(Zresidual, "type"))),
#                         outlier.return = TRUE, outlier.value = 3.5,
#                         category = NULL, outlier.set = list(), xlab = NULL,
#                         ...) {
#
#   sign.na <- function(x) {
#     sign.x <- sign(x)
#     sign.x[is.infinite(x)] <- 1
#     sign.x
#   }
#   as.character.na <- function(x) {
#     label.x <- as.character(x)
#     label.x[is.infinite(x)] <- "Inf"
#     label.x
#   }
#   args <- list(...)
#   var.call <- match.call()
#   unique.cats <- NULL
#   default.legend.title <- NULL
#   type <- attr(Zresidual, "type")
#   zero.id <- attr(Zresidual, "zero_id")
#   test.pv <- NULL # Initialize test.pv outside the loop
#
#   ## --- choose colors/pchs & legend defaults (fixed) ---
#   legend_colors <- NULL
#   legend_pchs   <- NULL
#
#   if (!is.null(category)) {
#     # map per-observation by category
#     unique.cats <- unique(category)
#     default.legend.title <- deparse(substitute(category))
#
#     pal  <- if (!is.null(args[["col"]])) args[["col"]] else
#       if (length(unique.cats) == 1) "red" else rainbow(max(length(unique.cats), 2))[seq_along(unique.cats)]
#     pchs <- if (!is.null(args[["pch"]])) args[["pch"]] else
#       if (length(unique.cats) == 1) 1 else (1:25)[seq_along(unique.cats)]
#
#     col <- pal[match(category, unique.cats)]
#     pch <- pchs[match(category, unique.cats)]
#     legend_colors <- pal
#     legend_pchs   <- pchs
#
#   } else if (is.null(type)) {
#     unique.cats <- NULL
#     col <- "black"
#     pch <- 1
#     legend_colors <- NULL
#     legend_pchs   <- NULL
#
#   } else if (type == "survival") {
#     unique.cats <- c("Uncensored", "Censored")
#     censored <- attr(Zresidual, "censored.status")
#     col <- c("blue","red")[censored + 1]
#     pch <- c(3,2)[censored + 1]
#     legend_colors <- c("blue", "red")
#     legend_pchs   <- c(3, 2)
#
#   } else {
#     if (type == "hurdle") {
#       legend_labels <- c("count", "zero")
#       legend_colors <- c("blue", "red")
#       legend_pchs   <- c(3, 2)
#       n <- NROW(Zresidual)
#       is_zero <- seq_len(n) %in% attr(Zresidual, "zero_id")
#       col <- c("blue", "red")[is_zero + 1]
#       pch <- c(3, 2)[is_zero + 1]
#       unique.cats <- legend_labels
#
#     } else if (type %in% c("count")) {
#       unique.cats <- type
#       col <- "blue"; pch <- 3
#       legend_colors <- "blue"; legend_pchs <- 3
#
#     } else if (type %in% c("zero")) {
#       unique.cats <- type
#       col <- "red"; pch <- 2
#       legend_colors <- "red"; legend_pchs <- 2
#
#     } else {
#       col <- "red"; pch <- 2
#       legend_colors <- "red"; legend_pchs <- 2
#     }
#
#     # If user passed col/pch in ..., honor them without altering legend labels
#     if (!is.null(args[["col"]])) col <- args[["col"]]
#     if (!is.null(args[["pch"]])) pch <- args[["pch"]]
#   }
#
#   ## --- build legend args only when needed (fixed) ---
#   legend.args <- NULL
#   if (!is.null(unique.cats)) {
#     default.legend <- list(
#       legend = unique.cats,
#       col = legend_colors,
#       pch = legend_pchs,
#       cex = 0.6, xpd = TRUE, bty = "n",
#       title = if (!hasArg("title")) default.legend.title else title,
#       horiz = FALSE, y.intersp = 1
#     )
#     # Only allow args that legend() knows, but keep title if present
#     user_legend_overrides <- args[!names(args) %in% c("col", "pch")]
#     legend.args <- modifyList(default.legend, user_legend_overrides)
#     legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
#   }
#
#
#   for (j in irep) {
#     #   plot.new()
#     par(mar = c(5, 4, 4, 6) + 0.1)
#     id.nan <- which(is.nan(Zresidual[, j]))
#     id.infinity <- which(is.infinite(Zresidual[, j]))
#     id.outlier <- which(abs(Zresidual[, j]) > outlier.value |
#                           is.infinite(Zresidual[, j]))
#     if (length(id.infinity) > 0L) {
#       value.notfinite <- as.character.na(Zresidual[, j][id.infinity])
#       max.non.infinity <- max(abs(Zresidual[, j]), na.rm = TRUE)
#       Zresidual[, j][id.infinity] <- sign.na(Zresidual[, j][id.infinity]) *
#         (max.non.infinity + 0.1)
#       message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
#     }
#     if (length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
#     ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[, j]), na.rm = TRUE))
#     test.legend <- NULL
#     if (X != "index") {
#       test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
#       current_test_pv <- NULL # Initialize for each iteration
#       for (a in normality.test) {
#         test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
#         current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
#       }
#
#       test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
#                           cex = 0.6, bty = "n", xpd = TRUE, adj = c(0, 0.5))
#     }
#     default.legend <- list(legend = unique.cats, col = legend_colors,pch = legend_pchs,
#                            #  col = unique(col),pch = unique(pch),
#                            cex = 0.6, xpd = TRUE, bty = "n",
#                            title = if (!hasArg("title")) default.legend.title else title,
#                            horiz = FALSE, y.intersp = 1)
#     legend.args <- modifyList(default.legend, args[!names(args) %in% c("col", "pch")])
#     legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
#     default.outlier <- list(pos = 4, labels = id.outlier, cex = 0.8,
#                             col = "darkolivegreen4", add = TRUE, inches = FALSE,
#                             circles = rep((par("usr")[2] - par("usr")[1]) * 0.03,
#                                           length(id.outlier)),
#                             fg = "darkolivegreen4", font = 2)
#     outlier.args <- modifyList(default.outlier, outlier.set)
#     text.args <- outlier.args[names(outlier.args) %in% formalArgs(text.default)]
#     symbols.args <- outlier.args[names(outlier.args) %in% formalArgs(symbols)]
#
#     current_xlab <- if (!is.null(xlab)) {
#       if (startsWith(as.character(xlab), "tex(")) {
#         if (requireNamespace("latex2exp", quietly = TRUE)) {
#           latex_string <- sub("tex\\((.*)\\)", "\\1", as.character(xlab))
#           latex2exp::TeX(latex_string)
#         } else {
#           warning("The 'latex2exp' package is not installed. LaTeX formatting for xlab will not be used.")
#           sub("tex\\((.*)\\)", "\\1", as.character(xlab))
#         }
#       } else {
#         xlab
#       }
#     } else if (X == "index") {
#       "Index"
#     } else if (X == "lp") {
#       fitted.value <- attr(Zresidual, "linear.pred")
#       if (!is.null(fitted.value)) "Linear Predictor" else "Fitted Value"
#     } else if (X != "index" & X != "lp") {
#       fitted.value <- attr(Zresidual, "covariates")
#       if (!is.null(fitted.value)) {
#         cov.name <- variable.names(fitted.value)
#         i <- if (X == "covariate") 1 else which(cov.name == X)
#         colnames(fitted.value)[i]
#       } else {
#         if (X == "covariate") "Covariate" else X
#       }
#     } else {
#       "x"
#     }
#
#     default.plot <- modifyList(list(col = col, pch = pch, ylab = ylab,
#                                     ylim = c(-ylim0, ylim0 + 1),
#                                     main = main.title,
#                                     xlab = current_xlab,font.lab = 2), args)
#     if (X == "index") {
#       do.call(plot.default, c(list(x = Zresidual[, j]), default.plot))
#
#       plot_limits <- par("usr")
#       plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
#       if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])
#
#       if (!is.null(unique.cats)) {
#         do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
#
#       }
#
#       #do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
#       if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
#       if (isTRUE(outlier.return)) {
#         if (!identical(id.outlier, integer(0))) {
#           # Recalculate circles based on current par("usr")
#           outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
#           # Update symbols.args with the recalculated circles
#           symbols.args$circles <- outlier.circles
#           symbols.args$add <- NULL
#           effective_symbols_args <- c(list(x = id.outlier, y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
#           do.call(symbols, effective_symbols_args)
#           # Adjust text position for the last outlier
#           text_x <- id.outlier
#           text_y <- Zresidual[, j][id.outlier]
#           text_pos <- outlier.args$pos # Default position
#
#           y_median <- median(Zresidual[, j], na.rm = TRUE)
#
#           for (i in seq_along(id.outlier)) {
#             if (!is.na(text_y[i])) {
#               if (text_y[i] < y_median) {
#                 text_pos[i] <- 3 # Position above for lower outliers
#               } else {
#                 text_pos[i] <- 1 # Position below for upper outliers
#               }
#             }
#           }
#           text.args_no_pos <- text.args[names(text.args) != "pos"]
#           do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
#           #do.call(text, c(list(id.outlier, Zresidual[, j][id.outlier]), text.args))
#         }
#       }
#     }
#     if (X == "lp") {
#       fitted.value <- attr(Zresidual, "linear.pred")
#       do.call(plot, c(list(x = fitted.value, y = Zresidual[, j]), default.plot))
#
#       plot_limits <- par("usr")
#       plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
#       if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])
#
#       if (!is.null(legend.args)) {
#         do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
#       }
#
#       if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
#       if (isTRUE(outlier.return)) {
#         if (!identical(id.outlier, integer(0))) {
#           outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
#           symbols.args$circles <- outlier.circles
#           symbols.args$add <- NULL
#           effective_symbols_args <- c(list(x = fitted.value[id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
#           do.call(symbols, effective_symbols_args)
#           #do.call(symbols, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), symbols.args))
#
#           text_x <- fitted.value[id.outlier]
#           text_y <- Zresidual[, j][id.outlier]
#           text_pos <- outlier.args$pos # Default position
#           y_median <- median(Zresidual[, j], na.rm = TRUE)
#           for (i in seq_along(id.outlier)) {
#             if (!is.na(text_y[i])) {
#               if (text_y[i] < y_median) {
#                 text_pos[i] <- 3 # Position above for lower outliers
#               } else {
#                 text_pos[i] <- 1 # Position below for upper outliers
#               }
#             }
#           }
#           text.args_no_pos <- text.args[names(text.args) != "pos"]
#           do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
#           #  do.call(text, c(list(fitted.value[id.outlier], Zresidual[, j][id.outlier]), text.args))
#         }
#       }
#     }
#     if (X != "index" & X != "lp") {
#       fitted.value <- attr(Zresidual, "covariates")
#       if (!is.null(fitted.value)) {
#         cov.name <- variable.names(fitted.value)
#         i <- if (X == "covariate") 1 else which(cov.name == X)
#         do.call(plot, c(list(x = fitted.value[, i], y = Zresidual[, j]), default.plot))
#
#         plot_limits <- par("usr")
#         plot_lim_convert <- setNames(c(1:2, 3:4, 1:4), c("x", "x", "y", "y", "xy", "xy", "xy", "xy"))
#         if(!is.null(args[["log"]])) plot_limits[names(plot_lim_convert)==args[["log"]]] <- 10^(plot_limits[names(plot_lim_convert)==args[["log"]]])
#
#         if (!is.null(legend.args)) {
#           do.call(legend, c(list(x = plot_limits[2], y = plot_limits[4]), legend.args))
#         }
#         if (!is.null(test.legend)) do.call(legend, c(list(x = plot_limits[2]+0.5 - (par("usr")[2] - par("usr")[1]) * 0.05, y = plot_limits[4] * 0.7), test.legend))
#         if (isTRUE(outlier.return)) {
#           if (!identical(id.outlier, integer(0))) {
#             # Recalculate circles based on current par("usr")
#             outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
#             # Update symbols.args with the recalculated circles
#             symbols.args$circles <- outlier.circles
#             symbols.args$add <- NULL
#             effective_symbols_args <- c(list(x = fitted.value[, i][id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
#             do.call(symbols, effective_symbols_args)
#             #do.call(symbols, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), symbols.args))
#             text_x <- fitted.value[, i][id.outlier]
#             text_y <- Zresidual[, j][id.outlier]
#             text_pos <- outlier.args$pos
#             y_median <- median(Zresidual[, j], na.rm = TRUE)
#             for (i in seq_along(id.outlier)) {
#               if (!is.na(text_y[i])) {
#                 if (text_y[i] < y_median) {
#                   text_pos[i] <- 3 # Position above for lower outliers
#                 } else {
#                   text_pos[i] <- 1 # Position below for upper outliers
#                 }
#               }
#             }
#             text.args_no_pos <- text.args[names(text.args) != "pos"]
#             do.call(text, c(list(x = text_x, y = text_y, pos = text_pos), text.args_no_pos))
#             # do.call(text, c(list(fitted.value[, i][id.outlier], Zresidual[, j][id.outlier]), text.args))
#           }
#         }
#       }
#     }
#
#     if (length(id.infinity) > 0L) {
#       text(id.infinity + 5, Zresidual[, j][id.infinity],
#            labels = value.notfinite, col = 2)
#     }
#     hlines <- c(1.96, 3)
#     hlines2 <- -hlines
#     abline(h = c(hlines, hlines2), lty = 3, col = "grey")
#     if (outlier.return & !identical(id.outlier, integer(0))) {
#       cat("Outlier Indices(", attr(Zresidual, "type"), "):", id.outlier, "\n",
#           sep = " ")
#       invisible(list(outliers = id.outlier))
#     }
#   }
#   par(mar = c(5, 4, 4, 2) + 0.1)
# }
#
