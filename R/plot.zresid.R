#' Plot Z-Residuals for Bayesian and Frequentist Count / Hurdle / Zero / Survival Models
#'
#' @description
#' Produces diagnostic scatterplots of Z-residuals for a wide range of model
#' types, including hurdle models, zero models, count models, and survival models.
#'
#' This method supports plotting Z-residuals against index, covariates, linear predictors,
#' or any user-specified x-axis vector (e.g., time).
#'
#' @usage
#' \method{plot}{zresid}(
#'   x,
#'   info = NULL,
#'   irep = 1:ncol(x),
#'   ylab = "Z-Residual",
#'   normality.test = c("SW", "AOV", "BL"),
#'   k.test = 10,
#'   x_axis_var = "index",
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
#' A numeric matrix of Z-residuals (dimensions: \eqn{n \times m}) of class `"zresid"`.
#' The method may use attributes such as `"type"`, `"zero_id"`, `"censored.status"`,
#' `"covariates"`, and `"linear.pred"` if present.
#'
#' @param info
#' Optional metadata (typically the output of `Zcov()` or an equivalent list).
#' When provided, this method can fill legacy attributes (e.g., `"covariates"`, `"linear.pred"`,
#' `"zero_id"`, `"censored.status"`, `"type"`) without changing the plotting behavior.
#'
#' @param irep
#' A vector specifying which residual column(s) to plot. Defaults to all columns.
#'
#' @param ylab Label for the y-axis.
#'
#' @param normality.test
#' Character vector specifying which normality tests to compute per iteration:
#' `"SW"` (Shapiro-Wilk), `"AOV"` (ANOVA-based), `"BL"` (Bartlett-based).
#' Note: if `x_axis_var` is non-numeric (e.g., `Date`, `POSIXct`, factor), some tests
#' may return `NA` or be skipped depending on helper implementations.
#'
#' @param k.test
#' Bin size used by normality tests when grouping a continuous x-axis predictor.
#'
#' @param x_axis_var
#' Specifies x-axis values. Supported inputs:
#' \itemize{
#'   \item A keyword: `"index"`, `"lp"`, or `"covariate"`.
#'   \item A covariate name present in `attr(x, "covariates")`.
#'   \item A length-\eqn{n} vector (numeric / `Date` / `POSIXct` / factor / character), e.g., time.
#'   \item A function `function(z, info) ...` returning a length-\eqn{n} vector; optionally it can return
#'         `list(values = <vector>, label = <character>)` to set an automatic x-axis label.
#' }
#'
#' @param main.title Plot title. Defaults to a type-based informative title if available.
#'
#' @param outlier.return
#' If `TRUE`, prints outlier indices to console and returns them invisibly.
#'
#' @param outlier.value
#' Threshold above which a residual is flagged as an outlier. Defaults to `3.5`.
#'
#' @param category
#' Optional vector categorizing observations (length \eqn{n}), used for coloring/shaping points.
#'
#' @param outlier.set
#' A named list of arguments passed to `symbols()` / `text()` for marking and labeling outliers.
#'
#' @param xlab
#' Label for the x-axis. If given as `tex("...")`, it will be interpreted via `latex2exp` if installed.
#'
#' @param my.mar
#' Plot margins passed to `par(mar = ...)`.
#'
#' @param add_lowess
#' Logical. If `TRUE`, add a LOWESS smooth (only applicable when x is numeric / `Date` / `POSIXct`;
#' otherwise the smooth may be skipped).
#'
#' @param ...
#' Additional graphical arguments passed to `plot()`, `legend()`, `symbols()`, and `text()`.
#'
#' @details
#' S3 method for objects of class `"zresid"`.
#'
#' Model-type behavior:
#' \itemize{
#'   \item Hurdle models: zeros and counts can be colored differently.
#'   \item Survival models: censored and uncensored observations can be separated.
#' }
#'
#' Outliers are defined as \eqn{|Z| > outlier.value} or non-finite values.
#'
#' @return
#' Invisibly returns (when `outlier.return = TRUE`) a list with:
#' \describe{
#'   \item{outliers}{Integer vector of outlier indices.}
#' }
#' Otherwise returns `NULL`. Always produces a scatterplot.
#'
#' @note
#' This function temporarily modifies graphical parameters and restores them on exit.
#'
#' @exportS3Method graphics::plot zresid
#' @export plot.zresid
plot.zresid <- function(x, info = NULL, irep = 1:ncol(x), ylab = "Z-Residual",
                        normality.test = c("SW", "AOV", "BL"), k.test = 10,
                        x_axis_var = "index",
                        main.title = ifelse(is.null(attr(x, "type")),
                                            "Z-residual Scatterplot",
                                            paste("Z-residual Scatterplot -",
                                                  attr(x, "type"))),
                        outlier.return = TRUE, outlier.value = 3.5,
                        category = NULL, outlier.set = list(), xlab = NULL,
                        my.mar=c(5,4,4,6)+0.1, add_lowess = FALSE, ...) {
  
  Zresidual <- x
  x_axis_expr <- substitute(x_axis_var)
  xlab_from_xaxis <- NULL
  
  # ---- Compatibility layer: support "new" format where metadata is separated (e.g., from Zcov()) ----
  info0 <- info
  if (is.null(info0)) {
    info0 <- attr(Zresidual, "info")
    if (is.null(info0)) info0 <- attr(Zresidual, "zcov")
    if (is.null(info0)) info0 <- attr(Zresidual, "Zcov")
  }
  
  .fill_attr_if_missing <- function(obj, nm, value) {
    if (is.null(attr(obj, nm)) && !is.null(value)) attr(obj, nm) <- value
    obj
  }
  
  .map_type_from_info <- function(info) {
    if (!is.null(info$type)) return(info$type)
    kind <- info$y_type_kind
    if (is.null(kind)) return(NULL)
    if (identical(kind, "censor")) return("survival")
    if (identical(kind, "hurdle")) return("hurdle")
    if (identical(kind, "trunc"))  return("count")
    if (identical(kind, "plain"))  return(NULL)
    NULL
  }
  
  if (!is.null(info0)) {
    cov0 <- info0$covariates
    if (is.null(cov0)) cov0 <- info0$covariate
    lp0 <- info0$linear_pred
    if (is.null(lp0)) lp0 <- info0$linear.pred
    Zresidual <- .fill_attr_if_missing(Zresidual, "covariates", cov0)
    Zresidual <- .fill_attr_if_missing(Zresidual, "linear.pred", lp0)
    
    if (is.null(attr(Zresidual, "type"))) {
      attr(Zresidual, "type") <- .map_type_from_info(info0)
    }
    
    if (is.null(attr(Zresidual, "zero_id"))) {
      zid <- NULL
      if (!is.null(info0$extra) && !is.null(info0$extra$zero_id)) zid <- info0$extra$zero_id
      if (is.null(zid) && !is.null(info0$zero_id)) zid <- info0$zero_id
      if (!is.null(zid)) attr(Zresidual, "zero_id") <- zid
    }
    
    if (is.null(attr(Zresidual, "censored.status"))) {
      if (!is.null(info0$censored.status)) {
        attr(Zresidual, "censored.status") <- info0$censored.status
      } else if (identical(info0$y_type_kind, "censor")) {
        yt <- info0$y_type
        if (is.null(yt)) yt <- attr(Zresidual, "y_type")
        if (!is.null(yt)) attr(Zresidual, "censored.status") <- as.integer(yt == 0L)
      }
    }
  }
  
  # If caller didn't set main.title, recompute after we possibly filled attr(type)
  if (missing(main.title)) {
    main.title <- ifelse(
      is.null(attr(Zresidual, "type")),
      "Z-residual Scatterplot",
      paste("Z-residual Scatterplot -", attr(Zresidual, "type"))
    )
  }
  
  # ---- x_axis_var: allow keyword / covariate name / length-n vector / function ----
  X <- x_axis_var
  if (is.character(X) && length(X) > 1L) X <- X[1L]  # safety for old default style
  
  # If user passes a function, evaluate it to get x vector.
  # It can return a vector length n, or list(values=..., label=...)
  if (is.function(X)) {
    xr <- X(Zresidual, info0)
    if (is.list(xr) && !is.null(xr$values)) {
      xlab_from_xaxis <- xr$label
      X <- xr$values
    } else {
      X <- xr
    }
  }
  
  op_mar <- graphics::par("mar")
  op_xpd <- graphics::par("xpd")
  on.exit(graphics::par(mar = op_mar, xpd = op_xpd), add = TRUE)
  
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
  
  # Convert Date/POSIXct to numeric for lowess checks (plot itself can keep Date/POSIXct)
  .as_numeric_for_smooth <- function(v) {
    if (inherits(v, "POSIXt") || inherits(v, "Date")) return(as.numeric(v))
    v
  }
  
  add_lowess_line <- function(xv, yv, args) {
    if (!isTRUE(add_lowess)) return(invisible(NULL))
    
    xv_num <- .as_numeric_for_smooth(xv)
    if (!is.numeric(xv_num)) return(invisible(NULL))
    
    log_opt <- args[["log"]]
    log_x <- !is.null(log_opt) && grepl("x", log_opt)
    
    ok <- is.finite(xv_num) & is.finite(yv)
    if (log_x) ok <- ok & (xv_num > 0)
    
    if (sum(ok) < 3L) return(invisible(NULL))
    
    if (log_x) {
      lx <- log10(xv_num[ok])
      lw <- stats::lowess(x = lx, y = yv[ok])
      graphics::lines(x = 10^lw$x, y = lw$y, col = "red", lwd = 3)
    } else {
      lw <- stats::lowess(x = xv_num[ok], y = yv[ok])
      # If original x is Date/POSIXct, plot uses numeric internally, so numeric lines will align
      graphics::lines(lw$x, lw$y, col = "red", lwd = 3)
    }
    invisible(NULL)
  }
  
  draw_legends_outside_right <- function(legend.args, test.legend,
                                         mar_right_min = 8, gap_factor = 0.8) {
    if (is.null(legend.args) && is.null(test.legend)) return(invisible(NULL))
    
    mar <- graphics::par("mar")
    if (mar[4] < mar_right_min) graphics::par(mar = c(mar[1], mar[2], mar[3], mar_right_min))
    
    old_xpd <- graphics::par("xpd")
    on.exit(graphics::par(xpd = old_xpd), add = TRUE)
    graphics::par(xpd = TRUE)
    
    usr <- graphics::par("usr")
    x_anchor <- usr[2] + (usr[2] - usr[1]) * 0.02
    y_top    <- usr[4]
    
    lg1 <- NULL
    if (!is.null(legend.args)) {
      la <- legend.args
      la$plot  <- NULL
      la$inset <- NULL
      la$x <- NULL
      la$y <- NULL
      la$xpd <- TRUE
      
      lg1 <- do.call(graphics::legend, c(
        list(x = x_anchor, y = y_top, xjust = 0, yjust = 1, plot = FALSE),
        la
      ))
      
      do.call(graphics::legend, c(
        list(x = x_anchor, y = y_top, xjust = 0, yjust = 1),
        la
      ))
    }
    
    if (!is.null(test.legend)) {
      tl <- test.legend
      tl$plot  <- NULL
      tl$inset <- NULL
      tl$x <- NULL
      tl$y <- NULL
      tl$xpd <- TRUE
      tl$bg <- "white"
      tl$box.col <- NA
      
      x2 <- x_anchor
      y2 <- y_top
      
      if (!is.null(lg1)) {
        gap <- graphics::strheight("M", units = "user") * gap_factor
        x2  <- lg1$rect$left
        y2  <- lg1$rect$top - lg1$rect$h - gap
      }
      
      do.call(graphics::legend, c(
        list(x = x2, y = y2, xjust = 0, yjust = 1),
        tl
      ))
    }
    
    invisible(NULL)
  }
  
  #####################
  unique.cats <- NULL
  default.legend.title <- NULL
  type <- attr(Zresidual, "type")
  zero.id <- attr(Zresidual, "zero_id")
  
  legend_colors <- NULL
  legend_pchs   <- NULL
  
  if (!is.null(category)) {
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
    
    if (!is.null(args[["col"]])) col <- args[["col"]]
    if (!is.null(args[["pch"]])) pch <- args[["pch"]]
  }
  
  legend.args <- NULL
  if (!is.null(unique.cats) && length(unique.cats) > 0 &&
      !is.null(legend_colors) && length(legend_colors) > 0 &&
      !is.null(legend_pchs)   && length(legend_pchs)   > 0) {
    default.legend <- list(
      legend = unique.cats,
      col = legend_colors,
      pch = legend_pchs,
      cex = 0.9, xpd = TRUE, bty = "n",
      title = if (!hasArg("title")) default.legend.title else title,
      horiz = FALSE, y.intersp = 1
    )
    user_legend_overrides <- args[!names(args) %in% c("col", "pch")]
    legend.args <- modifyList(default.legend, user_legend_overrides)
    legend.args <- legend.args[names(legend.args) %in% formalArgs(legend)]
  }
  
  # --- validate x_axis_var ---
  choices <- c("index", "covariate", "lp")
  n_obs   <- NROW(Zresidual)
  
  is_uservec <- (length(X) == n_obs)
  is_keyword <- is.character(X) && length(X) == 1 && X %in% choices
  is_covname <- is.character(X) && length(X) == 1 && !is_keyword
  
  if (!is_uservec && !is_keyword && !is_covname) {
    stop("x_axis_var must be one of: ",
         "'index','lp','covariate', a covariate name in attr(z,'covariates'), ",
         "a length-n vector (e.g., time), or a function(z, info) returning a length-n vector.")
  }
  
  # Helper: xlab rendering (supports tex("..."))
  .resolve_xlab <- function(lbl) {
    if (is.null(lbl)) return(NULL)
    if (startsWith(as.character(lbl), "tex(")) {
      if (requireNamespace("latex2exp", quietly = TRUE)) {
        latex_string <- sub("tex\\((.*)\\)", "\\1", as.character(lbl))
        return(latex2exp::TeX(latex_string))
      } else {
        warning("The 'latex2exp' package is not installed. LaTeX formatting for xlab will not be used.")
        return(sub("tex\\((.*)\\)", "\\1", as.character(lbl)))
      }
    }
    lbl
  }
  
  for (j in irep) {
    par(mar = my.mar)
    
    id.nan <- which(is.nan(Zresidual[, j]))
    id.infinity <- which(is.infinite(Zresidual[, j]))
    id.outlier <- which(abs(Zresidual[, j]) > outlier.value | is.infinite(Zresidual[, j]))
    
    if (length(id.infinity) > 0L) {
      value.notfinite <- as.character.na(Zresidual[, j][id.infinity])
      max.non.infinity <- max(abs(Zresidual[, j]), na.rm = TRUE)
      Zresidual[, j][id.infinity] <- sign.na(Zresidual[, j][id.infinity]) * (max.non.infinity + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
    
    ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[, j]), na.rm = TRUE))
    test.legend <- NULL
    
    # --------------------------
    # USER VECTOR MODE (time/any)
    # --------------------------
    if (is_uservec && !(is.character(X) && length(X) == 1 && X %in% choices)) {
      xvec <- X
      
      # tests (may fail for non-numeric; will be NA)
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
      
      # xlab default: user-provided > function label > expression label > "X"
      xlab_auto <- if (!is.null(xlab)) xlab else if (!is.null(xlab_from_xaxis)) xlab_from_xaxis else {
        # best-effort label from expression
        lab <- paste(deparse(x_axis_expr), collapse = "")
        if (identical(lab, "x_axis_var")) "X" else lab
      }
      current_xlab <- .resolve_xlab(xlab_auto)
      
      default.plot <- modifyList(list(col = col, pch = pch, ylab = ylab,
                                      ylim = c(-ylim0, ylim0 + 1),
                                      main = main.title,
                                      xlab = current_xlab, font.lab = 2), args)
      do.call(plot, c(list(x = xvec, y = Zresidual[, j]), default.plot))
      add_lowess_line(xvec, Zresidual[, j], args)
      
      draw_legends_outside_right(legend.args, test.legend)
      
      if (isTRUE(outlier.return) && !identical(id.outlier, integer(0))) {
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
      
      if (length(id.infinity) > 0L) {
        text(id.infinity + 5, Zresidual[, j][id.infinity],
             labels = value.notfinite, col = 2)
      }
      hlines <- c(1.96, 3)
      abline(h = c(hlines, -hlines), lty = 3, col = "grey")
      
      if (outlier.return && !identical(id.outlier, integer(0))) {
        cat("Outlier Indices(", attr(Zresidual, "type"), "):", id.outlier, "\n", sep = " ")
        invisible(list(outliers = id.outlier))
      }
      next
    }
    
    # --------------------------
    # KEYWORD / COVARIATE MODE
    # --------------------------
    Xk <- X
    if (is.character(Xk) && length(Xk) > 1L) Xk <- Xk[1L]
    
    if (Xk != "index") {
      test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
      current_test_pv <- NULL
      for (a in normality.test) {
        test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, Xk, k.test)
        current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
      }
      test.legend <- list(legend = c(expression(bold("P-value:")), current_test_pv),
                          cex = 1, bty = "n", xpd = TRUE, adj = c(0, 0.5))
    }
    
    default.outlier <- list(pos = 4, labels = id.outlier, cex = 0.8,
                            col = "darkolivegreen4", add = TRUE, inches = FALSE,
                            circles = rep((par("usr")[2] - par("usr")[1]) * 0.03,
                                          length(id.outlier)),
                            fg = "darkolivegreen4", font = 2)
    outlier.args <- modifyList(default.outlier, outlier.set)
    text.args <- outlier.args[names(outlier.args) %in% formalArgs(text.default)]
    symbols.args <- outlier.args[names(outlier.args) %in% formalArgs(symbols)]
    
    current_xlab <- if (!is.null(xlab)) {
      .resolve_xlab(xlab)
    } else if (Xk == "index") {
      "Index"
    } else if (Xk == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      if (!is.null(fitted.value)) "Linear Predictor" else "Fitted Value"
    } else {
      fitted.value <- attr(Zresidual, "covariates")
      if (!is.null(fitted.value)) {
        cov.name <- variable.names(fitted.value)
        i <- if (Xk == "covariate") 1 else which(cov.name == Xk)
        colnames(fitted.value)[i]
      } else {
        if (Xk == "covariate") "Covariate" else Xk
      }
    }
    
    default.plot <- modifyList(list(col = col, pch = pch, ylab = ylab,
                                    ylim = c(-ylim0, ylim0 + 1),
                                    main = main.title,
                                    xlab = current_xlab, font.lab = 2), args)
    
    if (Xk == "index") {
      do.call(plot.default, c(list(x = Zresidual[, j]), default.plot))
      draw_legends_outside_right(legend.args, test.legend)
      
      if (isTRUE(outlier.return) && !identical(id.outlier, integer(0))) {
        outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
        symbols.args$circles <- outlier.circles
        symbols.args$add <- NULL
        effective_symbols_args <- c(list(x = id.outlier, y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
        do.call(symbols, effective_symbols_args)
        
        text_x <- id.outlier
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
    
    if (Xk == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      do.call(plot, c(list(x = fitted.value, y = Zresidual[, j]), default.plot))
      add_lowess_line(fitted.value, Zresidual[, j], args)
      draw_legends_outside_right(legend.args, test.legend)
      
      if (isTRUE(outlier.return) && !identical(id.outlier, integer(0))) {
        outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
        symbols.args$circles <- outlier.circles
        symbols.args$add <- NULL
        effective_symbols_args <- c(list(x = fitted.value[id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
        do.call(symbols, effective_symbols_args)
        
        text_x <- fitted.value[id.outlier]
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
    
    if (Xk != "index" && Xk != "lp") {
      fitted.value <- attr(Zresidual, "covariates")
      if (!is.null(fitted.value)) {
        cov.name <- variable.names(fitted.value)
        i <- if (Xk == "covariate") 1 else which(cov.name == Xk)
        do.call(plot, c(list(x = fitted.value[, i], y = Zresidual[, j]), default.plot))
        add_lowess_line(fitted.value[, i], Zresidual[, j], args)
        draw_legends_outside_right(legend.args, test.legend)
        
        if (isTRUE(outlier.return) && !identical(id.outlier, integer(0))) {
          outlier.circles <- rep((par("usr")[2] - par("usr")[1]) * 0.03, length(id.outlier))
          symbols.args$circles <- outlier.circles
          symbols.args$add <- NULL
          effective_symbols_args <- c(list(x = fitted.value[, i][id.outlier], y = Zresidual[, j][id.outlier], add = TRUE), symbols.args)
          do.call(symbols, effective_symbols_args)
          
          text_x <- fitted.value[, i][id.outlier]
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
    }
    
    if (length(id.infinity) > 0L) {
      text(id.infinity + 5, Zresidual[, j][id.infinity],
           labels = value.notfinite, col = 2)
    }
    hlines <- c(1.96, 3)
    abline(h = c(hlines, -hlines), lty = 3, col = "grey")
    
    if (outlier.return && !identical(id.outlier, integer(0))) {
      cat("Outlier Indices(", attr(Zresidual, "type"), "):", id.outlier, "\n", sep = " ")
      invisible(list(outliers = id.outlier))
    }
  }
}
