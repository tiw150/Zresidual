#' Scatterplot diagnostics for Z-residuals
#'
#' Produces scatterplots of Z-residuals against observation index, linear
#' predictors, model covariates, or a user-supplied x-axis variable. Optional
#' metadata supplied through \code{info} (or stored in attributes of
#' \code{x}) are used to fill legacy plotting attributes such as
#' \code{"covariates"}, \code{"linear.pred"}, \code{"zero_id"}, and
#' \code{"censored.status"}.
#'
#' Depending on the metadata available, the plot can distinguish zero versus
#' positive observations for hurdle-type models, or censored versus uncensored
#' observations for survival models.
#'
#' @param x A numeric matrix of Z-residuals, typically returned by
#'   \code{\link{Zresidual}}, with one column per residual replicate.
#' @param zcov Optional metadata, typically returned by \code{\link{Zcov}}.
#' @param info Legacy alias for \code{zcov}.
#' @param irep Integer vector specifying which column(s) of \code{x} to plot.
#' @param ylab Label for the y-axis.
#' @param normality.test Character vector specifying which diagnostic p-values
#'   to display. Supported values are \code{"SW"}, \code{"AOV"}, and
#'   \code{"BL"}.
#' @param k.test Integer controlling grouping used by the diagnostic tests.
#' @param x_axis_var Variable used on the x-axis. It may be one of
#'   \code{"index"}, \code{"lp"}, \code{"covariate"}, a covariate name stored
#'   in \code{attr(x, "covariates")}, a length-\eqn{n} vector, or a function
#'   returning such a vector.
#' @param main.title Main title of the plot. If omitted, a default title is
#'   constructed from \code{attr(x, "type")}, when available.
#' @param outlier.return Logical; if \code{TRUE}, mark observations with
#'   \code{|Z| > outlier.value} (and non-finite residuals) and invisibly return
#'   their indices.
#' @param outlier.value Numeric threshold used to define outliers.
#' @param category Optional grouping variable of length \eqn{n} used to modify
#'   point appearance.
#' @param outlier.set A named list of graphical arguments passed to
#'   \code{\link[graphics]{symbols}} and \code{\link[graphics]{text}} when
#'   annotating outliers.
#' @param xlab Label for the x-axis. If \code{NULL}, an automatic label is used.
#' @param my.mar Numeric vector passed to \code{\link[graphics]{par}(mar = ...)}.
#' @param add_lowess Logical; if \code{TRUE}, add a LOWESS smooth when the x-axis
#'   is numeric.
#' @param ... Additional graphical arguments passed to plotting functions.
#'
#' @return Invisibly returns a list with component \code{outliers}, containing
#' the indices of observations flagged as outliers for the plotted replicate.
#' The main effect of the function is the diagnostic scatterplot.
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 30
#'   x <- rnorm(n)
#'   t_event <- rexp(n, rate = exp(0.3 * x))
#'   t_cens  <- rexp(n, rate = 0.4)
#'   status  <- as.integer(t_event <= t_cens)
#'   time    <- pmin(t_event, t_cens)
#'   dat <- data.frame(time = time, status = status, x = x)
#'
#'   fit <- survival::coxph(survival::Surv(time, status) ~ x, data = dat)
#'   z <- Zresidual(fit, data=dat, nrep = 1, seed = 1)
#'   info <- Zcov(fit, data = dat)
#'
#'   plot(z, info = info, x_axis_var = "index")
#'   plot(z, info = info, x_axis_var = "lp")
#' }
#'
#' @seealso \code{\link{Zresidual}}, \code{\link{Zcov}}
#'
#' @method plot zresid
#' @export
plot.zresid <- function(x, zcov = NULL, info = NULL, irep = 1, ylab = "Z-Residual",
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
  info0 <- if (!is.null(zcov)) zcov else info
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
    # Prefer explicit component type from Zcov()
    if (!is.null(info$type) && length(info$type) >= 1L && nzchar(as.character(info$type)[1])) {
      return(as.character(info$type)[1])
    }
    
    kind <- info$y_type_kind
    if (is.null(kind) || length(kind) == 0L) return(NULL)
    kind <- as.character(kind[1])
    
    if (identical(kind, "censor")) return("survival")
    if (identical(kind, "trunc"))  return("count")
    if (identical(kind, "hurdle")) return("hurdle")
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
    xlog <- graphics::par("xlog") 
    
    if (xlog) {
      x_anchor <- 10^(usr[2] + (usr[2] - usr[1]) * 0.02)
    } else {
      x_anchor <- usr[2] + (usr[2] - usr[1]) * 0.02
    }
    y_top <- usr[4]
    
    lg1 <- NULL
    if (!is.null(legend.args)) {
      la <- legend.args
      la$plot <- FALSE
      la$x <- x_anchor
      la$y <- y_top
      la$xjust <- 0
      la$yjust <- 1
      
      lg1 <- do.call(graphics::legend, la)
      
      la$plot <- TRUE
      do.call(graphics::legend, la)
    }
    
    if (!is.null(test.legend)) {
      tl <- test.legend
      tl$xpd <- TRUE
      tl$bg <- "white"
      tl$box.col <- NA
      
      if (!is.null(lg1)) {
        gap <- graphics::strheight("M", units = "user") * gap_factor
        x2 <- x_anchor 
        y2 <- lg1$rect$top - lg1$rect$h - gap
      } else {
        x2 <- x_anchor
        y2 <- y_top
      }
      
      do.call(graphics::legend, c(
        list(x = x2, y = y2, xjust = 0, yjust = 1),
        tl
      ))
    }
    invisible(NULL)
  }
  
  #####################
  #####################
  unique.cats <- NULL
  default.legend.title <- NULL
  type <- attr(Zresidual, "type")
  
  # FIX: normalize type to a single string or NULL (avoid length-0 errors)
  if (is.null(type) || length(type) == 0L || !nzchar(type[1])) type <- NULL else type <- as.character(type[1])
  
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
    
  } else if (identical(type, "survival")) {
    unique.cats <- c("Uncensored", "Censored")
    censored <- attr(Zresidual, "censored.status")
    col <- c("blue","red")[censored + 1]
    pch <- c(3,2)[censored + 1]
    legend_colors <- c("blue", "red")
    legend_pchs   <- c(3, 2)
    
  } else {
    if (identical(type, "hurdle")) {
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
      
    } else if (type %in% c("zero", "logistic")) {
      
      n <- NROW(Zresidual)
      
      # y_type must come from `info` (preferred). fallback only for old objects.
      yt <- NULL
      fam_info <- NULL
      if (!is.null(info) && is.list(info)) {
        if (!is.null(info$y_type) && length(info$y_type) == n) yt <- info$y_type
        if (!is.null(info$family)) fam_info <- as.character(info$family)[1]
      }
      if (is.null(yt)) {
        yt <- attr(Zresidual, "y_type")
        if (!is.null(yt) && length(yt) != n) yt <- NULL
      }
      if (is.null(yt)) {
        stop("plot.zresid (zero/logistic): cannot find y_type. Pass `info = Zcov(...)` aligned with this Zresidual.",
             call. = FALSE)
      }
      
      is_hurdle <- !is.null(fam_info) && grepl("^hurdle(_|$)", fam_info)
      
      if (is_hurdle) {
        # hurdle: yt==0 => zero, yt==1 => count
        legend_labels <- c("count", "zero")
        legend_colors <- c("blue", "red")
        legend_pchs   <- c(3, 2)
        codes <- c(1L, 0L)  # order must match legend_labels above
        
      } else {
        # bernoulli: be faithful to binary code 0/1
        legend_labels <- c("0", "1")
        legend_colors <- c("blue", "red")
        legend_pchs   <- c(3, 2)
        codes <- c(0L, 1L)
      }
      
      idx <- match(as.integer(yt), codes)
      if (anyNA(idx)) {
        stop("plot.zresid (zero/logistic): y_type contains values outside expected codes.", call. = FALSE)
      }
      
      col <- legend_colors[idx]
      pch <- legend_pchs[idx]
      unique.cats <- legend_labels
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
        message("Outlier Indices : ", paste(id.outlier, collapse = ", "))
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
      message("Outlier Indices : ", paste(id.outlier, collapse = ", "))
      invisible(list(outliers = id.outlier))
    }
  }
}
