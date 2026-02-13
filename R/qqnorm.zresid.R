#' Normal Q-Q Plot for Z-Residuals with Outlier Detection and Normality Diagnostics
#'
#' Produces a normal Q-Q plot for Z-residuals, with optional Shapiro–Wilk normality testing,
#' automatic handling of non-finite values, optional axis compression for very large residuals,
#' and visual annotation of detected outliers (indices refer to the original observation order).
#'
#' This diagnostic is designed for checking the normality assumption of Z-residuals obtained
#' from Bayesian predictive model checks (posterior, LOOCV, ISCV, etc.).
#'
#' @usage
#' \method{qqnorm}{zresid}(
#'   y,
#'   info = NULL,
#'   irep = 1,
#'   diagnosis.test = "SW",
#'   main.title = ifelse(is.null(attr(y, "type")),
#'                       "Normal Q-Q Plot",
#'                       paste("Normal Q-Q Plot -", attr(y, "type"))),
#'   xlab = "Theoretical Quantiles",
#'   ylab = "Sample Quantiles",
#'   outlier.return = TRUE,
#'   outlier.value = 3.5,
#'   outlier.set = list(),
#'   my.mar = c(5, 4, 4, 6) + 0.1,
#'   legend.settings = list(),
#'   clip.extreme = TRUE,
#'   clip.threshold = 6,
#'   ...
#' )
#'
#' @param y A numeric matrix of Z-residuals where each column corresponds to an iteration/draw.
#'   The function may use attribute \code{"type"} to build default titles.
#'
#' @param info Optional metadata (typically output of \code{Zcov()} or an equivalent list).
#'   When provided (or attached via \code{attr(y,"info")}/\code{attr(y,"zcov")}/\code{attr(y,"Zcov")}),
#'   this method can fill legacy attributes (e.g., \code{"type"}) for compatibility.
#'
#' @param irep Integer or integer vector selecting which column(s) of \code{y} to plot.
#'
#' @param diagnosis.test Character string indicating the normality test to perform.
#'   Currently supports \code{"SW"} (Shapiro–Wilk). If helper \code{sw.test.zresid()} exists,
#'   it will be used; otherwise a fallback \code{shapiro.test()} is applied per column.
#'
#' @param main.title Main title of the plot. If missing, it is automatically constructed using
#'   \code{attr(y, "type")}.
#'
#' @param xlab,ylab Axis labels for the Q-Q plot.
#'
#' @param outlier.return Logical; if \code{TRUE}, the function prints and returns outlier indices.
#'
#' @param outlier.value Numeric threshold used to classify an observation as an outlier based on
#'   \code{|Z| > outlier.value}. Default is \code{3.5}.
#'
#' @param outlier.set Optional named list of graphical parameters passed to \code{symbols()} and
#'   \code{text()} for customizing outlier circles/labels.
#'
#' @param my.mar Numeric vector passed to \code{par(mar = ...)}.
#'
#' @param legend.settings Optional named list to override default legend appearance.
#'
#' @param clip.extreme Logical. If \code{TRUE}, very large residuals are visually compressed
#'   onto the plot boundary (with optional axis break marks if available).
#'
#' @param clip.threshold Numeric. Threshold for extreme-value clipping (default 6).
#'
#' @param ... Additional graphical arguments passed to \code{plot()} (for Q-Q points).
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{outliers}{Integer vector of detected outlier indices (original observation indices).}
#' }
#'
#' @method qqnorm zresid
#' @export qqnorm.zresid
qqnorm.zresid <- function(y, info = NULL, irep = 1, diagnosis.test = "SW",
                          main.title = ifelse(is.null(attr(y, "type")),
                                              "Normal Q-Q Plot",
                                              paste("Normal Q-Q Plot -", attr(y, "type"))),
                          xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
                          outlier.return = TRUE, outlier.value = 3.5, outlier.set = list(),
                          my.mar = c(5, 4, 4, 6) + 0.1, legend.settings = list(),
                          clip.extreme = TRUE, clip.threshold = 6, ...) {
  
  Zresidual <- y
  
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
    if (is.null(attr(Zresidual, "type"))) {
      attr(Zresidual, "type") <- .map_type_from_info(info0)
    }
  }
  
  if (missing(main.title)) {
    main.title <- ifelse(is.null(attr(Zresidual, "type")),
                         "Normal Q-Q Plot",
                         paste("Normal Q-Q Plot -", attr(Zresidual, "type")))
  }
  
  # ---- helpers ----
  .safe_sw_pvalue <- function(zcol) {
    zf <- zcol[is.finite(zcol)]
    if (length(zf) < 3) return(NA_real_)
    out <- try(stats::shapiro.test(zf)$p.value, silent = TRUE)
    if (inherits(out, "try-error")) NA_real_ else out
  }
  
  .get_sw_vec <- function(Z) {
    # Prefer package helper if present
    if (exists("sw.test.zresid", mode = "function")) {
      tt <- try(sw.test.zresid(Z), silent = TRUE)
      if (!inherits(tt, "try-error")) return(tt)
    }
    apply(Z, 2, .safe_sw_pvalue)
  }
  
  .draw_legends_outside_right <- function(leg1, leg2, gap_factor = 0.8) {
    old_xpd <- graphics::par("xpd")
    on.exit(graphics::par(xpd = old_xpd), add = TRUE)
    graphics::par(xpd = TRUE)
    
    usr <- graphics::par("usr")
    x_anchor <- usr[2] + (usr[2] - usr[1]) * 0.02
    y_top <- usr[4]
    
    lg1 <- NULL
    if (!is.null(leg1)) {
      lg1 <- do.call(graphics::legend, c(list(x = x_anchor, y = y_top, xjust = 0, yjust = 1, plot = FALSE), leg1))
      do.call(graphics::legend, c(list(x = x_anchor, y = y_top, xjust = 0, yjust = 1), leg1))
    }
    
    if (!is.null(leg2)) {
      x2 <- x_anchor
      y2 <- y_top
      if (!is.null(lg1)) {
        gap <- graphics::strheight("M", units = "user") * gap_factor
        x2 <- lg1$rect$left
        y2 <- lg1$rect$top - lg1$rect$h - gap
      }
      do.call(graphics::legend, c(list(x = x2, y = y2, xjust = 0, yjust = 1), leg2))
    }
  }
  
  # Optional axis break (plotrix::axis.break) if available
  .axis_break <- NULL
  if (exists("axis.break", mode = "function")) {
    .axis_break <- get("axis.break", mode = "function")
  } else if (requireNamespace("plotrix", quietly = TRUE)) {
    .axis_break <- plotrix::axis.break
  }
  
  # ---- input checks ----
  irep <- unique(as.integer(irep))
  irep <- irep[irep >= 1 & irep <= ncol(Zresidual)]
  if (!length(irep)) stop("irep does not select any valid column in `y`.")
  
  # ---- normality test vector (column-wise) ----
  test_vec <- NULL
  diagnosis.test <- if (missing(diagnosis.test) || is.null(diagnosis.test)) "SW" else diagnosis.test
  if (identical(diagnosis.test, "SW")) {
    test_vec <- .get_sw_vec(Zresidual)
  }
  
  for (col_id in irep) {
    graphics::par(mar = my.mar)
    
    z <- Zresidual[, col_id]
    
    # Track non-finite before replacement
    id_nonfinite <- which(!is.finite(z) & !is.na(z))
    id_nan <- which(is.nan(z))
    
    if (length(id_nonfinite) > 0L) {
      # Replace Inf/-Inf with large finite values near max finite magnitude
      finite_abs_max <- suppressWarnings(max(abs(z[is.finite(z)]), na.rm = TRUE))
      if (!is.finite(finite_abs_max)) finite_abs_max <- 0
      z[id_nonfinite] <- sign(z[id_nonfinite]) * (finite_abs_max + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id_nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
    
    # Outliers defined on ORIGINAL index (after non-finite replacement for plotting)
    id.outlier <- which(abs(z) > outlier.value & !is.na(z))
    
    # Q-Q data (order-statistic based)
    qq <- stats::qqnorm(z, plot.it = FALSE)
    x_values <- qq$x
    y_values <- qq$y
    
    # Map original indices -> QQ positions
    ord <- order(z, na.last = NA)          # indices of sorted z (drops NA)
    pos_out <- match(id.outlier, ord)      # positions in QQ vectors
    pos_out <- pos_out[!is.na(pos_out)]
    x_out <- if (length(pos_out)) x_values[pos_out] else numeric(0)
    y_out <- if (length(pos_out)) y_values[pos_out] else numeric(0)
    
    # Determine if we need compression for extremes
    ymax <- suppressWarnings(max(y_values, na.rm = TRUE))
    ymin <- suppressWarnings(min(y_values, na.rm = TRUE))
    need_clip <- isTRUE(clip.extreme) && (is.finite(ymax) && (ymax > clip.threshold) ||
                                            is.finite(ymin) && (ymin < -clip.threshold))
    
    # Build y for plotting (possibly clipped)
    y_plot <- y_values
    if (need_clip) {
      upper <- clip.threshold + 0.5
      lower <- -clip.threshold - 0.5
      y_plot[y_plot >  clip.threshold] <- upper
      y_plot[y_plot < -clip.threshold] <- lower
    }
    
    # --- Base plot (use plot(), not qqnorm()) so we control y clipping correctly ---
    default_plot_args <- list(x = x_values, y = y_plot,
                              main = main.title, xlab = xlab, ylab = ylab,
                              pch = 1)
    user_args <- list(...)
    # allow user to override pch/cex/etc
    plot_args <- modifyList(default_plot_args, user_args)
    do.call(graphics::plot, plot_args)
    
    # Reference lines
    stats::qqline(z, col = "black", lwd = 1.5)
    graphics::abline(a = 0, b = 1, col = "green")
    
    # Optional axis break marks / labels for clipped extremes
    if (need_clip) {
      # annotate actual extremes on y-axis
      usr <- graphics::par("usr")
      upper <- clip.threshold + 0.5
      lower <- -clip.threshold - 0.5
      
      # Add ticks at clipped positions with labels = true max/min
      if (is.finite(ymax) && ymax > clip.threshold) {
        graphics::axis(2, at = upper, labels = sprintf("%.1f", ymax), las = 1)
        graphics::axis(4, at = upper, labels = sprintf("%.1f", ymax), las = 1)
        if (!is.null(.axis_break)) {
          .axis_break(2, clip.threshold, style = "slash")
          .axis_break(4, clip.threshold, style = "slash")
        }
      }
      if (is.finite(ymin) && ymin < -clip.threshold) {
        graphics::axis(2, at = lower, labels = sprintf("%.1f", ymin), las = 1)
        graphics::axis(4, at = lower, labels = sprintf("%.1f", ymin), las = 1)
        if (!is.null(.axis_break)) {
          .axis_break(2, -clip.threshold, style = "slash")
          .axis_break(4, -clip.threshold, style = "slash")
        }
      }
    }
    
    # --- Legends (outside right) ---
    test.pv <- if (!is.null(test_vec) && length(test_vec) >= col_id) test_vec[col_id] else NA_real_
    
    default.legend1 <- list(
      legend = c("qqline", "45\u00B0 line"),
      col = c("black", "green"),
      lty = c(1, 1),
      lwd = 1.5,
      cex = 0.6,
      bty = "n",
      seg.len = 1.5,
      adj = 0
    )
    
    pv_str <- if (is.finite(test.pv)) sprintf("%.2f", test.pv) else "NA"
    default.legend2 <- list(
      legend = c(expression(bold("P-value:")),
                 paste0("Z-", diagnosis.test, " = ", pv_str)),
      cex = 0.6,
      bty = "n",
      adj = 0
    )
    
    legend.args1 <- modifyList(default.legend1, legend.settings)
    legend.args2 <- modifyList(default.legend2, legend.settings)
    
    .draw_legends_outside_right(legend.args1, legend.args2)
    
    # --- Outlier annotation (circles + labels) ---
    if (isTRUE(outlier.return) && length(id.outlier)) {
      # default outlier settings (need par("usr") after plot)
      default.outlier <- list(
        pos = 4,
        labels = id.outlier,
        cex = 0.8,
        col = "red",
        add = TRUE,
        inches = FALSE,
        circles = rep((graphics::par("usr")[2] - graphics::par("usr")[1]) * 0.027, length(id.outlier)),
        fg = "red"
      )
      outlier.args <- modifyList(default.outlier, outlier.set)
      
      text.args <- outlier.args[names(outlier.args) %in% names(formals(graphics::text.default))]
      symbols.args <- outlier.args[names(outlier.args) %in% names(formals(graphics::symbols))]
      
      # y positions must match plotted y (clipped if necessary)
      y_out_plot <- y_out
      if (need_clip && length(y_out_plot)) {
        upper <- clip.threshold + 0.5
        lower <- -clip.threshold - 0.5
        y_out_plot[y_out_plot >  clip.threshold] <- upper
        y_out_plot[y_out_plot < -clip.threshold] <- lower
      }
      
      # draw emphasis points
      graphics::points(x_out, y_out_plot, pch = 16)
      
      # circles + labels
      do.call(graphics::symbols, c(list(x = x_out, y = y_out_plot), symbols.args))
      do.call(graphics::text,    c(list(x = x_out, y = y_out_plot), text.args))
    }
    
    if (isTRUE(outlier.return)) {
      cat("Outlier Indices :", id.outlier, "\n")
      invisible(list(outliers = id.outlier))
    }
  }
}
