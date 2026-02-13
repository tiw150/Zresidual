#' Boxplot of Z-Residuals
#'
#' Produces a boxplot of Z-residuals grouped by binned x-axis values (e.g., fitted
#' values / linear predictors, covariates, or any user-specified vector such as time).
#' This diagnostic plot supports count-data models (e.g., Bayesian hurdle, zero-truncated)
#' by visualizing residual distribution, detecting outliers, and reporting normality
#' diagnostics using Shapiro-Wilk, ANOVA, or Bartlett-type tests for Z-residuals.
#'
#' @usage
#' \method{boxplot}{zresid}(
#'   x,
#'   info = NULL,
#'   irep = 1,
#'   x_axis_var = "lp",
#'   num.bin = 10,
#'   normality.test = c("SW", "AOV", "BL"),
#'   k.test = 10,
#'   main.title = paste("Z-residual Boxplot -", attr(x, "type")),
#'   outlier.return = FALSE,
#'   outlier.value = 3.5,
#'   ...
#' )
#'
#' @param x
#' A numeric matrix of Z-residuals with one column per MCMC iteration (or replicate),
#' of class `"zresid"`. The method may use attributes such as `"type"`, `"zero_id"`,
#' `"covariates"`, and `"linear.pred"` if present.
#'
#' @param info
#' Optional metadata (typically output of `Zcov()` or an equivalent list). When provided,
#' this method can fill legacy attributes (e.g., `"covariates"`, `"linear.pred"`, `"zero_id"`,
#' `"type"`) to remain compatible with older plot workflows.
#'
#' @param irep
#' Integer vector indicating which columns of `x` to plot. Default is `1`.
#'
#' @param x_axis_var
#' Specifies x-axis values for binning. Supported inputs:
#' \itemize{
#'   \item A keyword: `"lp"` or `"covariate"`.
#'   \item A covariate name present in `attr(x, "covariates")`.
#'   \item A length-\eqn{n} vector (numeric / `Date` / `POSIXct` / factor / character), e.g., time.
#'   \item A function `function(z, info) ...` returning a length-\eqn{n} vector; optionally it can return
#'         `list(values = <vector>, label = <character>)` to set an automatic x-axis label.
#' }
#'
#' @param num.bin
#' Integer. Number of bins used for grouping numeric x-axis values. Defaults to `10`.
#' For factor-like x-axis values, binning is not applied (levels are used directly).
#'
#' @param normality.test
#' Character vector specifying which normality tests to report:
#' \itemize{
#'   \item `"SW"` - Shapiro-Wilk test for Z-residuals.
#'   \item `"AOV"` - ANOVA-based test for variance/mean structure.
#'   \item `"BL"` - Bartlett-type test for variance homogeneity.
#' }
#' Defaults to `c("SW","AOV","BL")`.
#'
#' @param k.test
#' Integer. Number of groups to use for ANOVA/Bartlett-type tests. Default is `10`.
#'
#' @param main.title
#' Character. Main title of the plot. Default includes model type automatically if available.
#'
#' @param outlier.return
#' Logical. If `TRUE`, returns the index of Z-residual values exceeding the threshold
#' defined by `outlier.value`. Default is `FALSE`.
#'
#' @param outlier.value
#' Numeric. Threshold for defining outliers based on absolute Z-residual magnitude.
#' Default is `3.5`.
#'
#' @param ...
#' Additional graphical parameters passed to `plot()` or `legend()`.
#'
#' @details
#' The function generates boxplots of Z-residuals across binned x-axis values, which helps
#' detect lack of fit, heteroscedasticity, and model misspecification.
#'
#' Infinite and non-finite residuals are automatically replaced with a maximal finite value
#' (with preserved sign), and a warning message is displayed.
#'
#' Normality diagnostics are displayed in the plot legend. Internally, the function calls:
#' `sw.test.zresid()`, `aov.test.zresid()`, and `bartlett.test.zresid()`.
#' If `x_axis_var` is non-numeric (e.g., `Date`, `POSIXct`, factor), some tests may return `NA`.
#'
#' @return
#' If `outlier.return = TRUE`, returns a list containing:
#' \itemize{
#'   \item `outliers` - vector of indices where `|Zresidual| > outlier.value`.
#' }
#' Otherwise, returns `NULL` (invisible) and produces a diagnostic plot.
#'
#' @examples
#' \dontrun{
#' boxplot.zresid(zres)                       # default: binned linear predictor
#' boxplot.zresid(zres, x_axis_var = "age")   # plot against a covariate
#' boxplot.zresid(zres, x_axis_var = timevec) # plot against time (length n)
#' box.out <- boxplot.zresid(zres, outlier.return = TRUE)
#' }
#'
#' @seealso
#' \code{\link[graphics]{plot}}, \code{\link[graphics]{boxplot}},
#' \code{\link[graphics]{legend}},
#' \code{sw.test.zresid}, \code{aov.test.zresid}, \code{bartlett.test.zresid}
#'
#' @method boxplot zresid
#' @export
boxplot.zresid <- function(x, info = NULL, irep = 1,
                           x_axis_var = "lp",
                           num.bin = 10,
                           normality.test = c("SW", "AOV", "BL"), k.test = 10,
                           main.title = ifelse(is.null(attr(x, "type")),
                                               "Z-residual Boxplot",
                                               paste("Z-residual Boxplot -", attr(x, "type"))),
                           outlier.return = FALSE, outlier.value = 3.5,
                           ...) {
  
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
  }
  
  # If caller didn't set main.title, recompute after we possibly filled attr(type)
  if (missing(main.title)) {
    main.title <- ifelse(
      is.null(attr(Zresidual, "type")),
      "Z-residual Boxplot",
      paste("Z-residual Boxplot -", attr(Zresidual, "type"))
    )
  }
  
  args <- list(...)
  
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
  
  # Normalize various x types for binning:
  # - factor/character/logical: factor
  # - Date/POSIXct: numeric internally for cut(), but bins are factors anyway
  .normalize_x_for_binning <- function(v) {
    if (is.null(v)) return(v)
    if (inherits(v, "POSIXt") || inherits(v, "Date")) return(as.numeric(v))
    if (is.logical(v)) return(as.integer(v))
    if (is.character(v)) {
      vv <- suppressWarnings(as.numeric(v))
      return(if (all(!is.na(vv))) vv else factor(v))
    }
    if (is.factor(v)) return(droplevels(v))
    v
  }
  
  calc.bin <- function(fitted.value, num.bin) {
    fv <- fitted.value
    
    if (is.factor(fv)) return(droplevels(fv))
    
    # if it's non-numeric (e.g., character), treat as factor bins
    if (!is.numeric(fv)) return(droplevels(as.factor(fv)))
    
    # numeric binning
    bin <- droplevels(cut(fv, num.bin))
    less2_factor <- which(tapply(bin, bin, length) <= 2)
    is.bins2 <- (nlevels(bin) - length(less2_factor)) > 2
    
    if (!is.bins2) {
      # fallback: log if possible, else rank transform
      if (all(fv > 0, na.rm = TRUE)) {
        fv2 <- log(fv)
        message("Too few effective bins; fitted values converted to log for binning.")
      } else {
        fv2 <- stats::rank(fv, na.last = "keep", ties.method = "average")
        message("Too few effective bins and non-positive values exist; using rank transform for binning.")
      }
      bin <- droplevels(cut(fv2, num.bin))
    }
    
    bin
  }
  
  # ---- x_axis_var: allow keyword / covariate name / length-n vector / function ----
  X <- x_axis_var
  if (is.character(X) && length(X) > 1L) X <- X[1L]
  
  if (is.function(X)) {
    xr <- X(Zresidual, info0)
    if (is.list(xr) && !is.null(xr$values)) {
      xlab_from_xaxis <- xr$label
      X <- xr$values
    } else {
      X <- xr
    }
  }
  
  choices <- c("lp", "covariate")
  n_obs   <- NROW(Zresidual)
  
  is_uservec <- (length(X) == n_obs)
  is_keyword <- is.character(X) && length(X) == 1 && X %in% choices
  is_covname <- is.character(X) && length(X) == 1 && !is_keyword
  
  if (!is_uservec && !is_keyword && !is_covname) {
    stop("x_axis_var must be one of: 'lp','covariate', a covariate name in attr(z,'covariates'), ",
         "a length-n vector (e.g., time), or a function(z, info) returning a length-n vector.")
  }
  
  # best-effort xlab for user vector mode
  xlab_auto_uservec <- function() {
    if (!is.null(args$xlab)) return(args$xlab)
    if (!is.null(xlab_from_xaxis)) return(xlab_from_xaxis)
    lab <- paste(deparse(x_axis_expr), collapse = "")
    if (identical(lab, "x_axis_var")) "X" else lab
  }
  
  for (j in irep) {
    
    id.nan <- which(is.nan(Zresidual[, j]))
    id.infinity <- which(is.infinite(Zresidual[, j]))
    id.outlier <- which(abs(Zresidual[, j]) > outlier.value | is.infinite(Zresidual[, j]))
    
    if (length(id.infinity) > 0L) {
      value.notfinite <- as.character.na(Zresidual[, j][id.infinity])
      max.non.infinity <- max(abs(Zresidual[, j][-id.infinity]), na.rm = TRUE)
      Zresidual[, j][id.infinity] <- sign.na(Zresidual[, j][id.infinity]) * (max.non.infinity + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
    
    ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[, j]), na.rm = TRUE))
    
    # normality tests (robust to non-numeric x; failures become NA)
    test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
    current_test_pv <- NULL
    for (a in normality.test) {
      testfun <- get(paste0(test.list[a], ".test.zresid"))
      tt <- try(testfun(Zresidual, X, k.test), silent = TRUE)
      pv_str <- if (inherits(tt, "try-error")) "NA" else sprintf("%3.2f", tt[j])
      current_test_pv <- c(current_test_pv, paste(a, "-", pv_str))
    }
    
    test.legend <- modifyList(list(
      legend = c(expression(bold("P-value (s):")), current_test_pv),
      cex = 1, bty = "n", xpd = TRUE, adj = c(0, 0.5)
    ), args)
    test.legend <- test.legend[names(test.legend) %in% formalArgs(legend)]
    
    default.plot <- modifyList(list(
      ylab = "Z-Residual",
      ylim = c(-ylim0, ylim0 + 1),
      main = main.title
    ), args)
    
    par(mar = c(5, 4, 4, 6) + 0.1)
    
    # --------------------------
    # USER VECTOR MODE
    # --------------------------
    if (is_uservec && !(is.character(X) && length(X) == 1 && X %in% choices)) {
      user_xv <- .normalize_x_for_binning(X)
      bin <- calc.bin(user_xv, num.bin)
      
      do.call(plot, c(
        modifyList(list(
          bin,
          Zresidual[, j],
          xlab = xlab_auto_uservec()
        ), default.plot)
      ))
      
      plot_limits <- par("usr")
      do.call(legend, c(list(
        x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01,
        y = plot_limits[4]
      ), test.legend))
      
      if (outlier.return) {
        cat("Outlier Indices:", id.outlier, "\n")
        invisible(list(outliers = id.outlier))
      }
      next
    }
    
    # --------------------------
    # KEYWORD / COVARIATE MODE
    # --------------------------
    if (identical(X, "lp")) {
      fitted.value <- attr(Zresidual, "linear.pred")
      if (is.null(fitted.value)) stop("attr(x,'linear.pred') is missing; provide `info=` or attach the attribute.")
      fv <- .normalize_x_for_binning(fitted.value)
      
      do.call(plot, c(
        modifyList(list(
          calc.bin(fv, num.bin),
          Zresidual[, j],
          xlab = "Linear Predictor"
        ), default.plot)
      ))
      
      plot_limits <- par("usr")
      do.call(legend, c(list(
        x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01,
        y = plot_limits[4]
      ), test.legend))
      
    } else {
      covs <- attr(Zresidual, "covariates")
      if (is.null(covs)) stop("attr(x,'covariates') is missing; provide `info=` or attach the attribute.")
      
      if (identical(X, "covariate")) {
        i <- 1
        cat("To plot against other covariates, set x_axis_var to a covariate name. Available covariates:\n",
            paste(variable.names(covs), collapse = ", "), "\n")
      } else if (X %in% variable.names(covs)) {
        cov.name <- variable.names(covs)
        i <- which(cov.name == X)
      } else {
        stop(paste0("x_axis_var must be one of covariate names: ", paste(variable.names(covs), collapse = ", "), "."))
      }
      
      xv <- .normalize_x_for_binning(covs[, i])
      do.call(plot, c(
        modifyList(list(
          calc.bin(xv, num.bin),
          Zresidual[, j],
          xlab = colnames(covs)[i]
        ), default.plot)
      ))
      
      plot_limits <- par("usr")
      do.call(legend, c(list(
        x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01,
        y = plot_limits[4]
      ), test.legend))
    }
    
    if (outlier.return) {
      cat("Outlier Indices:", id.outlier, "\n")
      invisible(list(outliers = id.outlier))
    }
  }
  
  par(mar = c(5, 4, 4, 2) + 0.1)
}
