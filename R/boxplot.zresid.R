#' Boxplot of Z-Residuals
#'
#' Produces a boxplot of Z-residuals grouped by binned fitted values or a selected
#' covariate. This diagnostic plot supports other count-data models (e.g. Bayesian
#' hurdle, zero-truncated) by visualizing residual distribution, detecting outliers,
#' and evaluating normality assumptions using Shapiro-Wilk, ANOVA, or Bartlett-type
#' tests for Z-residuals.
#'
#' @param x A matrix of Z-residuals with one column per MCMC iteration.
#'   Must contain attributes:
#'   \describe{
#'     \item{\code{"type"}}{Model type used to generate residuals (e.g., hurdle, truncated).}
#'     \item{\code{"fitted.value"}}{Vector of fitted values for the model.}
#'     \item{\code{"covariates"}}{Optional data frame of covariates.}
#'     \item{\code{"zero_id"}}{Indices of zero observations (if applicable).}
#'   }
#'
#' @param irep Integer vector indicating which columns of
#'   \code{Zresidual} to plot. Default is \code{1}.
#'
#' @param x_axis_var Character string specifying the x-axis variable:
#'   \itemize{
#'     \item \code{"fitted.value"} (default): Bin fitted values.
#'     \item \code{"covariate"}: Display a list of covariate names.
#'     \item A specific covariate name present in \code{attr(Zresidual, "covariates")}.
#'   }
#'
#' @param num.bin Integer. Number of bins for grouping fitted values or selected
#'   covariate. Defaults to \code{10}.
#'
#' @param normality.test Character vector specifying which normality tests to report:
#'   \itemize{
#'     \item \code{"SW"} - Shapiro-Wilk test for Z-residuals.
#'     \item \code{"AOV"} - ANOVA-based test for variance/mean structure.
#'     \item \code{"BL"} - Bartlett-type test for variance homogeneity.
#'   }
#'   Defaults to \code{c("SW","AOV","BL")}.
#'
#' @param k.test Integer. Number of groups to use for ANOVA/Bartlett-type tests.
#'   Default is \code{10}.
#'
#' @param main.title Character. Main title of the plot. Default includes the model
#'   type automatically.
#'
#' @param outlier.return Logical. If \code{TRUE}, returns the index of Z-residual
#'   values exceeding the threshold defined by \code{outlier.value}. Default is
#'   \code{FALSE}.
#'
#' @param outlier.value Numeric. Threshold for defining outliers based on absolute
#'   Z-residual magnitude. Default is \code{3.5}.
#'
#' @param ... Additional graphical parameters passed to \code{\link[graphics]{plot}}
#'   or \code{\link[graphics]{legend}}.
#'
#' @details
#' The function generates boxplots of Z-residuals across binned fitted values (or a
#' selected covariate), which helps detect lack of fit, heteroscedasticity,
#' and model misspecification.
#'
#' Infinite and non-finite residuals are automatically replaced with a maximal
#' finite value (with preserved sign), and a warning message is displayed.
#'
#' When \code{x_axis_var="covariate"}, users may supply any covariate name available in the
#' \code{"covariates"} attribute. If a covariate contains too few unique bins,
#' fitted values are transformed using \code{log()} to stabilize binning, with a
#' message provided.
#'
#' Normality diagnostics are displayed in the plot legend. Internally, the function
#' calls:
#' \code{sw.test.zresid()}, \code{aov.test.zresid()}, and
#' \code{bartlett.test.zresid()}.
#'
#' @return
#' If \code{outlier.return = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{outliers} - vector of indices where \code{|Zresidual| > outlier.value}.
#' }
#' Otherwise, the function returns \code{NULL} (invisible) and produces a diagnostic plot.
#'
#' @examples
#' \dontrun{
#' # Assuming 'zres' is a Z-residual matrix produced by a model-fitting function:
#' boxplot.zresid(zres)
#'
#' # Plot against a specific covariate
#' boxplot.zresid(zres, x_axis_var = "age")
#'
#' # Return outliers
#' box.out <- boxplot.zresid(zres, outlier.return = TRUE)
#' }
#'
#' @seealso
#' \code{\link[graphics]{plot}}, \code{\link[graphics]{boxplot}},
#' \code{\link[graphics]{legend}},
#' \code{sw.test.zresid}, \code{aov.test.zresid}, \code{bartlett.test.zresid}
#' @method boxplot zresid
#' @export
boxplot.zresid <- function(x, irep = 1,
                           x_axis_var = c("lp", "covariate"),
                           num.bin = 10,
                           normality.test = c("SW", "AOV", "BL"), k.test = 10,
                           main.title = paste("Z-residual Boxplot -", attr(x, "type")),
                           outlier.return = FALSE, outlier.value = 3.5,
                           ...) {
  Zresidual<-x
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

  calc.bin <- function(fitted.value, num.bin) {
    if (is.factor(fitted.value)) {
      bin <- fitted.value
    } else {
      bin <- droplevels(cut(fitted.value, num.bin))
      less2_factor <- which(tapply(bin, bin, length) <= 2)
      is.bins2 <- (nlevels(bin) - length(less2_factor)) > 2

      if (!is.bins2) {
        fitted.value <- log(fitted.value)
        message("Contrasts can be applied only to factors with 2 or more levels. Fitted values converted to log.")
        bin <- droplevels(cut(fitted.value, num.bin))
      }
    }
    bin
  }

  args <- list(...)

  if (missing(X)) X <- "lp"

  choices <- c("lp", "covariate")
  n_obs   <- NROW(Zresidual)

  is_uservec <- (length(X) == n_obs)                                   # user-supplied vector
  is_keyword <- is.character(X) && length(X) == 1 && X %in% choices     # "lp" or "covariate"
  is_covname <- is.character(X) && length(X) == 1 && !is_keyword        # e.g., "age"

  if (!is_uservec && !is_keyword && !is_covname) {
    stop("X must be: (1) a length-n vector, or (2) 'lp'/'covariate', or (3) a covariate name present in attr(Zresidual, 'covariates').")
  }

  if (is_uservec) {
    user_xv <- X
    if (is.logical(user_xv)) user_xv <- as.integer(user_xv)
    if (is.character(user_xv)) {
      xv_num <- suppressWarnings(as.numeric(user_xv))
      user_xv <- if (all(!is.na(xv_num))) xv_num else factor(user_xv)
    }
    if (is.factor(user_xv)) user_xv <- droplevels(user_xv)
  }

  for (j in irep) {

    type <- attr(Zresidual, "type")
    zero.id <- attr(Zresidual, "zero_id")
    id.nan <- which(is.nan(Zresidual[, j]))
    id.infinity <- which(is.infinite(Zresidual[, j]))
    id.outlier <- which(abs(Zresidual[, j]) > outlier.value | is.infinite(Zresidual[, j]))

    # Converting Inf/-Inf to maximum finite value
    if (length(id.infinity) > 0L) {
      value.notfinite <- as.character.na(Zresidual[, j][id.infinity])
      max.non.infinity <- max(abs(Zresidual[, j][-id.infinity]), na.rm = TRUE)
      Zresidual[, j][id.infinity] <- sign.na(Zresidual[, j][id.infinity]) * (max.non.infinity + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }

    if (length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")

    ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[, j]), na.rm = TRUE))

    test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
    current_test_pv <- NULL
    for (a in normality.test) {
      test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
      current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
    }
    test.legend <- modifyList(list(
      legend = c(expression(bold("P-value (s):")), current_test_pv),
      cex = 1, bty = "n", xpd = TRUE, adj = c(0, 0.5)
    ), args)
    test.legend <- test.legend[names(test.legend) %in% formalArgs(legend)]

    default.plot <- modifyList(list(ylab = "Z-Residual",
                                    ylim = c(-ylim0, ylim0 + 1),
                                    main = main.title), args)

    par(mar = c(5, 4, 4, 6) + 0.1)

    if (is_uservec) {
      xlab_user <- if (!is.null(args$xlab)) args$xlab else "X"
      do.call(plot, c(
        modifyList(list(
          calc.bin(user_xv, num.bin),
          Zresidual[, j],
          xlab = xlab_user
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
    # ------------------------------------------------------------------------

    if (X == "lp") {
      fitted.value <- attr(Zresidual, "linear.pred")
      do.call(plot, c(
        modifyList(list(
          calc.bin(fitted.value, num.bin),
          Zresidual[, j],
          xlab = "Fitted Value"
        ), default.plot)
      ))
      plot_limits <- par("usr")
      do.call(legend, c(list(
        x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01,
        y = plot_limits[4]
      ), test.legend))
    }

    if (X != "lp") {
      fitted.value <- attr(Zresidual, "covariates")
      if (X == "covariate") {
        i <- 1
        cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
            variable.names(fitted.value))
      } else if (X %in% variable.names(fitted.value)) {
        cov.name <- variable.names(fitted.value)
        i <- which(cov.name == X)
      } else {
        stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value), ". "))
      }

      do.call(plot, c(
        modifyList(list(
          calc.bin(fitted.value[, i], num.bin),
          Zresidual[, j],
          xlab = colnames(fitted.value)[i]
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




# boxplot.zresid <- function(Zresidual,irep=1,
#                            X = c("lp", "covariate"),
#                            num.bin = 10,
#                            normality.test = c("SW", "AOV", "BL"), k.test = 10,
#                            main.title=paste("Z-residual Boxplot -", attr(Zresidual, "type")),
#                            outlier.return = FALSE, outlier.value = 3.5,
#                            ...)
# {
#
#   sign.na <- function(x)
#   {
#     sign.x <- sign(x)
#     sign.x[is.infinite(x)] <- 1
#     sign.x
#   }
#
#   as.character.na <- function(x)
#   {
#     label.x <- as.character(x)
#     label.x[is.infinite(x)] <- "Inf"
#     label.x
#   }
#
#   calc.bin <- function(fitted.value, num.bin){
#     if (is.factor(fitted.value)) {
#       bin <- fitted.value
#     } else{
#       bin <- droplevels(cut(fitted.value, num.bin))
#       less2_factor<-which(tapply(bin,bin,length)<= 2)
#       is.bins2 <- (nlevels(bin) - length(less2_factor))>2
#
#       if(!is.bins2) {
#         fitted.value <- log(fitted.value)
#         message("Contrasts can be applied only to factors with 2 or more levels. Fitted values converted to log.")
#         bin <- droplevels(cut(fitted.value, num.bin))
#       }
#
#       bin
#     }
#   }
#
#   args <- list(...)
#
#   if (missing(X))
#     X = "lp"
#
#   for (j in irep) {
#
#
#     type <- attr(Zresidual, "type")
#     zero.id <- attr(Zresidual, "zero_id")
#     id.nan <- which(is.nan(Zresidual[,j]))
#     id.infinity <- which (is.infinite(Zresidual[,j]))
#     id.outlier <- which(abs(Zresidual[,j]) > outlier.value | is.infinite(Zresidual[,j]))
#
#     # Converting Inf/-Inf to maximum finite value
#     if (length(id.infinity) > 0L) {
#       value.notfinite <- as.character.na(Zresidual[,j][id.infinity])
#       max.non.infinity <- max(abs(Zresidual[,j][-id.infinity]), na.rm = T)
#       Zresidual[,j][id.infinity] <- sign.na(Zresidual[,j][id.infinity]) * (max.non.infinity + 0.1)
#       message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
#     }
#
#     if(length(id.nan) > 0L) message("NaNs exist! The model or the fitting process has a problem!")
#
#     ylim0 <- max(qnorm(c(0.9999)), max(abs(Zresidual[,j]), na.rm = T))
#
#     test.list <- c("SW" = "sw", "AOV" = "aov", "BL" = "bartlett")
#     current_test_pv <- NULL
#     for (a in normality.test) {
#       test <- get(paste0(test.list[a], ".test.zresid"))(Zresidual, X, k.test)
#       current_test_pv <- c(current_test_pv, paste(a, "-", sprintf("%3.2f", test[j])))
#     }
#
#     test.legend <- modifyList(list(legend = c(expression(bold("P-value (s):")), current_test_pv),
#                                    cex = 0.8, bty = "n", xpd = TRUE, adj = c(0, 0.5)), args)
#     test.legend <- test.legend[names(test.legend) %in% formalArgs(legend)]
#
#
#     default.plot <- modifyList(list(ylab = "Z-Residual", ylim = c(-ylim0, ylim0 + 1), main = main.title), args)
#
#     par(mar = c(5, 4, 4, 6) + 0.1)
#
#     if (X == "lp") {
#       fitted.value <- attr(Zresidual, "linear.pred")
#
#       do.call(plot, c(modifyList(list(calc.bin(fitted.value,num.bin), Zresidual[, j], xlab = "Fitted Value"), default.plot)))
#       plot_limits <- par("usr")
#       do.call(legend, c(list(x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01, y = plot_limits[4]), test.legend))
#     }
#
#     if (X != "lp") {
#       fitted.value <- attr(Zresidual, "covariates")
#       if(X == "covariate"){
#         i<-1
#         cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
#             variable.names(fitted.value))
#
#       } else if(X %in% variable.names(fitted.value)){
#         cov.name<-variable.names(fitted.value)
#         i<- which(cov.name==X)
#       } else{stop( paste0("X must be the one of covariate name: ", variable.names(fitted.value),". "))}
#
#       do.call(plot, c(modifyList(list(calc.bin(fitted.value[,i],num.bin), Zresidual[, j], xlab = colnames(fitted.value)[i]), default.plot)))
#       plot_limits <- par("usr")
#       do.call(legend, c(list(x = plot_limits[2] - (plot_limits[2] - plot_limits[1]) * 0.01, y = plot_limits[4]), test.legend))
#     }
#
#     if (outlier.return)
#     {
#       cat("Outlier Indices:", id.outlier, "\n")
#       invisible(list(outliers = id.outlier))
#     }
#
#   }
#   par(mar = c(5, 4, 4, 2) + 0.1)
# }
