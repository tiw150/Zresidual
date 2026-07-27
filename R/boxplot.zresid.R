#' Boxplot diagnostics for Z-residuals
#'
#' Produces ggplot2 boxplots of Z-residuals grouped by binned x-axis values.
#' The original binning rules and SW, AOV, and Bartlett diagnostics are
#' retained.
#'
#' @param x A numeric vector or matrix of Z-residuals, with one column per
#'   residual replication.
#' @param zcov Optional metadata returned by \code{\link{Zcov}}.
#' @param info Legacy alias for \code{zcov}.
#' @param irep Integer vector selecting residual columns.
#' @param x_axis_var \code{"lp"}, \code{"covariate"}, a covariate name, a
#'   length-n vector, or a function returning such a vector.
#' @param num.bin Number of bins for a numeric x-axis variable.
#' @param normality.test Any of \code{"SW"}, \code{"AOV"}, and \code{"BL"}.
#' @param k.test Grouping parameter passed to the diagnostic functions.
#' @param main.title Main plot title.
#' @param outlier.return Whether to report and return observations satisfying
#'   the absolute Z-residual outlier rule.
#' @param outlier.value Absolute Z-residual outlier threshold.
#' @param interactive Convert each plot with \code{plotly::ggplotly()}.
#' @param theme A ggplot2 theme.
#' @param boxplot.args Named list passed to \code{ggplot2::geom_boxplot()}.
#' @param ... Additional named boxplot aesthetics. Common base names
#'   \code{col}, \code{lwd}, and \code{cex} are translated where possible.
#'
#' @return Invisibly returns a list containing \code{outliers} and
#'   \code{plots}. Multiple selected replications are printed sequentially.
#'
#' @method boxplot zresid
#' @export
boxplot.zresid <- function(x,
                           zcov = NULL,
                           info = NULL,
                           irep = 1,
                           x_axis_var = "lp",
                           num.bin = 10,
                           normality.test = c("SW", "AOV", "BL"),
                           k.test = 10,
                           main.title = ifelse(
                             is.null(attr(x, "type")),
                             "Z-residual Boxplot",
                             paste("Z-residual Boxplot -", attr(x, "type"))
                           ),
                           outlier.return = FALSE,
                           outlier.value = 3.5,
                           interactive = FALSE,
                           theme = ggplot2::theme_bw(),
                           boxplot.args = list(),
                           ...) {
  original_x <- x
  Zresidual <- if (is.null(dim(x))) matrix(x, ncol = 1L) else as.matrix(x)
  x_axis_expression <- substitute(x_axis_var)
  xlab_from_function <- NULL

  info0 <- if (!is.null(zcov)) zcov else info
  if (is.null(info0)) {
    info0 <- attr(original_x, "info")
    if (is.null(info0)) info0 <- attr(original_x, "zcov")
    if (is.null(info0)) info0 <- attr(original_x, "Zcov")
  }
  if (is.null(info0)) info0 <- list()

  .first_nonnull <- function(...) {
    values <- list(...)
    for (value in values) if (!is.null(value)) return(value)
    NULL
  }

  covariates <- .first_nonnull(
    info0$covariates,
    info0$covariate,
    attr(original_x, "covariates")
  )
  linear_predictor <- .first_nonnull(
    info0$linear_pred,
    info0$linear.pred,
    attr(original_x, "linear.pred"),
    attr(original_x, "linear_pred")
  )

  residual_type <- attr(original_x, "type")
  if (is.null(residual_type) && !is.null(info0$type)) {
    residual_type <- info0$type
  }
  if (is.null(residual_type) && !is.null(info0$y_type_kind)) {
    residual_type <- switch(
      as.character(info0$y_type_kind)[1L],
      censor = "survival",
      hurdle = "hurdle",
      trunc = "count",
      plain = NULL,
      NULL
    )
  }
  if (missing(main.title)) {
    main.title <- if (is.null(residual_type)) {
      "Z-residual Boxplot"
    } else {
      paste("Z-residual Boxplot -", residual_type)
    }
  }

  if (NROW(Zresidual) < 1L || NCOL(Zresidual) < 1L) {
    stop("boxplot.zresid: x must contain at least one observation.", call. = FALSE)
  }
  if (!is.numeric(num.bin) || length(num.bin) != 1L || is.na(num.bin) ||
      num.bin < 2 || num.bin != as.integer(num.bin)) {
    stop("boxplot.zresid: num.bin must be one integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(k.test) || length(k.test) != 1L || is.na(k.test) ||
      k.test < 2 || k.test != as.integer(k.test)) {
    stop("boxplot.zresid: k.test must be one integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(outlier.value) || length(outlier.value) != 1L ||
      is.na(outlier.value) || outlier.value <= 0) {
    stop("boxplot.zresid: outlier.value must be one positive number.", call. = FALSE)
  }

  irep <- unique(as.integer(irep))
  irep <- irep[!is.na(irep) & irep >= 1L & irep <= NCOL(Zresidual)]
  if (!length(irep)) stop("irep does not select a valid residual column.", call. = FALSE)

  .normalize_x_for_binning <- function(value) {
    if (is.null(value)) return(value)
    if (inherits(value, "POSIXt") || inherits(value, "Date")) return(as.numeric(value))
    if (is.logical(value)) return(as.integer(value))
    if (is.character(value)) {
      numeric_value <- suppressWarnings(as.numeric(value))
      if (all(!is.na(numeric_value))) return(numeric_value)
      return(factor(value))
    }
    if (is.factor(value)) return(droplevels(value))
    value
  }

  .calculate_bins <- function(fitted_value, number_bins) {
    value <- fitted_value
    if (is.factor(value)) return(droplevels(value))
    if (!is.numeric(value)) return(droplevels(as.factor(value)))

    bins <- droplevels(cut(value, number_bins))
    sparse_bins <- which(tapply(bins, bins, length) <= 2)
    has_effective_bins <- (nlevels(bins) - length(sparse_bins)) > 2

    if (!has_effective_bins) {
      if (all(value > 0, na.rm = TRUE)) {
        transformed_value <- log(value)
        message("Too few effective bins; fitted values converted to log for binning.")
      } else {
        transformed_value <- base::rank(
          value,
          na.last = "keep",
          ties.method = "average"
        )
        message(
          "Too few effective bins and non-positive values exist; using rank transform for binning."
        )
      }
      bins <- droplevels(cut(transformed_value, number_bins))
    }
    bins
  }

  X <- x_axis_var
  if (is.character(X) && length(X) > 1L) X <- X[1L]
  if (is.function(X)) {
    function_result <- X(original_x, info0)
    if (is.list(function_result) && !is.null(function_result$values)) {
      xlab_from_function <- function_result$label
      X <- function_result$values
    } else {
      X <- function_result
    }
  }

  choices <- c("lp", "covariate")
  n_observations <- NROW(Zresidual)
  is_user_vector <- length(X) == n_observations &&
    !(is.character(X) && length(X) == 1L && X %in% choices)
  is_keyword <- is.character(X) && length(X) == 1L && X %in% choices
  is_covariate_name <- is.character(X) && length(X) == 1L && !is_keyword

  if (!is_user_vector && !is_keyword && !is_covariate_name) {
    stop(
      "x_axis_var must be 'lp', 'covariate', a covariate name, a length-n vector, or a function.",
      call. = FALSE
    )
  }

  dots <- list(...)
  automatic_xlab <- function() {
    if (!is.null(dots$xlab)) return(dots$xlab)
    if (!is.null(xlab_from_function)) return(xlab_from_function)
    label <- paste(deparse(x_axis_expression), collapse = "")
    if (identical(label, "x_axis_var")) "X" else label
  }

  if (is_user_vector) {
    grouping_value <- .normalize_x_for_binning(X)
    x_label <- automatic_xlab()
    test_x <- X
  } else if (identical(X, "lp")) {
    if (is.null(linear_predictor)) {
      stop("linear-predictor metadata is missing; provide zcov or info.", call. = FALSE)
    }
    if (length(linear_predictor) != n_observations) {
      stop("linear-predictor length does not match the residual rows.", call. = FALSE)
    }
    grouping_value <- .normalize_x_for_binning(linear_predictor)
    x_label <- "Linear Predictor"
    test_x <- "lp"
  } else {
    if (is.null(covariates)) {
      stop("covariate metadata is missing; provide zcov or info.", call. = FALSE)
    }
    if (!is.data.frame(covariates)) covariates <- as.data.frame(covariates)
    if (NROW(covariates) != n_observations) {
      stop("covariate rows do not match the residual rows.", call. = FALSE)
    }
    covariate_names <- names(covariates)
    if (identical(X, "covariate")) {
      selected_covariate <- covariate_names[1L]
      cat(
        "To plot against other covariates, set x_axis_var to a covariate name. Available covariates:\n",
        paste(covariate_names, collapse = ", "), "\n"
      )
    } else {
      selected_covariate <- X
      if (!(selected_covariate %in% covariate_names)) {
        stop(
          "x_axis_var must be one of covariate names: ",
          paste(covariate_names, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
    grouping_value <- .normalize_x_for_binning(covariates[[selected_covariate]])
    x_label <- selected_covariate
    test_x <- selected_covariate
  }

  bins <- .calculate_bins(grouping_value, as.integer(num.bin))
  if (length(bins) != n_observations) {
    stop("resolved grouping variable length does not match residual rows.", call. = FALSE)
  }

  normality.test <- unique(toupper(normality.test))
  normality.test <- normality.test[normality.test %in% c("SW", "AOV", "BL")]

  boxplot_arguments <- utils::modifyList(
    list(
      width = 0.75,
      colour = "#1F2937",
      fill = "#DCEBFA",
      linewidth = 0.6,
      outlier.colour = "#D62728",
      outlier.shape = 1,
      outlier.size = 2,
      na.rm = TRUE
    ),
    boxplot.args
  )
  if (!is.null(dots$col)) boxplot_arguments$colour <- dots$col
  if (!is.null(dots$colour)) boxplot_arguments$colour <- dots$colour
  if (!is.null(dots$fill)) boxplot_arguments$fill <- dots$fill
  if (!is.null(dots$lwd)) boxplot_arguments$linewidth <- dots$lwd
  if (!is.null(dots$linewidth)) boxplot_arguments$linewidth <- dots$linewidth

  y_label <- if (!is.null(dots$ylab)) dots$ylab else "Z-Residual"
  plot_title <- if (!is.null(dots$main)) dots$main else main.title

  plots <- vector("list", length(irep))
  names(plots) <- paste0("Replicate ", irep)
  outliers <- stats::setNames(vector("list", length(irep)), names(plots))

  for (position in seq_along(irep)) {
    column_id <- irep[position]
    residual <- Zresidual[, column_id]
    original_residual <- residual

    id_nan <- which(is.nan(residual))
    id_infinite <- which(is.infinite(residual))
    id_outlier <- which(abs(residual) > outlier.value | is.infinite(residual))
    outliers[[position]] <- id_outlier

    if (length(id_infinite)) {
      finite_values <- abs(residual[is.finite(residual)])
      finite_max <- if (length(finite_values)) max(finite_values) else 0
      residual[id_infinite] <- sign(residual[id_infinite]) * (finite_max + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id_nan)) {
      message("NaNs exist! The model or the fitting process has a problem!")
    }

    plot_data <- data.frame(
      observation = seq_len(n_observations),
      bin = bins,
      residual = residual,
      residual_original = original_residual,
      is_threshold_outlier = seq_len(n_observations) %in% id_outlier
    )

    diagnostic_values <- numeric(0L)
    diagnostic_labels <- character(0L)
    test_function_names <- c(
      SW = "sw.test.zresid",
      AOV = "aov.test.zresid",
      BL = "bartlett.test.zresid"
    )
    for (test_name in normality.test) {
      function_name <- test_function_names[[test_name]]
      test_result <- tryCatch(
        {
          test_function <- get(function_name, mode = "function")
          if (identical(test_name, "SW")) {
            test_function(Zresidual)
          } else if (identical(test_name, "AOV")) {
            test_function(Zresidual, X = test_x, k.anova = k.test, zcov = info0)
          } else {
            test_function(Zresidual, X = test_x, k.bl = k.test, zcov = info0)
          }
        },
        error = function(error) rep(NA_real_, NCOL(Zresidual))
      )
      p_value <- suppressWarnings(as.numeric(test_result[column_id]))
      diagnostic_values[test_name] <- p_value
      diagnostic_labels <- c(
        diagnostic_labels,
        paste0(
          test_name,
          " - ",
          if (is.finite(p_value)) sprintf("%.2f", p_value) else "NA"
        )
      )
    }

    diagnostic_label <- paste(
      c("P-value(s):", diagnostic_labels),
      collapse = "\n"
    )

    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$bin, y = .data$residual)
    )
    plot <- plot + do.call(ggplot2::geom_boxplot, boxplot_arguments)
    plot <- plot +
      ggplot2::geom_hline(
        yintercept = 0,
        colour = "grey60",
        linetype = 3,
        linewidth = 0.5
      ) +
      ggplot2::geom_text(
        data = data.frame(label = diagnostic_label),
        mapping = ggplot2::aes(x = Inf, y = Inf, label = .data$label),
        inherit.aes = FALSE,
        hjust = -0.12,
        vjust = 1,
        size = 3.5,
        lineheight = 1.15
      )

    finite_values <- abs(residual[is.finite(residual)])
    finite_extent <- if (length(finite_values)) max(finite_values) else 0
    y_limit <- max(stats::qnorm(0.9999), finite_extent)
    y_limits <- c(-y_limit, y_limit + 1)
    if (!is.null(dots$ylim) && length(dots$ylim) == 2L) y_limits <- dots$ylim

    plot <- plot +
      ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
      ggplot2::labs(
        title = plot_title,
        x = x_label,
        y = y_label
      ) +
      theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.title.position = "panel",
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
        plot.margin = ggplot2::margin(8, 90, 8, 8)
      )

    attr(plot, "zresid_outliers") <- id_outlier
    attr(plot, "zresid_diagnostics") <- diagnostic_values
    attr(plot, "zresid_boxplot_data") <- plot_data
    plots[[position]] <- plot

    if (isTRUE(outlier.return)) {
      message("Outlier Indices : ", paste(id_outlier, collapse = ", "))
    }

    if (isTRUE(interactive)) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        stop("interactive = TRUE requires the optional package 'plotly'.", call. = FALSE)
      }
      widget <- plotly::ggplotly(plot)
      plots[[position]] <- widget
      print(widget)
    } else {
      print(plot)
    }
  }

  result <- list(outliers = outliers, plots = plots)
  return(invisible(result))
}
