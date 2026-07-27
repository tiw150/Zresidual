#' Scatterplot diagnostics for Z-residuals
#'
#' Draws ggplot2 scatterplots of Z-residuals against observation index, a
#' linear predictor, a model covariate, or a user-supplied x-axis variable.
#' Multiple residual replications are displayed in facets.
#'
#' @param x A numeric vector or matrix of Z-residuals. Matrix columns are
#'   residual replications.
#' @param zcov Optional metadata returned by `Zcov()`.
#' @param irep Integer vector selecting residual columns.
#' @param ylab Label for the y-axis.
#' @param normality.test Any of `"SW"`, `"AOV"`, and `"BL"`.
#' @param k.test Number of groups used by AOV and Bartlett diagnostics.
#' @param x_axis_var `"index"`, `"lp"`, `"covariate"`, a covariate name,
#'   a length-n vector, or a function of `x` and `zcov`.
#' @param main.title Plot title.
#' @param outlier.return Whether to mark and label outliers.
#' @param outlier.value Absolute residual threshold for outliers.
#' @param category Optional grouping variable of length n.
#' @param outlier.set Named list controlling outlier colour, size, label size,
#'   and whether labels are shown. Supported names are `colour`, `size`,
#'   `label_colour`, `label_size`, and `label`.
#' @param xlab Optional x-axis label.
#' @param my.mar Retained for source compatibility; margins are controlled by
#'   the ggplot theme and this argument is ignored.
#' @param add_lowess Whether to add the original `stats::lowess()` smooth.
#' @param reference.lines Numeric horizontal reference lines. Their negatives
#'   are also drawn.
#' @param facet Whether multiple selected replications are shown in facets.
#'   When `FALSE`, replications are overlaid and distinguished by linetype only
#'   for the LOWESS curves; faceting is recommended.
#' @param x_scale Either `"identity"` or `"log10"`.
#' @param interactive If `TRUE`, convert the ggplot with `plotly::ggplotly()`.
#' @param point.args Named list of arguments passed to `geom_point()`.
#' @param smooth.args Named list of arguments passed to the LOWESS
#'   `geom_line()` layer.
#' @param theme A complete or partial ggplot2 theme.
#' @param legend.position Legend position passed to `ggplot2::theme()`. One of
#'   `"inside"`, `"bottom"`, `"top"`, `"left"`, `"right"`, or `"none"`.
#' @param legend.nrow Optional positive integer giving the number of rows in
#'   the combined colour-and-shape legend. Use `2L` with
#'   `legend.position = "bottom"` for a compact two-row legend.
#' @param ... Additional named arguments for the main `geom_point()` layer.
#'
#' @return A ggplot object, invisibly. With `interactive = TRUE`, returns a
#'   plotly htmlwidget. Diagnostics and outlier indices are stored in the
#'   `zresid_diagnostics` and `zresid_outliers` attributes of the ggplot.
#'
#' @method plot zresid
#' @export
plot.zresid <- function(x,
                        zcov = NULL,
                        irep = 1L,
                        ylab = "Z-Residual",
                        normality.test = c("SW", "AOV", "BL"),
                        k.test = 10L,
                        x_axis_var = "index",
                        main.title = NULL,
                        outlier.return = TRUE,
                        outlier.value = 3.5,
                        category = NULL,
                        outlier.set = list(),
                        xlab = NULL,
                        my.mar = c(5, 4, 4, 6) + 0.1,
                        add_lowess = FALSE,
                        reference.lines = c(1.96, 3),
                        facet = TRUE,
                        x_scale = c("identity", "log10"),
                        interactive = FALSE,
                        point.args = list(),
                        smooth.args = list(),
                        theme = ggplot2::theme_bw(),
                        legend.position = "inside",
                        legend.nrow = NULL,
                        ...) {
  x_scale <- match.arg(x_scale)
  
  if (!is.numeric(outlier.value) || length(outlier.value) != 1L ||
      is.na(outlier.value) || outlier.value <= 0) {
    stop("plot.zresid: outlier.value must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(k.test) || length(k.test) != 1L || is.na(k.test) ||
      k.test < 2 || k.test != as.integer(k.test)) {
    stop("plot.zresid: k.test must be one integer >= 2.", call. = FALSE)
  }
  if (!is.logical(interactive) || length(interactive) != 1L || is.na(interactive)) {
    stop("plot.zresid: interactive must be TRUE or FALSE.", call. = FALSE)
  }
  valid_legend_positions <- c(
    "inside",
    "bottom",
    "top",
    "left",
    "right",
    "none"
  )
  if (!is.character(legend.position) ||
      length(legend.position) != 1L ||
      is.na(legend.position) ||
      !(legend.position %in% valid_legend_positions)) {
    stop(
      "plot.zresid: legend.position must be one of ",
      paste(shQuote(valid_legend_positions), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.null(legend.nrow) &&
      (!is.numeric(legend.nrow) ||
       length(legend.nrow) != 1L ||
       is.na(legend.nrow) ||
       legend.nrow < 1 ||
       legend.nrow != as.integer(legend.nrow))) {
    stop(
      "plot.zresid: legend.nrow must be NULL or one positive integer.",
      call. = FALSE
    )
  }
  if (!is.null(legend.nrow)) {
    legend.nrow <- as.integer(legend.nrow)
  }
  
  prepared <- .zresid_prepare_plot(
    x = x,
    zcov = zcov,
    irep = irep,
    x_axis_var = x_axis_var,
    x_axis_expression = substitute(x_axis_var),
    category = category,
    category_label = paste(deparse(substitute(category)), collapse = ""),
    outlier_value = outlier.value,
    normality_test = normality.test,
    k_test = as.integer(k.test),
    add_lowess = add_lowess,
    x_scale = x_scale
  )
  
  residual_type <- attr(prepared$residuals, "type")
  if (is.null(main.title)) {
    main.title <- if (is.null(residual_type) || !length(residual_type)) {
      "Z-residual Scatterplot"
    } else {
      paste("Z-residual Scatterplot -", as.character(residual_type)[1L])
    }
  }
  
  dots <- list(...)
  old_names <- c(col = "colour", pch = "shape", cex = "size")
  for (old_name in names(old_names)) {
    if (!is.null(dots[[old_name]]) && is.null(dots[[old_names[[old_name]]]])) {
      dots[[old_names[[old_name]]]] <- dots[[old_name]]
      dots[[old_name]] <- NULL
    }
  }
  point_arguments <- utils::modifyList(list(size = 2, na.rm = TRUE), point.args)
  point_arguments <- utils::modifyList(point_arguments, dots)

  plot <- ggplot2::ggplot(
    prepared$data,
    ggplot2::aes(
      x = .data$x_value,
      y = .data$z_display,
      colour = .data$group,
      shape = .data$group,
      text = .data$tooltip
    )
  )
  plot <- plot + do.call(ggplot2::geom_point, point_arguments)
  
  plot <- plot +
    ggplot2::scale_colour_manual(
      values = prepared$groups$colours,
      name = prepared$groups$title,
      drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = prepared$groups$shapes,
      name = prepared$groups$title,
      drop = FALSE
    )
  
  if (!isTRUE(prepared$groups$show)) {
    plot <- plot + ggplot2::guides(colour = "none", shape = "none")
  } else if (!is.null(legend.nrow)) {
    plot <- plot + ggplot2::guides(
      colour = ggplot2::guide_legend(
        nrow = legend.nrow,
        byrow = TRUE
      ),
      shape = ggplot2::guide_legend(
        nrow = legend.nrow,
        byrow = TRUE
      )
    )
  }
  
  reference.lines <- unique(abs(reference.lines[is.finite(reference.lines)]))
  if (length(reference.lines)) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = c(reference.lines, -reference.lines),
      linetype = 3,
      colour = "grey60"
    )
  }
  
  lowess_data <- data.frame()
  if (isTRUE(add_lowess) &&
      (is.numeric(prepared$data$x_value) ||
       inherits(prepared$data$x_value, "Date") ||
       inherits(prepared$data$x_value, "POSIXt"))) {
    lowess_parts <- lapply(
      split(prepared$data, prepared$data$replication, drop = TRUE),
      function(current_data) {
        x_numeric <- as.numeric(current_data$x_value)
        y_numeric <- current_data$z_display
        keep <- is.finite(x_numeric) & is.finite(y_numeric)
        if (identical(x_scale, "log10")) {
          keep <- keep & x_numeric > 0
        }
        if (sum(keep) < 3L) return(NULL)

        smooth_x <- if (identical(x_scale, "log10")) {
          log10(x_numeric[keep])
        } else {
          x_numeric[keep]
        }
        lowess_result <- stats::lowess(
          x = smooth_x,
          y = y_numeric[keep]
        )
        plotted_x <- if (identical(x_scale, "log10")) {
          10^lowess_result$x
        } else {
          lowess_result$x
        }
        if (inherits(current_data$x_value, "Date")) {
          plotted_x <- as.Date(plotted_x, origin = "1970-01-01")
        }

        data.frame(
          replication = current_data$replication[1L],
          x_value = plotted_x,
          z_display = lowess_result$y
        )
      }
    )
    lowess_parts <- lowess_parts[lengths(lowess_parts) > 0L]
    if (length(lowess_parts)) lowess_data <- do.call(rbind, lowess_parts)
  }
  
  if (NROW(lowess_data)) {
    smooth_arguments <- utils::modifyList(
      list(
        data = lowess_data,
        mapping = ggplot2::aes(x = .data$x_value, y = .data$z_display, group = .data$replication),
        inherit.aes = FALSE,
        colour = "red",
        linewidth = 1.2,
        na.rm = TRUE
      ),
      smooth.args
    )
    plot <- plot + do.call(ggplot2::geom_line, smooth_arguments)
  }
  
  outlier_defaults <- list(
    colour = "darkolivegreen4",
    size = 4,
    label_colour = "darkolivegreen4",
    label_size = 3,
    label = TRUE
  )
  outlier_options <- utils::modifyList(outlier_defaults, outlier.set)
  
  if (isTRUE(outlier.return) && any(prepared$data$is_outlier, na.rm = TRUE)) {
    outlier_data <- prepared$data[prepared$data$is_outlier %in% TRUE, , drop = FALSE]
    plot <- plot + ggplot2::geom_point(
      data = outlier_data,
      mapping = ggplot2::aes(x = .data$x_value, y = .data$z_display),
      inherit.aes = FALSE,
      shape = 1,
      colour = outlier_options$colour,
      size = outlier_options$size,
      stroke = 1,
      na.rm = TRUE
    )
    if (isTRUE(outlier_options$label)) {
      plot <- plot + ggplot2::geom_text(
        data = outlier_data,
        mapping = ggplot2::aes(
          x = .data$x_value,
          y = .data$z_display,
          label = .data$observation
        ),
        inherit.aes = FALSE,
        colour = outlier_options$label_colour,
        size = outlier_options$label_size,
        fontface = "bold",
        nudge_y = 0.25,
        check_overlap = TRUE,
        na.rm = TRUE
      )
    }
  }
  
  infinite_data <- prepared$data[prepared$data$is_infinite %in% TRUE, , drop = FALSE]
  if (NROW(infinite_data)) {
    plot <- plot + ggplot2::geom_text(
      data = infinite_data,
      mapping = ggplot2::aes(
        x = .data$x_value,
        y = .data$z_display,
        label = .data$infinite_label
      ),
      inherit.aes = FALSE,
      colour = "red",
      nudge_x = 0.1,
      hjust = 0,
      na.rm = TRUE
    )
  }
  
  if (NROW(prepared$annotation)) {
    plot <- plot + ggplot2::geom_text(
      data = prepared$annotation,
      mapping = ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      hjust = -0.12,
      vjust = 1,
      size = 3.5,
      lineheight = 1.15
    )
  }
  
  if (isTRUE(facet) && length(unique(prepared$data$replication)) > 1L) {
    plot <- plot +
      ggplot2::facet_wrap(
        stats::as.formula("~ replication")
      )
  }
  
  if (identical(x_scale, "log10")) {
    if (!is.numeric(prepared$data$x_value)) {
      stop("plot.zresid: x_scale = 'log10' requires a numeric x-axis.", call. = FALSE)
    }
    if (any(prepared$data$x_value <= 0, na.rm = TRUE)) {
      stop("plot.zresid: x_scale = 'log10' requires positive x values.", call. = FALSE)
    }
    plot <- plot + ggplot2::scale_x_log10()
  }
  
  finite_y <- prepared$data$z_display[is.finite(prepared$data$z_display)]
  y_extent <- max(
    abs(finite_y),
    abs(reference.lines),
    na.rm = TRUE
  )
  if (!is.finite(y_extent) || y_extent <= 0) y_extent <- 1
  y_padding <- max(0.15, 0.04 * y_extent)

  plot <- plot +
    ggplot2::coord_cartesian(
      ylim = c(
        -y_extent - y_padding,
        y_extent + y_padding
      ),
      clip = "off"
    ) +
    ggplot2::labs(
      title = main.title,
      x = .zresid_resolve_label(
        xlab %zresid_or% prepared$x_info$label
      ),
      y = ylab,
      colour = prepared$groups$title,
      shape = prepared$groups$title
    ) +
    theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        hjust = 0.5
      ),
      plot.title.position = "panel",
      axis.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      
      legend.position = legend.position,
      legend.position.inside = c(1, 0.58),
      legend.justification = if (identical(legend.position, "inside")) {
        c(0, 0.5)
      } else {
        c(0.5, 0.5)
      },
      legend.background = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(
        t = 0,
        r = 0,
        b = 0,
        l = 4,
        unit = "pt"
      ),
      
      # The right margin contains the diagnostic annotation
      # and, by default, the inside category legend.
      plot.margin = ggplot2::margin(
        t = 8,
        r = 90,
        b = 8,
        l = 8,
        unit = "pt"
      )
    )

  attr(plot, "zresid_diagnostics") <- prepared$diagnostics
  attr(plot, "zresid_outliers") <- prepared$outliers
  attr(plot, "zresid_plot_data") <- prepared$data

  outlier_indices <- unique(unlist(prepared$outliers, use.names = FALSE))
  if (isTRUE(outlier.return) && length(outlier_indices)) {
    message("Outlier indices: ", paste(outlier_indices, collapse = ", "))
  }
  if (any(prepared$data$is_infinite)) {
    message("Infinite Z-residuals exist; finite display positions were used.")
  }
  if (any(prepared$data$is_nan)) {
    message("NaN Z-residuals exist and were omitted from plotted layers.")
  }

  if (isTRUE(interactive)) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      stop(
        "interactive = TRUE requires the optional package 'plotly'. Install it with install.packages('plotly').",
        call. = FALSE
      )
    }
    widget <- plotly::ggplotly(plot, tooltip = "text")
    print(widget)
    return(invisible(widget))
  }

  print(plot)
  invisible(plot)
}
