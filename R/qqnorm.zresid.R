#' Normal Q-Q plot for Z-residuals
#'
#' Produces a ggplot2 normal Q-Q plot for one or more columns of a
#' \code{"zresid"} object. The Shapiro-Wilk diagnostic, outlier annotation,
#' extreme-value compression, reference lines, and axis-break marks retain the
#' behavior of the original plotting method.
#'
#' @param y A numeric matrix of Z-residuals, typically returned by
#'   \code{\link{Zresidual}}, with one column per residual replicate.
#' @param zcov Optional metadata, typically returned by \code{\link{Zcov}}.
#' @param info Legacy alias for \code{zcov}.
#' @param irep Integer vector specifying which column(s) of \code{y} to plot.
#' @param diagnosis.test Character string specifying the normality diagnostic.
#'   Currently \code{"SW"} is supported.
#' @param main.title Main title of the plot.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param outlier.return Logical; mark and return observations satisfying the
#'   outlier rule when \code{TRUE}.
#' @param outlier.value Numeric absolute-residual outlier threshold.
#' @param outlier.set Named list controlling outlier annotation. Supported
#'   ggplot settings are \code{colour}, \code{size}, \code{label_colour},
#'   \code{label_size}, and \code{label}.
#' @param outlier.label.xpad Numeric padding added to the right x-axis limit.
#' @param my.mar Retained for source compatibility; ggplot themes control plot
#'   margins.
#' @param legend.settings Named list modifying the line legend. Supported names
#'   include \code{position}, \code{text_size}, and \code{title_size}.
#' @param clip.extreme Logical; visually compress extreme residuals.
#' @param clip.threshold Numeric clipping threshold.
#' @param interactive Logical; convert each ggplot using
#'   \code{plotly::ggplotly()} when \code{TRUE}.
#' @param theme A ggplot2 theme.
#' @param ... Additional named arguments passed to the QQ-point
#'   \code{geom_point()} layer.
#'
#' @return Invisibly returns a list containing \code{outliers} and
#'   \code{plots}. A single selected replication produces one ggplot in
#'   \code{plots}; multiple selected replications are printed sequentially,
#'   matching the original method.
#'
#' @method qqnorm zresid
#' @export
qqnorm.zresid <- function(y,
                          zcov = NULL,
                          info = NULL,
                          irep = 1,
                          diagnosis.test = "SW",
                          main.title = ifelse(
                            is.null(attr(y, "type")),
                            "Normal Q-Q Plot",
                            paste("Normal Q-Q Plot -", attr(y, "type"))
                          ),
                          xlab = "Theoretical Quantiles",
                          ylab = "Sample Quantiles",
                          outlier.return = TRUE,
                          outlier.value = 3.5,
                          outlier.set = list(),
                          outlier.label.xpad = 0.5,
                          my.mar = c(5, 4, 4, 6) + 0.1,
                          legend.settings = list(),
                          clip.extreme = TRUE,
                          clip.threshold = 6,
                          interactive = FALSE,
                          theme = ggplot2::theme_bw(),
                          ...) {
  Zresidual <- if (is.null(dim(y))) matrix(y, ncol = 1L) else as.matrix(y)

  info0 <- if (!is.null(zcov)) zcov else info
  if (is.null(info0)) {
    info0 <- attr(y, "info")
    if (is.null(info0)) info0 <- attr(y, "zcov")
    if (is.null(info0)) info0 <- attr(y, "Zcov")
  }

  .map_type_from_info <- function(info_object) {
    if (is.null(info_object)) return(NULL)
    if (!is.null(info_object$type)) return(info_object$type)
    switch(
      as.character(info_object$y_type_kind)[1L],
      censor = "survival",
      hurdle = "hurdle",
      trunc = "count",
      plain = NULL,
      NULL
    )
  }

  residual_type <- attr(y, "type")
  if (is.null(residual_type)) residual_type <- .map_type_from_info(info0)
  if (missing(main.title)) {
    main.title <- if (is.null(residual_type)) {
      "Normal Q-Q Plot"
    } else {
      paste("Normal Q-Q Plot -", residual_type)
    }
  }

  .safe_sw_pvalue <- function(z_column) {
    finite_z <- z_column[is.finite(z_column)]
    if (length(finite_z) < 3L) return(NA_real_)
    result <- try(stats::shapiro.test(finite_z)$p.value, silent = TRUE)
    if (inherits(result, "try-error")) NA_real_ else result
  }

  .get_sw_vector <- function(z_matrix) {
    if (exists("sw.test.zresid", mode = "function")) {
      result <- try(sw.test.zresid(z_matrix), silent = TRUE)
      if (!inherits(result, "try-error")) return(result)
    }
    apply(z_matrix, 2L, .safe_sw_pvalue)
  }

  outlier.label.xpad <- suppressWarnings(as.numeric(outlier.label.xpad)[1L])
  if (!is.finite(outlier.label.xpad) || outlier.label.xpad < 0) {
    outlier.label.xpad <- 0.5
  }
  if (!is.numeric(outlier.value) || length(outlier.value) != 1L ||
      is.na(outlier.value) || outlier.value <= 0) {
    stop("qqnorm.zresid: outlier.value must be one positive number.", call. = FALSE)
  }
  if (!is.numeric(clip.threshold) || length(clip.threshold) != 1L ||
      is.na(clip.threshold) || clip.threshold <= 0) {
    stop("qqnorm.zresid: clip.threshold must be one positive number.", call. = FALSE)
  }

  irep <- unique(as.integer(irep))
  irep <- irep[!is.na(irep) & irep >= 1L & irep <= NCOL(Zresidual)]
  if (!length(irep)) stop("irep does not select any valid column in y.", call. = FALSE)

  diagnosis.test <- if (missing(diagnosis.test) || is.null(diagnosis.test)) {
    "SW"
  } else {
    diagnosis.test
  }
  test_vector <- if (identical(diagnosis.test, "SW")) {
    .get_sw_vector(Zresidual)
  } else {
    NULL
  }

  dots <- list(...)
  old_names <- c(col = "colour", pch = "shape", cex = "size")
  for (old_name in names(old_names)) {
    if (!is.null(dots[[old_name]]) && is.null(dots[[old_names[[old_name]]]])) {
      dots[[old_names[[old_name]]]] <- dots[[old_name]]
      dots[[old_name]] <- NULL
    }
  }
  point_arguments <- utils::modifyList(
    list(shape = 1, size = 2, colour = "black", na.rm = TRUE),
    dots
  )

  outlier_defaults <- list(
    colour = "red",
    size = 4,
    label_colour = "red",
    label_size = 3,
    label = TRUE
  )
  outlier_options <- utils::modifyList(outlier_defaults, outlier.set)

  plots <- vector("list", length(irep))
  names(plots) <- paste0("Replicate ", irep)
  outliers <- stats::setNames(vector("list", length(irep)), names(plots))

  for (position in seq_along(irep)) {
    column_id <- irep[position]
    z <- Zresidual[, column_id]

    id_nonfinite <- which(!is.finite(z) & !is.na(z))
    id_nan <- which(is.nan(z))
    if (length(id_nonfinite)) {
      finite_max <- suppressWarnings(max(abs(z[is.finite(z)]), na.rm = TRUE))
      if (!is.finite(finite_max)) finite_max <- 0
      z[id_nonfinite] <- sign(z[id_nonfinite]) * (finite_max + 0.1)
      message("Non-finite Zresiduals exist! The model or the fitting process has a problem!")
    }
    if (length(id_nan)) {
      message("NaNs exist! The model or the fitting process has a problem!")
    }

    id_outlier <- which(abs(z) > outlier.value & !is.na(z))
    outliers[[position]] <- id_outlier

    qq_values <- stats::qqnorm(z, plot.it = FALSE)
    x_values <- qq_values$x
    y_values <- qq_values$y

    # Preserve the original observation-to-QQ-position mapping.
    kept_indices <- which(!is.na(z))
    outlier_positions <- match(id_outlier, kept_indices)
    outlier_positions <- outlier_positions[!is.na(outlier_positions)]

    ymax <- suppressWarnings(max(y_values, na.rm = TRUE))
    ymin <- suppressWarnings(min(y_values, na.rm = TRUE))
    need_clip <- isTRUE(clip.extreme) &&
      ((is.finite(ymax) && ymax > clip.threshold) ||
       (is.finite(ymin) && ymin < -clip.threshold))

    display_threshold <- clip.threshold

    upper <- display_threshold + 0.5
    lower <- -display_threshold - 0.5
    y_plot <- y_values
    if (need_clip) {
      y_plot[y_plot > clip.threshold] <- upper
      y_plot[y_plot < -clip.threshold] <- lower
    }

    y_limits <- range(y_plot, finite = TRUE)
    y_padding <- max(0.15, 0.04 * diff(y_limits))
    y_limits <- y_limits + c(-y_padding, y_padding)

    plot_data <- data.frame(
      theoretical = x_values,
      sample_original = y_values,
      sample_display = y_plot,
      qq_position = seq_along(x_values),
      tooltip = paste0(
        "Replicate: ", column_id,
        "<br>Theoretical quantile: ", signif(x_values, 5),
        "<br>Sample quantile: ", signif(y_values, 5)
      )
    )

    x_limits <- range(x_values, finite = TRUE)
    if (!all(is.finite(x_limits))) x_limits <- c(-3, 3)
    x_limits[2L] <- x_limits[2L] + outlier.label.xpad

    quartiles_y <- stats::quantile(z, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    quartiles_x <- stats::qnorm(c(0.25, 0.75))
    qq_slope <- diff(quartiles_y) / diff(quartiles_x)
    qq_intercept <- quartiles_y[1L] - qq_slope * quartiles_x[1L]

    plot <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        x = .data$theoretical,
        y = .data$sample_display,
        text = .data$tooltip
      )
    )
    plot <- plot + do.call(ggplot2::geom_point, point_arguments)

    .clip_reference_line <- function(
    intercept,
    slope,
    x_limits,
    y_limits
    ) {
      candidate_x <- c(
        x_limits,
        (y_limits - intercept) / slope
      )
      
      candidate_y <- intercept + slope * candidate_x

      keep <- candidate_x >= x_limits[1L] &
        candidate_x <= x_limits[2L] &
        candidate_y >= y_limits[1L] &
        candidate_y <= y_limits[2L] &
        is.finite(candidate_x) &
        is.finite(candidate_y)

      candidate_x <- candidate_x[keep]
      candidate_y <- candidate_y[keep]

      order_x <- order(candidate_x)

      data.frame(
        x = candidate_x[order_x][1L],
        xend = candidate_x[order_x][length(order_x)],
        y = candidate_y[order_x][1L],
        yend = candidate_y[order_x][length(order_x)]
      )
    }
    
    qqline_segment <- .clip_reference_line(
      intercept = qq_intercept,
      slope = qq_slope,
      x_limits = x_limits,
      y_limits = y_limits
    )
    
    degree45_segment <- .clip_reference_line(
      intercept = 0,
      slope = 1,
      x_limits = x_limits,
      y_limits = y_limits
    )
    
    reference_data <- rbind(
      transform(
        qqline_segment,
        reference = "qqline"
      ),
      transform(
        degree45_segment,
        reference = "45\u00B0 line"
      )
    )
    
    reference_data$reference <- factor(
      reference_data$reference,
      levels = c("qqline", "45\u00B0 line")
    )
    
    plot <- plot +
      ggplot2::geom_segment(
        data = reference_data,
        mapping = ggplot2::aes(
          x = .data$x,
          xend = .data$xend,
          y = .data$y,
          yend = .data$yend,
          colour = .data$reference,
          linetype = .data$reference
        ),
        inherit.aes = FALSE,
        linewidth = 0.9,
        show.legend = TRUE
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          "qqline" = "#E41A1C",
          "45\u00B0 line" = "#0066FF"
        ),
        breaks = c(
          "qqline",
          "45\u00B0 line"
        ),
        name = NULL
      ) +
      ggplot2::scale_linetype_manual(
        values = c(
          "qqline" = "solid",
          "45\u00B0 line" = "dotdash"
        ),
        breaks = c(
          "qqline",
          "45\u00B0 line"
        ),
        name = NULL
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(order = 1),
        linetype = ggplot2::guide_legend(order = 1)
      )

    if (isTRUE(outlier.return) && length(outlier_positions)) {
      outlier_data <- plot_data[outlier_positions, , drop = FALSE]
      plot <- plot +
        ggplot2::geom_point(
          data = outlier_data,
          mapping = ggplot2::aes(x = .data$theoretical, y = .data$sample_display),
          inherit.aes = FALSE,
          shape = 16,
          colour = "black",
          size = 2,
          na.rm = TRUE
        ) +
        ggplot2::geom_point(
          data = outlier_data,
          mapping = ggplot2::aes(x = .data$theoretical, y = .data$sample_display),
          inherit.aes = FALSE,
          shape = 1,
          colour = outlier_options$colour,
          size = outlier_options$size,
          stroke = 1,
          na.rm = TRUE
        )
      if (isTRUE(outlier_options$label)) {
        outlier_data$observation <- id_outlier[seq_len(NROW(outlier_data))]
        plot <- plot + ggplot2::geom_text(
          data = outlier_data,
          mapping = ggplot2::aes(
            x = .data$theoretical,
            y = .data$sample_display,
            label = .data$observation
          ),
          inherit.aes = FALSE,
          colour = outlier_options$label_colour,
          size = outlier_options$label_size,
          hjust = 0,
          vjust = 0,
          nudge_x = 0.08,
          nudge_y = 0.16,
          check_overlap = TRUE,
          na.rm = TRUE
        )
      }
    }

    p_value <- if (!is.null(test_vector) && length(test_vector) >= column_id) {
      suppressWarnings(as.numeric(test_vector[column_id]))
    } else {
      NA_real_
    }
    p_value_string <- if (is.finite(p_value)) sprintf("%.2f", p_value) else "NA"
    test_label <- paste0("P-value:\nZ-", diagnosis.test, " = ", p_value_string)
    test_data <- data.frame(label = test_label)
    plot <- plot + ggplot2::geom_text(
      data = test_data,
      mapping = ggplot2::aes(x = Inf, y = Inf, label = .data$label),
      inherit.aes = FALSE,
      hjust = -0.12,
      vjust = 1,
      size = legend.settings$text_size %zresid_or% 3.5,
      lineheight = 1.15
    )

    y_breaks <- ggplot2::waiver()
    y_labels <- ggplot2::waiver()
    if (need_clip) {
      base_breaks <- pretty(range(y_plot, finite = TRUE))
      base_breaks <- base_breaks[base_breaks >= min(y_plot, na.rm = TRUE) &
                                   base_breaks <= max(y_plot, na.rm = TRUE)]
      y_breaks <- base_breaks
      y_labels <- format(base_breaks, trim = TRUE)

      if (is.finite(ymax) && ymax > clip.threshold) {
        y_breaks <- c(y_breaks, upper)
        y_labels <- c(y_labels, sprintf("%.1f", ymax))
      }
      if (is.finite(ymin) && ymin < -clip.threshold) {
        y_breaks <- c(y_breaks, lower)
        y_labels <- c(y_labels, sprintf("%.1f", ymin))
      }
      order_breaks <- order(y_breaks)
      y_breaks <- y_breaks[order_breaks]
      y_labels <- y_labels[order_breaks]

      x_range <- diff(x_limits)
      slash_dx <- 0.018 * x_range
      slash_dy <- 0.12
      break_rows <- list()
      if (is.finite(ymax) && ymax > clip.threshold) {
        break_rows[[length(break_rows) + 1L]] <- data.frame(
          x = c(x_limits[1L] - slash_dx, x_limits[2L] - slash_dx),
          xend = c(x_limits[1L] + slash_dx, x_limits[2L] + slash_dx),
          y = display_threshold - slash_dy,
          yend = display_threshold + slash_dy
        )
      }
      if (is.finite(ymin) && ymin < -clip.threshold) {
        break_rows[[length(break_rows) + 1L]] <- data.frame(
          x = c(x_limits[1L] - slash_dx, x_limits[2L] - slash_dx),
          xend = c(x_limits[1L] + slash_dx, x_limits[2L] + slash_dx),
          y = -display_threshold - slash_dy,
          yend = -display_threshold + slash_dy
        )
      }

      break_positions <- numeric(0L)
      if (is.finite(ymax) && ymax > clip.threshold) {
        break_positions <- c(break_positions, display_threshold)
      }
      if (is.finite(ymin) && ymin < -clip.threshold) {
        break_positions <- c(break_positions, -display_threshold)
      }
      break_positions <- sort(unique(break_positions))
      border_gap <- slash_dy * 1.35

      vertical_start <- c(
        y_limits[1L],
        break_positions + border_gap
      )
      vertical_end <- c(
        break_positions - border_gap,
        y_limits[2L]
      )
      keep_border <- vertical_end > vertical_start
      vertical_start <- vertical_start[keep_border]
      vertical_end <- vertical_end[keep_border]

      vertical_border <- data.frame(
        x = rep(x_limits, each = length(vertical_start)),
        xend = rep(x_limits, each = length(vertical_start)),
        y = rep(vertical_start, times = 2L),
        yend = rep(vertical_end, times = 2L)
      )
      horizontal_border <- data.frame(
        x = x_limits[1L],
        xend = x_limits[2L],
        y = y_limits,
        yend = y_limits
      )

      plot <- plot +
        ggplot2::geom_segment(
          data = vertical_border,
          mapping = ggplot2::aes(
            x = .data$x,
            xend = .data$xend,
            y = .data$y,
            yend = .data$yend
          ),
          inherit.aes = FALSE,
          linewidth = 0.7,
          colour = "black"
        ) +
        ggplot2::geom_segment(
          data = horizontal_border,
          mapping = ggplot2::aes(
            x = .data$x,
            xend = .data$xend,
            y = .data$y,
            yend = .data$yend
          ),
          inherit.aes = FALSE,
          linewidth = 0.7,
          colour = "black"
        )

      if (length(break_rows)) {
        break_data <- do.call(rbind, break_rows)
        plot <- plot + ggplot2::geom_segment(
          data = break_data,
          mapping = ggplot2::aes(
            x = .data$x,
            xend = .data$xend,
            y = .data$y,
            yend = .data$yend
          ),
          inherit.aes = FALSE,
          linewidth = 0.7,
          colour = "black"
        )
      }
    }

    legend_position <- legend.settings$position %zresid_or% "right"
    plot <- plot +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        labels = y_labels,
        expand = ggplot2::expansion(mult = 0)
      ) +
      ggplot2::coord_cartesian(
        xlim = x_limits,
        ylim = y_limits,
        clip = "off"
      ) +
      ggplot2::labs(
        title = main.title,
        x = xlab,
        y = ylab,
        colour = NULL,
        linetype = NULL
      ) +
      theme +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        plot.title.position = "panel",
        axis.title = ggplot2::element_text(face = "bold"),
        legend.position = legend_position,
        legend.justification = c(0, 0.45),
        legend.title = ggplot2::element_text(
          size = legend.settings$title_size %zresid_or% ggplot2::rel(1)
        ),
        legend.text = ggplot2::element_text(
          size = legend.settings$text_size %zresid_or% ggplot2::rel(0.9)
        ),
        plot.margin = ggplot2::margin(8, 8, 8, 8)
      )

    if (need_clip) {
      plot <- plot + ggplot2::theme(panel.border = ggplot2::element_blank())
    }

    attr(plot, "zresid_outliers") <- id_outlier
    attr(plot, "zresid_sw_pvalue") <- p_value
    attr(plot, "zresid_clipped") <- need_clip
    attr(plot, "zresid_display_break") <- if (need_clip) display_threshold else NA_real_
    attr(plot, "zresid_qq_data") <- plot_data
    plots[[position]] <- plot

    if (isTRUE(outlier.return) && length(id_outlier)) {
      message("Outlier Indices : ", paste(id_outlier, collapse = ", "))
    }

    if (isTRUE(interactive)) {
      if (!requireNamespace("plotly", quietly = TRUE)) {
        stop(
          "interactive = TRUE requires the optional package 'plotly'.",
          call. = FALSE
        )
      }
      widget <- plotly::ggplotly(plot, tooltip = "text")
      plots[[position]] <- widget
      print(widget)
    } else {
      print(plot)
    }
  }

  result <- list(
    outliers = outliers,
    plots = plots
  )
  return(invisible(result))
}
