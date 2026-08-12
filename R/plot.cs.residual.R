#' Cox-Snell residual plot for survival models
#'
#' Draws a ggplot2 Cox-Snell residual diagnostic based on the Fleming-Harrington
#' estimate of the cumulative hazard. Under a correctly specified model, the
#' curve should be close to the 45-degree reference line.
#'
#' @importFrom survival Surv survfit
#'
#' @param x Numeric Cox-Snell residual vector or one-column matrix carrying a
#'   length-n `censored.status` event-indicator attribute.
#' @param ylab Y-axis label.
#' @param main.title Main plot title.
#' @param outlier.return If `TRUE`, mark and return observations satisfying the
#'   original `abs(x) > 3.5` rule.
#' @param point.args Named arguments for the KM-point `geom_point()` layer.
#' @param curve.args Named arguments for the cumulative-hazard `geom_step()`.
#' @param reference.args Named arguments for the `y = x` reference segment.
#' @param theme A ggplot2 theme.
#' @param ... Additional plotting arguments. Common base arguments `xlab`,
#'   `ylab`, `main`, `xlim`, `ylim`, `col`, `pch`, and `cex` are supported.
#'
#' @return Invisibly returns `NULL`, or `list(outliers = ...)` when
#'   `outlier.return = TRUE`, matching the original method.
#'
#' @method plot cs.residual
#' @export
plot.cs.residual <- function(x,
                             ylab = "Cumulative Hazard Function",
                             main.title = "Cox-Snell Residuals Scatterplot",
                             outlier.return = TRUE,
                             point.args = list(),
                             curve.args = list(),
                             reference.args = list(),
                             theme = ggplot2::theme_bw(),
                             ...) {
  original_attributes <- attributes(x)
  cs_residual <- as.numeric(x)
  n <- length(cs_residual)

  if (!n) {
    stop("plot.cs.residual: x contains no observations.", call. = FALSE)
  }

  status <- original_attributes$censored.status
  if (is.null(status) || length(status) != n) {
    stop(
      "plot.cs.residual: x must contain a length-n 'censored.status' attribute.",
      call. = FALSE
    )
  }
  status <- as.integer(status)
  if (anyNA(status) || any(!(status %in% c(0L, 1L)))) {
    stop("plot.cs.residual: 'censored.status' must contain only 0 and 1.", call. = FALSE)
  }

  id_nonfinite <- which(!is.finite(cs_residual))
  if (length(id_nonfinite)) {
    finite_values <- abs(cs_residual[is.finite(cs_residual)])
    finite_max <- if (length(finite_values)) max(finite_values) else 0
    replacement_sign <- sign(cs_residual[id_nonfinite])
    replacement_sign[!is.finite(replacement_sign) | replacement_sign == 0] <- 1
    cs_residual[id_nonfinite] <- replacement_sign * (finite_max + 0.1)
    message(
      "Non-finite Cox-Snell residuals exist! The model or the fitting process has a problem!"
    )
  }

  is_outlier <- abs(cs_residual) > 3.5
  outlier_indices <- which(is_outlier)

  km <- survival::survfit(
    survival::Surv(cs_residual, status) ~ 1,
    type = "fleming"
  )
  cumulative_hazard <- -log(km$surv)

  ordered_indices <- order(cs_residual)
  ordered_status <- status[ordered_indices]
  if (length(ordered_status) != length(km$time)) {
    time_match <- match(km$time, cs_residual[ordered_indices])
    valid_match <- !is.na(time_match)
    matched_status <- rep(NA_integer_, length(km$time))
    matched_status[valid_match] <- ordered_status[time_match[valid_match]]
    ordered_status <- matched_status
  }

  group <- factor(
    ifelse(ordered_status == 1L, "Uncensored", "Censored"),
    levels = c("Uncensored", "Censored")
  )
  plot_data <- data.frame(
    time = km$time,
    cumulative_hazard = cumulative_hazard,
    group = group
  )

  dots <- list(...)
  if (!is.null(dots$xlab)) {
    x_label <- dots$xlab
  } else {
    x_label <- "Cox-Snell Residuals"
  }
  if (!is.null(dots$ylab)) ylab <- dots$ylab
  if (!is.null(dots$main)) main.title <- dots$main

  old_names <- c(col = "colour", pch = "shape", cex = "size")
  for (old_name in names(old_names)) {
    if (!is.null(dots[[old_name]]) && is.null(dots[[old_names[[old_name]]]])) {
      dots[[old_names[[old_name]]]] <- dots[[old_name]]
    }
  }
  point_allowed <- c("alpha", "colour", "fill", "shape", "size", "stroke", "na.rm")
  dot_point_args <- dots[intersect(names(dots), point_allowed)]
  point_args <- utils::modifyList(list(size = 2, na.rm = TRUE), point.args)
  point_args <- utils::modifyList(point_args, dot_point_args)

  curve_args <- utils::modifyList(
    list(
      data = plot_data,
      mapping = ggplot2::aes(
        x = .data$time,
        y = .data$cumulative_hazard
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.7,
      na.rm = TRUE
    ),
    curve.args
  )

  x_range <- range(plot_data$time, finite = TRUE)
  y_range <- range(plot_data$cumulative_hazard, finite = TRUE)
  if (!all(is.finite(x_range)) || !all(is.finite(y_range))) {
    stop("plot.cs.residual: cumulative hazard could not be computed.", call. = FALSE)
  }
  if (!is.null(dots$xlim) && length(dots$xlim) == 2L) x_range <- dots$xlim
  if (!is.null(dots$ylim) && length(dots$ylim) == 2L) y_range <- dots$ylim

  reference_min <- max(min(x_range), min(y_range))
  reference_max <- min(max(x_range), max(y_range))
  if (!is.finite(reference_min) || !is.finite(reference_max) ||
      reference_max < reference_min) {
    reference_min <- min(x_range)
    reference_max <- max(x_range)
  }
  reference_data <- data.frame(
    x = reference_min,
    xend = reference_max,
    y = reference_min,
    yend = reference_max
  )
  reference_args <- utils::modifyList(
    list(
      data = reference_data,
      mapping = ggplot2::aes(
        x = .data$x,
        xend = .data$xend,
        y = .data$y,
        yend = .data$yend
      ),
      inherit.aes = FALSE,
      colour = "red",
      linetype = "dashed",
      linewidth = 0.8,
      show.legend = FALSE
    ),
    reference.args
  )

  plot <- ggplot2::ggplot()
  plot <- plot + do.call(ggplot2::geom_step, curve_args)
  plot <- plot + do.call(ggplot2::geom_segment, reference_args)
  point_layer_args <- utils::modifyList(
    list(
      data = plot_data,
      mapping = ggplot2::aes(
        x = .data$time,
        y = .data$cumulative_hazard,
        colour = .data$group,
        shape = .data$group
      ),
      inherit.aes = FALSE
    ),
    point_args
  )
  plot <- plot + do.call(ggplot2::geom_point, point_layer_args)

  plot <- plot +
    ggplot2::scale_colour_manual(
      values = c(
        Uncensored = "darkolivegreen4",
        Censored = "blue"
      ),
      breaks = c("Uncensored", "Censored"),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c(Uncensored = 2, Censored = 3),
      breaks = c("Uncensored", "Censored"),
      name = NULL
    )

  if (isTRUE(outlier.return) && length(outlier_indices)) {
    outlier_time <- cs_residual[outlier_indices]
    outlier_position <- match(outlier_time, km$time)
    outlier_data <- data.frame(
      observation = outlier_indices,
      time = outlier_time,
      cumulative_hazard = cumulative_hazard[outlier_position]
    )
    outlier_data <- outlier_data[is.finite(outlier_data$cumulative_hazard), , drop = FALSE]

    if (NROW(outlier_data)) {
      plot <- plot +
        ggplot2::geom_point(
          data = outlier_data,
          mapping = ggplot2::aes(x = .data$time, y = .data$cumulative_hazard),
          inherit.aes = FALSE,
          shape = 1,
          colour = "red",
          size = 4,
          stroke = 1,
          na.rm = TRUE
        ) +
        ggplot2::geom_text(
          data = outlier_data,
          mapping = ggplot2::aes(
            x = .data$time,
            y = .data$cumulative_hazard,
            label = .data$observation
          ),
          inherit.aes = FALSE,
          colour = "red",
          size = 3,
          nudge_y = 0.08,
          vjust = 0,
          check_overlap = TRUE,
          na.rm = TRUE
        )
    }
  }

  x_padding <- max(0.05, 0.04 * diff(x_range))
  y_padding <- max(0.05, 0.06 * diff(y_range))
  x_limits <- x_range + c(-x_padding, x_padding)
  y_limits <- y_range + c(-y_padding, y_padding)

  plot <- plot +
    ggplot2::coord_cartesian(
      xlim = x_limits,
      ylim = y_limits,
      clip = "off"
    ) +
    ggplot2::labs(
      title = main.title,
      x = x_label,
      y = ylab
    ) +
    theme +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.title.position = "panel",
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_rect(fill = "white", colour = "black"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  attr(plot, "cs_residual_outliers") <- outlier_indices
  attr(plot, "cs_residual_data") <- plot_data
  attr(plot, "cs_residual_km") <- km
  print(plot)

  if (isTRUE(outlier.return)) {
    if (!length(outlier_indices)) return(invisible(NULL))
    message("Outlier Indices: ", paste(outlier_indices, collapse = " "))
    return(invisible(list(outliers = outlier_indices)))
  }
  invisible(NULL)
}
