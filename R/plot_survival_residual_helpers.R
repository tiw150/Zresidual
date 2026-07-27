utils::globalVariables(".data")

.plot_survival_residual_ggplot <- function(
    x,
    residual_name,
    ylab,
    x_axis_var,
    main_title,
    outlier_return,
    is_outlier,
    point_args = list(),
    smooth_args = list(),
    theme = ggplot2::theme_bw(),
    dots = list()) {
  residual_attributes <- attributes(x)
  residual <- as.numeric(x)
  n <- length(residual)

  if (!n) {
    stop("The residual object contains no observations.", call. = FALSE)
  }

  status <- residual_attributes$censored.status
  if (is.null(status) || length(status) != n) {
    stop(
      "The residual object must contain a length-n 'censored.status' attribute.",
      call. = FALSE
    )
  }
  status <- as.integer(status)
  if (anyNA(status) || any(!(status %in% c(0L, 1L)))) {
    stop("'censored.status' must contain only 0 and 1.", call. = FALSE)
  }

  id_nonfinite <- which(!is.finite(residual))
  if (length(id_nonfinite)) {
    finite_values <- abs(residual[is.finite(residual)])
    finite_max <- if (length(finite_values)) max(finite_values) else 0
    replacement_sign <- sign(residual[id_nonfinite])
    replacement_sign[!is.finite(replacement_sign) | replacement_sign == 0] <- 1
    residual[id_nonfinite] <- replacement_sign * (finite_max + 0.1)
    message(
      "Non-finite ", tolower(residual_name),
      " residual exist! The model or the fitting process has a problem!"
    )
  }

  X <- x_axis_var
  if (is.character(X) && length(X) > 1L) X <- X[1L]
  if (!is.character(X) || length(X) != 1L) {
    stop(
      "x_axis_var must be 'index', 'lp', 'covariate', or a covariate name.",
      call. = FALSE
    )
  }

  add_lowess <- FALSE
  if (identical(X, "index")) {
    x_value <- seq_len(n)
    x_label <- "Index"
  } else if (identical(X, "lp")) {
    x_value <- residual_attributes$linear.pred
    if (is.null(x_value)) x_value <- residual_attributes$linear_pred
    if (is.null(x_value) || length(x_value) != n) {
      stop("The residual object is missing a length-n 'linear.pred' attribute.", call. = FALSE)
    }
    x_label <- "Linear Predictor"
    add_lowess <- TRUE
  } else {
    covariates <- residual_attributes$covariates
    if (is.null(covariates)) {
      stop("The residual object is missing the 'covariates' attribute.", call. = FALSE)
    }
    if (!is.data.frame(covariates)) covariates <- as.data.frame(covariates)
    if (NROW(covariates) != n) {
      stop("The covariate rows do not match the residual length.", call. = FALSE)
    }
    covariate_names <- names(covariates)
    if (identical(X, "covariate")) {
      selected <- covariate_names[1L]
      cat(
        "To plot against other covariates, set x_axis_var to the covariate name. Available covariates: ",
        paste(covariate_names, collapse = ", "), "\n",
        sep = ""
      )
    } else {
      selected <- X
      if (!(selected %in% covariate_names)) {
        stop(
          "x_axis_var must be one of the covariate names: ",
          paste(covariate_names, collapse = ", "), ".",
          call. = FALSE
        )
      }
    }
    x_value <- covariates[[selected]]
    x_label <- selected
    add_lowess <- TRUE
  }

  if (is.factor(x_value)) {
    if (add_lowess) {
      warning("LOWESS is not drawn for a factor x-axis.", call. = FALSE)
      add_lowess <- FALSE
    }
  } else if (!is.numeric(x_value)) {
    numeric_x <- suppressWarnings(as.numeric(x_value))
    if (anyNA(numeric_x)) {
      stop("The selected x-axis variable must be numeric or a factor.", call. = FALSE)
    }
    x_value <- numeric_x
  }

  group <- factor(
    ifelse(status == 1L, "Uncensored", "Censored"),
    levels = c("Uncensored", "Censored")
  )
  plot_data <- data.frame(
    observation = seq_len(n),
    x_value = x_value,
    residual = residual,
    group = group,
    stringsAsFactors = FALSE
  )

  default_point_args <- list(size = 2, na.rm = TRUE)
  old_names <- c(col = "colour", pch = "shape", cex = "size")
  for (old_name in names(old_names)) {
    if (!is.null(dots[[old_name]]) && is.null(dots[[old_names[[old_name]]]])) {
      dots[[old_names[[old_name]]]] <- dots[[old_name]]
    }
  }
  point_allowed <- c("alpha", "colour", "fill", "shape", "size", "stroke", "na.rm")
  dot_point_args <- dots[intersect(names(dots), point_allowed)]
  point_args <- utils::modifyList(default_point_args, point_args)
  point_args <- utils::modifyList(point_args, dot_point_args)

  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$x_value,
      y = .data$residual,
      colour = .data$group,
      shape = .data$group
    )
  )
  plot <- plot + do.call(ggplot2::geom_point, point_args)
  plot <- plot +
    ggplot2::scale_colour_manual(
      values = c(Uncensored = "darkolivegreen4", Censored = "blue"),
      breaks = c("Uncensored", "Censored"),
      name = NULL
    ) +
    ggplot2::scale_shape_manual(
      values = c(Uncensored = 2, Censored = 3),
      breaks = c("Uncensored", "Censored"),
      name = NULL
    )

  if (isTRUE(add_lowess)) {
    numeric_x <- as.numeric(x_value)
    keep <- is.finite(numeric_x) & is.finite(residual)
    if (sum(keep) >= 3L) {
      lowess_result <- stats::lowess(
        x = numeric_x[keep],
        y = residual[keep]
      )
      smooth_data <- data.frame(
        x_value = lowess_result$x,
        residual = lowess_result$y
      )
      default_smooth_args <- list(
        data = smooth_data,
        mapping = ggplot2::aes(x = .data$x_value, y = .data$residual),
        inherit.aes = FALSE,
        colour = "red",
        linewidth = 1.2,
        na.rm = TRUE
      )
      smooth_args <- utils::modifyList(default_smooth_args, smooth_args)
      plot <- plot + do.call(ggplot2::geom_line, smooth_args)
    }
  }

  outlier_indices <- integer(0L)
  if (isTRUE(outlier_return)) {
    if (is.null(is_outlier)) {
      stop(
        "outlier.return = TRUE requires a logical vector named 'is.outlier' in the calling environment.",
        call. = FALSE
      )
    }
    if (!is.logical(is_outlier) || length(is_outlier) != n || anyNA(is_outlier)) {
      stop("'is.outlier' must be a non-missing logical vector matching x.", call. = FALSE)
    }
    outlier_indices <- which(is_outlier)
    if (length(outlier_indices)) {
      outlier_data <- plot_data[outlier_indices, , drop = FALSE]
      plot <- plot +
        ggplot2::geom_point(
          data = outlier_data,
          mapping = ggplot2::aes(x = .data$x_value, y = .data$residual),
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
            x = .data$x_value,
            y = .data$residual,
            label = .data$observation
          ),
          inherit.aes = FALSE,
          colour = "red",
          size = 3,
          nudge_y = -0.12,
          vjust = 1,
          check_overlap = TRUE,
          na.rm = TRUE
        )
    }
  }

  finite_residual <- residual[is.finite(residual)]
  y_limits <- c(min(finite_residual), max(finite_residual) + 1)
  if (!is.null(dots$ylim) && length(dots$ylim) == 2L) y_limits <- dots$ylim
  if (!is.null(dots$xlab)) x_label <- dots$xlab
  if (!is.null(dots$ylab)) ylab <- dots$ylab
  if (!is.null(dots$main)) main_title <- dots$main

  plot <- plot +
    ggplot2::coord_cartesian(ylim = y_limits, clip = "off") +
    ggplot2::labs(
      title = main_title,
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

  attr(plot, "survival_residual_outliers") <- outlier_indices
  attr(plot, "survival_residual_data") <- plot_data
  attr(plot, "survival_residual_type") <- residual_name
  plot
}
