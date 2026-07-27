utils::globalVariables(".data")

`%zresid_or%` <- function(x, y) {
  if (!is.null(x)) x else y
}

.zresid_as_matrix <- function(x) {
  old_attributes <- attributes(x)
  out <- if (is.null(dim(x))) matrix(x, ncol = 1L) else as.matrix(x)

  keep <- setdiff(
    names(old_attributes),
    c("dim", "dimnames", "names", "class")
  )
  for (name in keep) {
    attr(out, name) <- old_attributes[[name]]
  }
  out
}

.zresid_first_nonempty <- function(...) {
  values <- list(...)
  for (value in values) {
    if (!is.null(value) && length(value) > 0L) return(value)
  }
  NULL
}

.zresid_get_zcov <- function(x, zcov) {
  if (!is.null(zcov)) return(zcov)
  attr(x, "zcov") %zresid_or% attr(x, "Zcov") %zresid_or% list()
}

.zresid_type_from_zcov <- function(zcov) {
  if (!is.list(zcov)) return(NULL)

  if (!is.null(zcov$type) && length(zcov$type) > 0L) {
    type <- as.character(zcov$type)[1L]
    if (!is.na(type) && nzchar(type)) return(type)
  }

  kind <- if (length(zcov$y_type_kind)) {
    as.character(zcov$y_type_kind)[1L]
  } else {
    NULL
  }

  switch(
    kind,
    censor = "survival",
    trunc = "count",
    hurdle = "hurdle",
    NULL
  )
}

.zresid_fill_metadata <- function(x, zcov) {
  if (!is.list(zcov)) return(x)

  covariates <- zcov$covariates %zresid_or% zcov$covariate
  linear_predictor <- zcov$linear_pred %zresid_or% zcov$linear.pred

  if (is.null(attr(x, "covariates")) && !is.null(covariates)) {
    attr(x, "covariates") <- covariates
  }
  if (is.null(attr(x, "linear_pred")) && !is.null(linear_predictor)) {
    attr(x, "linear_pred") <- linear_predictor
  }
  if (is.null(attr(x, "linear.pred")) && !is.null(linear_predictor)) {
    attr(x, "linear.pred") <- linear_predictor
  }
  if (is.null(attr(x, "type"))) {
    attr(x, "type") <- .zresid_type_from_zcov(zcov)
  }
  if (is.null(attr(x, "zero_id"))) {
    zero_id <- NULL
    if (is.list(zcov$extra)) zero_id <- zcov$extra$zero_id
    zero_id <- zero_id %zresid_or% zcov$zero_id
    if (!is.null(zero_id)) attr(x, "zero_id") <- zero_id
  }
  if (is.null(attr(x, "y_type")) && !is.null(zcov$y_type)) {
    attr(x, "y_type") <- zcov$y_type
  }
  if (is.null(attr(x, "censored.status"))) {
    censored <- zcov$censored.status
    if (is.null(censored) && identical(zcov$y_type_kind, "censor")) {
      if (!is.null(zcov$y_type)) censored <- as.integer(zcov$y_type == 0L)
    }
    if (!is.null(censored)) attr(x, "censored.status") <- censored
  }

  x
}

.zresid_resolve_label <- function(label) {
  if (is.null(label)) return(NULL)
  label_character <- as.character(label)[1L]

  if (!startsWith(label_character, "tex(")) return(label)

  latex_string <- sub("tex\\((.*)\\)", "\\1", label_character)
  if (requireNamespace("latex2exp", quietly = TRUE)) {
    return(latex2exp::TeX(latex_string))
  }

  warning(
    "The 'latex2exp' package is not installed; xlab is shown as plain text.",
    call. = FALSE
  )
  latex_string
}

.zresid_resolve_x <- function(x, x_axis_var, zcov, expression_label) {
  n <- NROW(x)
  keywords <- c("index", "lp", "covariate")

  resolve_again <- function(value, function_label = NULL) {
    if (is.function(value)) {
      result <- value(x, zcov)
      if (is.list(result) && !is.null(result$values)) {
        return(resolve_again(result$values, result$label))
      }
      return(resolve_again(result, function_label))
    }

    is_user_vector <- length(value) == n &&
      !(is.character(value) && length(value) == 1L && value %in% keywords)

    if (is_user_vector) {
      label <- function_label
      if (is.null(label)) {
        label <- paste(deparse(expression_label), collapse = "")
        if (identical(label, "x_axis_var")) label <- "X"
      }
      if (is.logical(value)) value <- as.integer(value)
      return(list(values = value, label = label, test_x = value, mode = "user_vector"))
    }

    if (!is.character(value) || length(value) != 1L) {
      stop(
        "x_axis_var must be 'index', 'lp', 'covariate', a covariate name, a length-n vector, or function(z, zcov).",
        call. = FALSE
      )
    }

    if (identical(value, "index")) {
      return(list(values = seq_len(n), label = "Index", test_x = "index", mode = "index"))
    }

    if (identical(value, "lp")) {
      linear_predictor <- .zresid_first_nonempty(
        zcov$linear_pred,
        zcov$linear.pred,
        attr(x, "linear_pred"),
        attr(x, "linear.pred")
      )
      if (is.null(linear_predictor)) {
        stop("plot.zresid: x_axis_var = 'lp' requires linear-predictor metadata.", call. = FALSE)
      }
      if (length(linear_predictor) != n) {
        stop("plot.zresid: linear-predictor length does not match the residual rows.", call. = FALSE)
      }
      return(list(
        values = as.numeric(linear_predictor),
        label = "Linear Predictor",
        test_x = "lp",
        mode = "lp"
      ))
    }

    covariates <- .zresid_first_nonempty(
      zcov$covariates,
      zcov$covariate,
      attr(x, "covariates")
    )
    if (is.null(covariates)) {
      stop("plot.zresid: covariates are missing; supply zcov or a length-n x-axis vector.", call. = FALSE)
    }
    if (!is.data.frame(covariates)) covariates <- as.data.frame(covariates)
    if (NROW(covariates) != n) {
      stop("plot.zresid: covariate rows do not match the residual rows.", call. = FALSE)
    }

    covariate_names <- names(covariates)
    if (!length(covariate_names)) stop("plot.zresid: covariates have no names.", call. = FALSE)

    selected <- value
    if (identical(value, "covariate")) {
      selected <- covariate_names[1L]
      message(
        "x_axis_var = 'covariate' uses the first covariate: ", selected,
        ". Available covariates: ", paste(covariate_names, collapse = ", ")
      )
    }
    if (!(selected %in% covariate_names)) {
      stop(
        "plot.zresid: covariate '", selected, "' was not found. Available covariates: ",
        paste(covariate_names, collapse = ", "),
        call. = FALSE
      )
    }

    values <- covariates[[selected]]
    if (is.logical(values)) values <- as.integer(values)
    list(values = values, label = selected, test_x = selected, mode = "covariate")
  }

  resolve_again(x_axis_var)
}

.zresid_make_group <- function(x, zcov, category, category_label) {
  n <- NROW(x)
  colours <- NULL
  shapes <- NULL
  legend_title <- NULL

  if (!is.null(category)) {
    if (length(category) != n) {
      stop("plot.zresid: length(category) must match the residual rows.", call. = FALSE)
    }
    group <- factor(category, levels = unique(category))
    levels_group <- levels(group)
    colours <- stats::setNames(
      grDevices::hcl.colors(max(2L, length(levels_group)), "Dark 3")[seq_along(levels_group)],
      levels_group
    )
    shapes <- stats::setNames(rep(1:25, length.out = length(levels_group)), levels_group)
    legend_title <- category_label
    return(list(group = group, colours = colours, shapes = shapes, title = legend_title, show = TRUE))
  }

  type <- attr(x, "type")
  type <- if (length(type)) as.character(type)[1L] else NULL

  if (identical(type, "survival")) {
    censored <- attr(x, "censored.status")
    if (!is.null(censored) && length(censored) == n) {
      group <- factor(ifelse(as.integer(censored) == 1L, "Censored", "Uncensored"),
                      levels = c("Uncensored", "Censored"))
      return(list(
        group = group,
        colours = c(Uncensored = "blue", Censored = "red"),
        shapes = c(Uncensored = 3, Censored = 2),
        title = NULL,
        show = TRUE
      ))
    }
  }

  if (identical(type, "hurdle")) {
    zero_id <- attr(x, "zero_id") %zresid_or% integer(0L)
    group <- factor(ifelse(seq_len(n) %in% zero_id, "zero", "count"),
                    levels = c("count", "zero"))
    return(list(
      group = group,
      colours = c(count = "blue", zero = "red"),
      shapes = c(count = 3, zero = 2),
      title = NULL,
      show = TRUE
    ))
  }

  if (!is.null(type) && type %in% c("zero", "logistic")) {
    y_type <- .zresid_first_nonempty(zcov$y_type, attr(x, "y_type"))
    if (is.null(y_type) || length(y_type) != n) {
      stop("plot.zresid: zero/logistic plotting requires y_type metadata.", call. = FALSE)
    }
    family_name <- if (length(zcov$family)) as.character(zcov$family)[1L] else ""
    if (grepl("^hurdle(_|$)", family_name)) {
      group <- factor(ifelse(as.integer(y_type) == 0L, "zero", "count"),
                      levels = c("count", "zero"))
      colours <- c(count = "blue", zero = "red")
      shapes <- c(count = 3, zero = 2)
    } else {
      group <- factor(as.character(as.integer(y_type)), levels = c("0", "1"))
      colours <- c("0" = "blue", "1" = "red")
      shapes <- c("0" = 3, "1" = 2)
    }
    return(list(group = group, colours = colours, shapes = shapes, title = NULL, show = TRUE))
  }

  if (identical(type, "count")) {
    return(list(
      group = factor(rep("count", n)),
      colours = c(count = "blue"),
      shapes = c(count = 3),
      title = NULL,
      show = TRUE
    ))
  }

  if (!is.null(type) && nzchar(type)) {
    group <- factor(rep(type, n))
    return(list(
      group = group,
      colours = stats::setNames("red", type),
      shapes = stats::setNames(2, type),
      title = NULL,
      show = TRUE
    ))
  }

  list(
    group = factor(rep("Observation", n)),
    colours = c(Observation = "black"),
    shapes = c(Observation = 1),
    title = NULL,
    show = FALSE
  )
}

.zresid_sanitize_column <- function(z, outlier_value) {
  original <- z
  infinite <- is.infinite(original)
  nan <- is.nan(original)
  finite_values <- abs(original[is.finite(original)])
  finite_max <- if (length(finite_values)) max(finite_values) else stats::qnorm(0.9999)

  display <- original
  display[infinite] <- sign(original[infinite]) * (finite_max + 0.1)
  limit <- max(stats::qnorm(0.9999), abs(display[is.finite(display)]), na.rm = TRUE)
  if (!is.finite(limit)) {
    stop("plot.zresid: residuals contain no values that can be plotted.", call. = FALSE)
  }

  list(
    original = original,
    display = display,
    infinite = infinite,
    nan = nan,
    outlier = abs(display) > outlier_value | infinite,
    infinite_label = ifelse(infinite, ifelse(original > 0, "Inf", "-Inf"), NA_character_),
    limit = limit
  )
}

.zresid_diagnostic_table <- function(x, irep, x_info, zcov, tests, k_test) {
  tests <- unique(toupper(tests))
  tests <- tests[tests %in% c("SW", "AOV", "BL")]
  if (!length(tests) || identical(x_info$mode, "index")) return(data.frame())

  result <- vector("list", length(tests))
  for (i in seq_along(tests)) {
    test <- tests[i]
    values <- tryCatch(
      switch(
        test,
        SW = sw.test.zresid(x),
        AOV = aov.test.zresid(x, X = x_info$test_x, zcov = zcov, k.anova = k_test),
        BL = bartlett.test.zresid(x, X = x_info$test_x, zcov = zcov, k.bl = k_test)
      ),
      error = function(e) rep(NA_real_, NCOL(x))
    )
    result[[i]] <- data.frame(
      replication_index = irep,
      test = test,
      p_value = suppressWarnings(as.numeric(values[irep])),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, result)
}

.zresid_prepare_plot <- function(x, zcov, irep, x_axis_var, x_axis_expression,
                                  category, category_label, outlier_value,
                                  normality_test, k_test, add_lowess, x_scale) {
  residuals <- .zresid_as_matrix(x)
  if (NROW(residuals) < 1L || NCOL(residuals) < 1L) {
    stop("plot.zresid: x must contain at least one row and one residual column.", call. = FALSE)
  }

  if (!is.numeric(irep) || anyNA(irep) || any(irep != as.integer(irep))) {
    stop("plot.zresid: irep must contain integer column indices.", call. = FALSE)
  }
  irep <- as.integer(irep)
  if (any(irep < 1L | irep > NCOL(residuals))) {
    stop("plot.zresid: irep contains an invalid residual column index.", call. = FALSE)
  }

  metadata <- .zresid_get_zcov(x, zcov)
  residuals <- .zresid_fill_metadata(residuals, metadata)
  x_info <- .zresid_resolve_x(residuals, x_axis_var, metadata, x_axis_expression)
  if (length(x_info$values) != NROW(residuals)) {
    stop("plot.zresid: resolved x-axis length does not match the residual rows.", call. = FALSE)
  }
  groups <- .zresid_make_group(residuals, metadata, category, category_label)

  plot_parts <- vector("list", length(irep))
  smooth_parts <- vector("list", length(irep))
  outlier_list <- stats::setNames(vector("list", length(irep)), paste0("Replicate ", irep))
  limits <- numeric(length(irep))

  for (position in seq_along(irep)) {
    column <- irep[position]
    cleaned <- .zresid_sanitize_column(residuals[, column], outlier_value)
    replication <- factor(
      rep(paste0("Replicate ", column), NROW(residuals)),
      levels = paste0("Replicate ", irep)
    )

    plot_parts[[position]] <- data.frame(
      observation = seq_len(NROW(residuals)),
      replication = replication,
      replication_index = column,
      x_value = x_info$values,
      z_original = cleaned$original,
      z_display = cleaned$display,
      group = groups$group,
      is_outlier = cleaned$outlier,
      is_infinite = cleaned$infinite,
      is_nan = cleaned$nan,
      infinite_label = cleaned$infinite_label,
      tooltip = paste0(
        "Observation: ", seq_len(NROW(residuals)),
        "<br>Replicate: ", column,
        "<br>x: ", as.character(x_info$values),
        "<br>Z: ", as.character(cleaned$original),
        "<br>Group: ", as.character(groups$group)
      ),
      stringsAsFactors = FALSE
    )
    outlier_list[[position]] <- which(cleaned$outlier)
    limits[position] <- cleaned$limit

    if (isTRUE(add_lowess) &&
        (is.numeric(x_info$values) || inherits(x_info$values, c("Date", "POSIXt")))) {
      x_numeric <- as.numeric(x_info$values)
      ok <- is.finite(x_numeric) & is.finite(cleaned$display)
      if (identical(x_scale, "log10")) ok <- ok & x_numeric > 0
      if (sum(ok) >= 3L) {
        smooth_x <- if (identical(x_scale, "log10")) log10(x_numeric[ok]) else x_numeric[ok]
        lowess_result <- stats::lowess(smooth_x, cleaned$display[ok])
        plotted_x <- if (identical(x_scale, "log10")) 10^lowess_result$x else lowess_result$x
        if (inherits(x_info$values, "Date")) plotted_x <- as.Date(plotted_x, origin = "1970-01-01")
        smooth_parts[[position]] <- data.frame(
          replication = factor(paste0("Replicate ", column), levels = paste0("Replicate ", irep)),
          x_value = plotted_x,
          z_display = lowess_result$y
        )
      }
    }
  }

  diagnostics <- .zresid_diagnostic_table(
    residuals, irep, x_info, metadata, normality_test, k_test
  )
  annotation <- data.frame()
  if (NROW(diagnostics)) {
    labels <- lapply(split(diagnostics, diagnostics$replication_index), function(value) {
      formatted <- ifelse(is.finite(value$p_value), sprintf("%.2f", value$p_value), "NA")
      paste(
        c("P-value:", paste0(value$test, " - ", formatted)),
        collapse = "\n"
      )
    })
    annotation <- data.frame(
      replication = factor(
        paste0("Replicate ", as.integer(names(labels))),
        levels = paste0("Replicate ", irep)
      ),
      label = unlist(labels, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  }

  list(
    data = do.call(rbind, plot_parts),
    smooth = if (any(lengths(smooth_parts))) do.call(rbind, smooth_parts[lengths(smooth_parts) > 0L]) else data.frame(),
    diagnostics = diagnostics,
    annotation = annotation,
    outliers = outlier_list,
    x_info = x_info,
    groups = groups,
    y_limit = max(limits),
    residuals = residuals
  )
}
