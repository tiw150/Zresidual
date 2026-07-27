#' Predictive quantities for survival::coxph with frailty
#'
#' @param fit A fitted coxph object, usually inheriting from coxph.penal.
#' @param data New data to evaluate.
#' @param traindata Training data used to fit the model. This is required for
#'   reconstructing the baseline cumulative hazard.
#' @param ... Additional arguments passed from Zresidual.
#'
#' @return A list containing log-survival probabilities, log-likelihood
#' contributions, discreteness indicators, and linear predictors.
#'
#' @export
log_pointpred_survival_coxph.penal <- function(
    fit,
    data,
    traindata,
    ...
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop(
      paste0(
        "log_pointpred_survival_coxph.penal requires ",
        "package 'survival'."
      ),
      call. = FALSE
    )
  }
  
  if (
    (missing(data) || is.null(data)) &&
    (missing(traindata) || is.null(traindata))
  ) {
    stop(
      paste0(
        "log_pointpred_survival_coxph.penal: at least one of ",
        "`data` or `traindata` must be provided."
      ),
      call. = FALSE
    )
  }

  if (missing(data) || is.null(data)) {
    data <- traindata
  }

  if (missing(traindata) || is.null(traindata)) {
    traindata <- data
  }
  
  data <- as.data.frame(data)
  traindata <- as.data.frame(traindata)
  
  get_frailty_group_var <- function(object) {
    term_labels <- attr(
      stats::terms(object),
      "term.labels"
    )

    frailty_labels <- term_labels[
      grepl(
        "(^|::)frailty\\s*\\(",
        term_labels
      )
    ]

    if (length(frailty_labels) != 1L) {
      stop(
        paste0(
          "coxph_frailty: expected exactly one ",
          "frailty() term."
        ),
        call. = FALSE
      )
    }

    matched <- regexec(
      "frailty\\s*\\(([^,\\)]+)",
      frailty_labels[1L]
    )

    pieces <- regmatches(
      frailty_labels[1L],
      matched
    )[[1L]]

    if (length(pieces) < 2L) {
      stop(
        paste0(
          "coxph_frailty: cannot parse the ",
          "frailty group variable."
        ),
        call. = FALSE
      )
    }

    group_name <- trimws(pieces[2L])
    gsub("^`|`$", "", group_name)
  }
  
  group_var_name <- get_frailty_group_var(fit)
  
  get_group <- function(dat, name, where) {
    if (!name %in% names(dat)) {
      stop(
        sprintf(
          "coxph_frailty: %s is missing group variable '%s'.",
          where,
          name
        ),
        call. = FALSE
      )
    }

    dat[[name]]
  }
  
  group_train_raw <- get_group(
    traindata,
    group_var_name,
    "traindata"
  )

  group_new_raw <- get_group(
    data,
    group_var_name,
    "data"
  )
  
  frailty_values <- fit$frail

  if (
    is.null(frailty_values) ||
    length(frailty_values) < 1L
  ) {
    stop(
      "coxph_frailty: fitted model contains no frailty estimates.",
      call. = FALSE
    )
  }

  frailty_levels <- NULL

  if (
    !is.null(names(frailty_values)) &&
    all(nzchar(names(frailty_values)))
  ) {
    frailty_levels <- names(frailty_values)
  }
  
  if (is.null(frailty_levels)) {
    group_train <- factor(
      as.character(group_train_raw)
    )
  } else {
    group_train <- factor(
      as.character(group_train_raw),
      levels = frailty_levels
    )
  }
  
  if (anyNA(group_train)) {
    stop(
      paste0(
        "coxph_frailty: `traindata` contains group levels ",
        "not represented in the fitted frailty effects."
      ),
      call. = FALSE
    )
  }
  
  group_new <- factor(
    as.character(group_new_raw),
    levels = levels(group_train)
  )
  
  if (anyNA(group_new)) {
    unseen_groups <- unique(
      as.character(
        group_new_raw[is.na(group_new)]
      )
    )

    stop(
      sprintf(
        paste0(
          "coxph_frailty: `data` contains frailty group ",
          "level(s) absent from `traindata`: %s."
        ),
        paste(unseen_groups, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  beta <- fit$coefficients
  
  if (is.null(beta)) {
    beta <- numeric(0)
  }
  
  if (anyNA(beta) || any(!is.finite(beta))) {
    stop(
      paste0(
        "coxph_frailty: fixed-effect coefficients contain ",
        "missing or non-finite values."
      ),
      call. = FALSE
    )
  }
  
  # Construct a terms object containing fixed effects only.
  #
  # The frailty term must not be evaluated through model.matrix() on
  # a one-row held-out dataset. In LOOCV, that term has only one
  # observed level and therefore cannot define contrasts.
  fixed_terms <- stats::delete.response(
    stats::terms(fit)
  )

  frailty_special <- survival::untangle.specials(
    fixed_terms,
    "frailty"
  )

  if (length(frailty_special$terms) > 0L) {
    fixed_terms <- stats::drop.terms(
      fixed_terms,
      dropx = frailty_special$terms,
      keep.response = FALSE
    )
  }

  make_fixed_matrix <- function(dat, where) {
    fixed_frame <- tryCatch(
      stats::model.frame(
        fixed_terms,
        data = dat,
        na.action = stats::na.pass,
        xlev = fit$xlevels
      ),
      error = function(e) {
        stop(
          sprintf(
            "coxph_frailty: failed to construct %s fixed-effect model frame: %s",
            where,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )
    
    if (any(!stats::complete.cases(fixed_frame))) {
      stop(
        sprintf(
          "coxph_frailty: %s contains missing fixed-effect variables.",
          where
        ),
        call. = FALSE
      )
    }
    
    fixed_matrix <- tryCatch(
      stats::model.matrix(
        fixed_terms,
        data = fixed_frame,
        contrasts.arg = fit$contrasts
      ),
      error = function(e) {
        stop(
          sprintf(
            "coxph_frailty: failed to construct %s fixed-effect model matrix: %s",
            where,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    assign_value <- attr(
      fixed_matrix,
      "assign"
    )

    if (!is.null(assign_value)) {
      fixed_matrix <- fixed_matrix[
        ,
        assign_value != 0L,
        drop = FALSE
      ]
    } else if (
      ncol(fixed_matrix) > 0L &&
      "(Intercept)" %in% colnames(fixed_matrix)
    ) {
      fixed_matrix <- fixed_matrix[
        ,
        colnames(fixed_matrix) != "(Intercept)",
        drop = FALSE
      ]
    }
    
    if (length(beta) == 0L) {
      return(
        matrix(
          numeric(0),
          nrow = nrow(dat),
          ncol = 0L
        )
      )
    }
    
    if (
      is.null(names(beta)) ||
      is.null(colnames(fixed_matrix))
    ) {
      if (ncol(fixed_matrix) != length(beta)) {
        stop(
          sprintf(
            paste0(
              "coxph_frailty: %s fixed-effect matrix has ",
              "%d columns, but the fitted model has %d ",
              "fixed-effect coefficients."
            ),
            where,
            ncol(fixed_matrix),
            length(beta)
          ),
          call. = FALSE
        )
      }

      return(fixed_matrix)
    }
    
    coefficient_index <- match(
      names(beta),
      colnames(fixed_matrix)
    )
    
    if (anyNA(coefficient_index)) {
      missing_columns <- names(beta)[
        is.na(coefficient_index)
      ]

      stop(
        sprintf(
          paste0(
            "coxph_frailty: %s fixed-effect matrix is missing ",
            "coefficient column(s): %s."
          ),
          where,
          paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    fixed_matrix[
      ,
      coefficient_index,
      drop = FALSE
    ]
  }
  
  fixed_train <- make_fixed_matrix(
    traindata,
    "training"
  )

  fixed_new <- make_fixed_matrix(
    data,
    "new-data"
  )

  fixed_lp <- function(matrix_value, coefficients) {
    if (length(coefficients) == 0L) {
      return(
        rep(
          0,
          nrow(matrix_value)
        )
      )
    }

    as.vector(
      matrix_value %*% coefficients
    )
  }
  
  lp_fixed_train <- fixed_lp(
    fixed_train,
    beta
  )

  lp_fixed_new <- fixed_lp(
    fixed_new,
    beta
  )

  frailty_train <- as.numeric(
    frailty_values[
      as.integer(group_train)
    ]
  )

  frailty_new <- as.numeric(
    frailty_values[
      as.integer(group_new)
    ]
  )

  frailty_train[is.na(frailty_train)] <- 0
  frailty_new[is.na(frailty_new)] <- 0

  linear_pred_train <- (
    lp_fixed_train +
      frailty_train
  )

  linear_pred_new <- (
    lp_fixed_new +
      frailty_new
  )

  relative_risk_train <- exp(
    linear_pred_train
  )
  
  response_formula <- stats::update.formula(
    stats::formula(fit),
    . ~ 1
  )

  make_surv_response <- function(dat, where) {
    response_frame <- tryCatch(
      stats::model.frame(
        response_formula,
        data = dat,
        na.action = stats::na.pass
      ),
      error = function(e) {
        stop(
          sprintf(
            "coxph_frailty: failed to construct %s response: %s",
            where,
            conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    response <- stats::model.response(
      response_frame
    )

    if (!inherits(response, "Surv")) {
      stop(
        sprintf(
          "coxph_frailty: %s response is not a Surv object.",
          where
        ),
        call. = FALSE
      )
    }

    response
  }
  
  as_counting_surv <- function(response) {
    response_matrix <- as.matrix(response)

    if (ncol(response_matrix) == 2L) {
      return(
        as.matrix(
          survival::Surv(
            rep(0, nrow(response_matrix)),
            response_matrix[, 1L],
            response_matrix[, 2L]
          )
        )
      )
    }

    if (ncol(response_matrix) == 3L) {
      return(response_matrix)
    }

    stop(
      "coxph_frailty: unsupported Surv response structure.",
      call. = FALSE
    )
  }

  response_train <- make_surv_response(
    traindata,
    "training"
  )

  response_new <- make_surv_response(
    data,
    "new-data"
  )

  y_train <- as_counting_surv(
    response_train
  )

  y_new <- as_counting_surv(
    response_new
  )

  get_cumulative_hazard <- function(y, relative_risk) {
    death <- (
      y[, ncol(y)] == 1
    )

    stop_time <- y[, ncol(y) - 1L]

    event_times <- sort(
      unique(stop_time)
    )

    if (length(event_times) < 1L) {
      return(
        data.frame(
          time = numeric(0),
          cumhaz = numeric(0)
        )
      )
    }

    number_events <- as.vector(
      rowsum(
        as.integer(death),
        stop_time
      )
    )

    number_at_risk <- rev(
      cumsum(
        rev(
          rowsum(
            relative_risk,
            stop_time
          )
        )
      )
    )

    time_increment <- if (length(event_times) >= 2L) {
      min(diff(event_times)) / 2
    } else {
      0.5
    }

    entry_times <- c(
      sort(unique(y[, 1L])),
      max(y[, 1L]) + time_increment
    )

    entry_index <- stats::approx(
      entry_times,
      seq_along(entry_times),
      event_times,
      method = "constant",
      rule = 2,
      f = 1
    )$y

    entry_risk_sum <- rev(
      cumsum(
        rev(
          rowsum(
            relative_risk,
            y[, 1L]
          )
        )
      )
    )

    number_at_risk <- (
      number_at_risk -
        c(entry_risk_sum, 0)[entry_index]
    )

    number_at_risk <- pmax(
      number_at_risk,
      .Machine$double.xmin
    )

    cumulative_hazard <- cumsum(
      number_events / number_at_risk
    )

    data.frame(
      time = event_times,
      cumhaz = cumulative_hazard
    )
  }
  
  baseline_hazard <- get_cumulative_hazard(
    y_train,
    relative_risk_train
  )
  
  if (nrow(baseline_hazard) < 1L) {
    stop(
      paste0(
        "coxph_frailty: failed to construct the baseline ",
        "cumulative hazard."
      ),
      call. = FALSE
    )
  }
  
  stop_time_new <- y_new[, 2L]
  status_new <- as.integer(
    y_new[, 3L]
  )

  n_new <- nrow(y_new)

  risk_multiplier_new <- exp(
    linear_pred_new
  )

  index_after <- findInterval(
    stop_time_new,
    baseline_hazard$time
  )

  hazard_after <- numeric(n_new)
  
  positive_after <- which(
    index_after > 0L
  )

  if (length(positive_after) > 0L) {
    hazard_after[positive_after] <- baseline_hazard$cumhaz[
      index_after[positive_after]
    ]
  }
  
  is_event <- (
    status_new == 1L
  )

  tolerance <- 1e-7
  exact_event <- rep(
    FALSE,
    n_new
  )

  event_on_grid <- which(
    is_event &
      index_after > 0L
  )

  if (length(event_on_grid) > 0L) {
    exact_event[event_on_grid] <- (
      abs(
        baseline_hazard$time[
          index_after[event_on_grid]
        ] -
          stop_time_new[event_on_grid]
      ) <= tolerance
    )
  }
  
  index_before <- index_after
  index_before[exact_event] <- (
    index_before[exact_event] - 1L
  )
  
  hazard_before <- numeric(n_new)
  
  positive_before <- which(
    index_before > 0L
  )
  
  if (length(positive_before) > 0L) {
    hazard_before[positive_before] <- baseline_hazard$cumhaz[
      index_before[positive_before]
    ]
  }
  
  log_survival_after <- (
    -risk_multiplier_new *
      hazard_after
  )
  
  log_survival_before <- (
    -risk_multiplier_new *
      hazard_before
  )
  
  log_surv <- log_survival_after
  log_surv[!is_event] <- -Inf
  
  log_like <- log_survival_after

  if (any(is_event)) {
    log_like[is_event] <- log_diff_exp(
      log_survival_before[is_event],
      log_survival_after[is_event]
    )
  }
  
  is_discrete <- as.integer(
    !is_event
  )
  
  list(
    log_surv = matrix(
      log_surv,
      nrow = 1L,
      ncol = n_new
    ),
    log_like = matrix(
      log_like,
      nrow = 1L,
      ncol = n_new
    ),
    is_discrete = matrix(
      is_discrete,
      nrow = 1L,
      ncol = n_new
    ),
    linear_pred = matrix(
      linear_pred_new,
      nrow = 1L,
      ncol = n_new
    )
  )
}