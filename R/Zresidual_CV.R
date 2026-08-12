#' Cross-validated Z-residuals
#'
#' Compute cross-validated Z-residuals for supported frequentist survival
#' models. Each fold is evaluated by refitting the model on its training data
#' and calling \code{Zresidual()} on the corresponding held-out data.
#'
#' Currently supported fitted-model classes are \code{"coxph"}, including
#' \code{"coxph.penal"} frailty models, and \code{"survreg"}.
#'
#' @param object A fitted survival model. Supported classes are \code{"coxph"}
#'   and \code{"survreg"}.
#' @param data A data frame containing all variables used to fit \code{object}.
#' @param nfolds Number of cross-validation folds. When both \code{nfolds} and
#'   \code{foldlist} are \code{NULL}, the backend uses
#'   \code{min(10, nrow(complete_case_data))}.
#' @param foldlist Optional list of held-out row indices. Indices refer to the
#'   complete-case data used internally. Supply either \code{nfolds} or
#'   \code{foldlist}, not both.
#' @param nrep Number of randomized Z-residual replicates for each observation.
#' @param randomized Logical. If \code{TRUE}, randomized residuals are generated
#'   for censored observations. If \code{FALSE}, one non-randomized residual is
#'   returned regardless of \code{nrep}.
#' @param type Optional model subtype passed to the model-specific backend.
#' @param log_pointpred Optional predictive backend function or function name.
#'   When \code{NULL}, the backend is resolved automatically.
#' @param seed Optional single integer seed used for fold generation and
#'   residual randomization.
#' @param ... Additional arguments passed to the model-specific CV backend and
#'   ultimately to \code{Zresidual()}.
#'
#' @return A matrix-like object with class
#'   \code{c("cvzresid", "zresid", "matrix", "array")}. Rows correspond to
#'   complete observations and columns correspond to residual replicates.
#'   Attributes include fold indices, original row indices, probability-scale
#'   residuals, covariates, linear predictors, and a \code{zcov} metadata list.
#'
#' @details
#' \code{Zresidual_CV()} is the public dispatcher. Cox calculations are handled
#' by \code{Zresidual_CV_survival_coxph()}, and parametric survival calculations
#' are handled by \code{Zresidual_CV_survival_survreg()}.
#'
#' For frailty Cox models, the current CV target is within-group prediction:
#' every frailty group appearing in a test fold must also appear in that fold's
#' training data.
#'
#' @examples
#' lung2 <- stats::na.omit(
#'   survival::lung[, c("time", "status", "age", "sex", "ph.ecog")]
#' )
#'
#' fit <- survival::survreg(
#'   survival::Surv(time, status == 2) ~ age + sex + ph.ecog,
#'   data = lung2,
#'   dist = "weibull"
#' )
#'
#' zcv <- Zresidual_CV(
#'   object = fit,
#'   data = lung2,
#'   nfolds = 2,
#'   nrep = 2,
#'   seed = 123
#' )
#'
#' print(zcv)
#'
#' @export
Zresidual_CV <- function(object,
                         data,
                         nfolds = NULL,
                         foldlist = NULL,
                         nrep = 1L,
                         randomized = TRUE,
                         type = NULL,
                         log_pointpred = NULL,
                         seed = NULL,
                         ...) {
  if (missing(object) || is.null(object)) {
    stop("Zresidual_CV: `object` must be provided.", call. = FALSE)
  }

  if (missing(data) || is.null(data)) {
    stop("Zresidual_CV: `data` must be provided.", call. = FALSE)
  }

  if (!is.data.frame(data)) {
    data <- tryCatch(
      as.data.frame(data),
      error = function(e) {
        stop("Zresidual_CV: `data` must be coercible to a data frame.", call. = FALSE)
      }
    )
  }

  if (nrow(data) < 2L) {
    stop("Zresidual_CV: `data` must contain at least two observations.", call. = FALSE)
  }

  nrep_value <- suppressWarnings(as.numeric(nrep))
  if (length(nrep_value) != 1L || !is.finite(nrep_value) ||
      nrep_value < 1 || nrep_value != floor(nrep_value)) {
    stop("Zresidual_CV: `nrep` must be a positive integer.", call. = FALSE)
  }
  nrep <- as.integer(nrep_value)

  if (length(randomized) != 1L || is.na(randomized) || !is.logical(randomized)) {
    stop("Zresidual_CV: `randomized` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(nfolds)) {
    nfolds_value <- suppressWarnings(as.numeric(nfolds))
    if (length(nfolds_value) != 1L || !is.finite(nfolds_value) ||
        nfolds_value < 2 || nfolds_value != floor(nfolds_value)) {
      stop("Zresidual_CV: `nfolds` must be an integer >= 2.", call. = FALSE)
    }
    nfolds <- as.integer(nfolds_value)
    if (nfolds > nrow(data)) {
      stop("Zresidual_CV: `nfolds` cannot exceed `nrow(data)`.", call. = FALSE)
    }
  }

  if (!is.null(foldlist)) {
    if (!is.null(nfolds)) {
      stop("Zresidual_CV: supply either `nfolds` or `foldlist`, not both.", call. = FALSE)
    }

    if (!is.list(foldlist) || length(foldlist) < 2L) {
      stop("Zresidual_CV: `foldlist` must be a list containing at least two folds.", call. = FALSE)
    }

    valid_elements <- vapply(
      foldlist,
      function(x) {
        is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
          all(x >= 1) && all(x == floor(x))
      },
      logical(1L)
    )

    if (!all(valid_elements)) {
      stop(
        "Zresidual_CV: every `foldlist` element must be a non-empty finite numeric index vector.",
        call. = FALSE
      )
    }
  }

  if (!is.null(type)) {
    if (!is.character(type) || length(type) != 1L || is.na(type) || !nzchar(type)) {
      stop("Zresidual_CV: `type` must be NULL or one non-empty character string.", call. = FALSE)
    }
  }

  if (!is.null(log_pointpred)) {
    valid_backend <- is.function(log_pointpred) ||
      (is.character(log_pointpred) && length(log_pointpred) == 1L &&
         !is.na(log_pointpred) && nzchar(log_pointpred))

    if (!valid_backend) {
      stop(
        "Zresidual_CV: `log_pointpred` must be NULL, a function, or one function name.",
        call. = FALSE
      )
    }
  }

  if (!is.null(seed)) {
    seed_value <- suppressWarnings(as.numeric(seed))
    if (length(seed_value) != 1L || !is.finite(seed_value) ||
        seed_value != floor(seed_value)) {
      stop("Zresidual_CV: `seed` must be NULL or one integer.", call. = FALSE)
    }
    seed <- as.integer(seed_value)
  }

  if (inherits(object, "coxph")) {
    return(
      Zresidual_CV_survival_coxph(
        object = object,
        data = data,
        nfolds = nfolds,
        foldlist = foldlist,
        nrep = nrep,
        randomized = randomized,
        type = type,
        log_pointpred = log_pointpred,
        seed = seed,
        ...
      )
    )
  }

  if (inherits(object, "survreg")) {
    return(
      Zresidual_CV_survival_survreg(
        object = object,
        data = data,
        nfolds = nfolds,
        foldlist = foldlist,
        nrep = nrep,
        randomized = randomized,
        type = type,
        log_pointpred = log_pointpred,
        seed = seed,
        ...
      )
    )
  }

  stop(
    sprintf(
      "Zresidual_CV: unsupported model class: %s.",
      paste(class(object), collapse = ", ")
    ),
    call. = FALSE
  )
}
