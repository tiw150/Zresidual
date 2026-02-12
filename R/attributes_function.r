#' Attach standardized metadata attributes to a pred_dist object (optional).
#'
#' If the required inputs are not provided and pred_dist does not already contain
#' the required attributes, this function returns NULL and emits a warning.
#' Otherwise, it returns pred_dist with standardized attributes attached/updated.
#'
#' Required metadata:
#' - family
#' - response
#' - covariates
#' - linear.pred
#'
#' Optional metadata (for zero-related structures):
#' - zero_id
#' - count_id
#'
#' @param pred_dist A list-like object that contains prediction outputs (e.g., lpmf_hat/lcdf_hat).
#' @param family Model/family label (e.g., "HurdleNB", "coxph"). Optional if already present in attributes.
#' @param response Observed response vector. Optional if already present in attributes.
#' @param covariates Data frame of covariates. Optional if already present in attributes.
#' @param linear_pred Numeric vector used for diagnostics. Optional if already present in attributes.
#' @param zero_id Integer indices for zeros/censored observations (optional).
#' @param count_id Integer indices for positive/count observations (optional).
#' @param extra Named list of additional attributes to attach (optional).
#' @param overwrite Logical; if TRUE, overwrite existing attributes with the same name.
#'
#' @return pred_dist with attributes attached, or NULL if required metadata is unavailable.
pred_dist_meta <- function(pred_dist,
                           family = NULL,
                           response = NULL,
                           covariates = NULL,
                           linear_pred = NULL,
                           zero_id = NULL,
                           count_id = NULL,
                           extra = NULL,
                           overwrite = TRUE) {
  if (is.null(pred_dist) || !is.list(pred_dist)) {
    warning("pred_dist_meta: pred_dist is NULL or not a list; returning NULL.")
    return(NULL)
  }
  
  attrs <- attributes(pred_dist)
  if (is.null(attrs)) attrs <- list()
  
  # Helper: use argument value if provided; otherwise fall back to existing attribute.
  pick <- function(arg, name) {
    if (!is.null(arg)) return(arg)
    if (!is.null(attrs[[name]])) return(attrs[[name]])
    NULL
  }
  
  family_use        <- pick(family, "family")
  response_use    <- pick(response, "response")
  covariates_use  <- pick(covariates, "covariates")
  linear_pred_use <- pick(linear_pred, "linear.pred")
  
  missing <- character(0)
  if (is.null(family_use))        missing <- c(missing, "family")
  if (is.null(response_use))    missing <- c(missing, "response")
  if (is.null(covariates_use))  missing <- c(missing, "covariates")
  if (is.null(linear_pred_use)) missing <- c(missing, "linear.pred")
  
  if (length(missing) > 0) {
    warning(sprintf(
      "pred_dist_meta: missing required metadata (%s); returning NULL.",
      paste(missing, collapse = ", ")
    ))
    return(NULL)
  }
  
  meta <- list(
    family = family_use,
    response = response_use,
    covariates = covariates_use,
    linear.pred = linear_pred_use
  )
  
  # Optional ids
  if (!is.null(zero_id))  meta$zero_id  <- zero_id
  if (!is.null(count_id)) meta$count_id <- count_id
  
  # Merge extra attributes (must be named)
  if (!is.null(extra)) {
    if (!is.list(extra) || is.null(names(extra)) || any(names(extra) == "")) {
      stop("pred_dist_meta: extra must be a named list.")
    }
    meta <- c(meta, extra)
  }
  
  # Attach attributes
  if (!overwrite) {
    existing <- names(attributes(pred_dist))
    meta <- meta[!names(meta) %in% existing]
  }
  
  attributes(pred_dist) <- c(attributes(pred_dist), meta)
  pred_dist
}