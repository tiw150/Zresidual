# =========================================================
# Top functions
# =========================================================

#' Posterior predictive checks for count models
#'
#' @description
#' Computes posterior predictive check summaries for a fitted count model using
#' randomized and mid-point probability integral transform residuals, a
#' chi-squared discrepancy, and an optional ANOVA-style residual diagnostic.
#'
#' @param fit A fitted model object.
#' @param data A data frame used for posterior predictive checking.
#' @param predcheck_pointpred Optional low-level backend function. This function
#'   must accept at least `fit`, `data`, and `draw_ids`, and return a list with
#'   components `support`, `family`, `y`, `n`, `ndraws`, `pmf`, `tail`, `rng`,
#'   and `moments`.
#' @param x Optional grouping variable for the ANOVA-style diagnostic. This can
#'   be either a column name in `data` or a vector of length `nrow(data)`.
#' @param ndraws Optional number of posterior draws to use. If `NULL`, all
#'   available draws are used.
#' @param seed Optional random seed used for draw subsampling and randomized
#'   residual generation.
#' @param k_anova Maximum number of bins used when `x` is numeric.
#' @param ... Additional arguments passed to `predcheck_pointpred`.
#'
#' @return
#' An object of class `"predcheck"` containing the predictive-check summary
#' statistics.
#'
#' @examples
#' \dontrun{
#' res <- ppc(fit = fit_nb, data = dat, x = "depth", ndraws = 500, seed = 1)
#' print(res)
#' }
#'
#' @export
ppc <- function(fit,
                data,
                predcheck_pointpred = NULL,
                x = NULL,
                ndraws = NULL,
                seed = NULL,
                k_anova = 10,
                ...) {
  .predictive_check_core(
    fit = fit,
    data = data,
    mode = "ppc",
    predcheck_pointpred = predcheck_pointpred,
    x = x,
    ndraws = ndraws,
    seed = seed,
    k_anova = k_anova,
    ...
  )
}

#' Holdout predictive checks for count models
#'
#' @description
#' Computes holdout predictive check summaries for a fitted count model on a new
#' dataset. The returned summaries use the same discrepancy measures as
#' `ppc()`, but are evaluated on `newdata`.
#'
#' @param fit A fitted model object.
#' @param newdata A data frame used for holdout predictive checking.
#' @param predcheck_pointpred Optional low-level backend function. This function
#'   must accept at least `fit`, `data`, and `draw_ids`, and return a list with
#'   components `support`, `family`, `y`, `n`, `ndraws`, `pmf`, `tail`, `rng`,
#'   and `moments`.
#' @param x Optional grouping variable for the ANOVA-style diagnostic. This can
#'   be either a column name in `newdata` or a vector of length `nrow(newdata)`.
#' @param ndraws Optional number of posterior draws to use. If `NULL`, all
#'   available draws are used.
#' @param seed Optional random seed used for draw subsampling and randomized
#'   residual generation.
#' @param k_anova Maximum number of bins used when `x` is numeric.
#' @param ... Additional arguments passed to `predcheck_pointpred`.
#'
#' @return
#' An object of class `"predcheck"` containing the predictive-check summary
#' statistics.
#'
#' @examples
#' \dontrun{
#' res <- hpc(fit = fit_nb, newdata = dat_hold, x = "depth", ndraws = 500, seed = 1)
#' print(res)
#' }
#'
#' @export
hpc <- function(fit,
                newdata,
                predcheck_pointpred = NULL,
                x = NULL,
                ndraws = NULL,
                seed = NULL,
                k_anova = 10,
                ...) {
  .predictive_check_core(
    fit = fit,
    data = newdata,
    mode = "hpc",
    predcheck_pointpred = predcheck_pointpred,
    x = x,
    ndraws = ndraws,
    seed = seed,
    k_anova = k_anova,
    ...
  )
}

#' Print a predictive-check summary
#'
#' @description
#' Prints the main summary statistics stored in an object of class
#' `"predcheck"`.
#'
#' @param x An object of class `"predcheck"`.
#' @param ... Reserved for future use.
#'
#' @return
#' The input object, invisibly.
#'
#' @export
print.predcheck <- function(x, ...) {
  cat("Predictive check\n")
  cat("  mode:   ", x$mode, "\n", sep = "")
  cat("  family: ", x$family, "\n", sep = "")
  cat("  n:      ", x$n, "\n", sep = "")
  cat("  ndraws: ", x$ndraws, "\n\n", sep = "")
  cat("  chisq:       ", format(x$chisq, digits = 4), "\n", sep = "")
  cat("  sw_rpp:      ", format(x$sw_rpp, digits = 4), "\n", sep = "")
  cat("  aov_rpp:     ", format(x$aov_rpp, digits = 4), "\n", sep = "")
  cat("  sw_mpp:      ", format(x$sw_mpp, digits = 4), "\n", sep = "")
  cat("  mean_sw_pv:  ", format(x$mean_sw_pv, digits = 4), "\n", sep = "")
  cat("  mean_aov_pv: ", format(x$mean_aov_pv, digits = 4), "\n", sep = "")
  invisible(x)
}
