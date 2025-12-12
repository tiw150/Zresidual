#' @title Print Analysis of Entropy Results
#' @concept tabulation
#' @description Standard print method for objects of class \code{anoen_results}.
#'   Displays the analysis of entropy table with options for formatted output (Markdown or gt).
#'
#' @param x An object of class \code{anoen_results} returned by \code{\link{anoen}}.
#' @param digits Integer; the number of decimal places to use. Default is 4.
#' @param test A character string specifying which hypothesis tests to display.
#'   \itemize{
#'     \item \code{"all"} (default): Displays both Chi-squared and Sequential F-test p-values.
#'     \item \code{"chi"}: Displays only the Asymptotic Chi-squared p-values.
#'     \item \code{"F"}: Displays only the recommended Sequential F-test p-values.
#'   }
#' @param format A character string specifying the output format.
#'   \itemize{
#'     \item \code{"text"} (default): Standard console output.
#'     \item \code{"md"}: Returns a Markdown table using \code{knitr::kable}.
#'     \item \code{"gt"}: Returns a \code{gt} table object with parsed math headers and spanning groups.
#'   }
#' @param ... Additional arguments passed to the formatting functions.
#'
#' @return Invisibly returns \code{x} (for text), or the formatted table object (for md/gt).
#' @export
print.anoen_results <- function(x, digits = 4, test = c("all", "chi", "F"),
                                format = c("text", "md", "gt"), ...) {

    # 1. Resolve arguments
    test <- match.arg(test)
    format <- match.arg(format)

    # 2. Check dependencies and FALLBACK if needed
    if (format == "md" && !requireNamespace("knitr", quietly = TRUE)) {
        warning("Package 'knitr' is not installed/available. Falling back to text output.")
        format <- "text"
    }

    if (format == "gt" && !requireNamespace("gt", quietly = TRUE)) {
        warning("Package 'gt' is not installed/available. Falling back to text output.")
        format <- "text"
    }

    # 3. Select columns based on 'test'
    df_print <- x$results
    cols_all <- colnames(df_print)

    # Define which indices to keep based on the test selection
    idx_keep <- seq_along(cols_all)

    if (test == "chi") {
        idx_keep <- which(cols_all != "f_pvalue")
    } else if (test == "F") {
        idx_keep <- which(cols_all != "chi_pvalue")
    }

    df_print <- df_print[, idx_keep]

    # Retrieve specific headers for the selected columns
    # Try attribute first, fall back to list component (Robust fix)
    headers_all <- attr(x$results, "latex_headers")
    if (is.null(headers_all) && !is.null(x$latex_headers)) {
        headers_all <- x$latex_headers
    }

    # Ensure we have headers matching columns length
    if (length(headers_all) >= max(idx_keep)) {
        headers_selected <- headers_all[idx_keep]
    } else {
        # Fallback if headers are completely missing/mismatched
        headers_selected <- colnames(df_print)
    }

    # --- FORMAT: GT (Grammar of Tables) ---
    if (format == "gt") {
        # Construct the named list for cols_label
        col_map <- stats::setNames(
            lapply(headers_selected, gt::md),
            colnames(df_print)
        )

        # Identify numeric columns for formatting
        num_cols <- names(df_print)[sapply(df_print, is.numeric)]

        gt_obj <- gt::gt(df_print) |>
            # 1. Apply LaTeX Headers
            gt::cols_label(.list = col_map) |>

            # 2. Add Spanning Headers (The Visual Groups)
            gt::tab_spanner(
                label = "Model Complexity",
                columns = dplyr::any_of(c("p", "Delta_par", "df_c"))
            ) |>
            gt::tab_spanner(
                label = "Entropy Info",
                columns = dplyr::any_of(c("d_adj", "I_H"))
            ) |>
            gt::tab_spanner(
                label = gt::md("Entropy-based $R^2$"),
                columns = dplyr::any_of(c("partialR2H", "compR2H", "R2H"))
            ) |>
            gt::tab_spanner(
                label = "Significance",
                columns = dplyr::any_of(c("chi_pvalue", "f_pvalue"))
            ) |>

            # 3. Formatting & Style
            gt::fmt_number(columns = dplyr::all_of(num_cols), decimals = digits) |>
            gt::sub_missing(missing_text = "") |>
            gt::tab_header(
                title = gt::md("**Analysis of Entropy (ANOEN)**"),
                subtitle = paste("Model Family:", x$model_type)
            )

        return(gt_obj)
    }

    # --- FORMAT: MARKDOWN (knitr::kable) ---
    if (format == "md") {
        # Replace NAs with placeholder before printing
        is_num <- sapply(df_print, is.numeric)
        df_print[is_num] <- lapply(df_print[is_num], round, digits)
        df_print[is.na(df_print)] <- ""

        return(knitr::kable(df_print, col.names = headers_selected, digits = digits, ...))
    }

    # --- FORMAT: TEXT (Default Console Output) ---
    cat("\nAnalysis of Entropy (ANOEN) Table\n")
    cat("Model Family:", x$model_type, "\n")
    if (!is.null(x$normalization)) {
        cat("Normalization:", x$normalization, "\n")
    }

    # Round numeric columns
    is_num <- sapply(df_print, is.numeric)
    df_print[is_num] <- lapply(df_print[is_num], round, digits)
    df_print[is.na(df_print)] <- "."

    print(df_print, row.names = FALSE)
    return(invisible(x))
}
