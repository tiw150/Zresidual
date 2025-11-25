#' Goodness-of-fit test for censored Z-residuals
#'
#' Perform a normality goodness-of-fit test for censored Z-residuals using
#' \code{\link[EnvStats]{gofTestCensored}}. This is typically used to assess
#' whether censored Z-residuals are approximately standard normal under the
#' fitted survival model.
#'
#' Infinite residual values (both positive and negative) are truncated to
#' large finite values (\code{+/- 1e10}) before the test is applied.
#' The censoring indicator is taken from the \code{"censored.status"}
#' attribute of \code{censored.Zresidual}.
#'
#' @importFrom EnvStats gofTestCensored
#'
#' @param censored.Zresidual Numeric vector (or one-column matrix) of censored
#'   Z-residuals. It must carry an attribute \code{"censored.status"} giving
#'   the censoring indicator (1 = event, 0 = censored) for each observation.
#'
#' @return A single numeric value: the p-value from
#'   \code{EnvStats::gofTestCensored} using the Shapiro–Francia test
#'   (\code{test = "sf"}) for a right-censored normal distribution
#'   (\code{distribution = "norm"}). Larger p-values indicate no strong
#'   evidence against normality of the censored Z-residuals.
#'
#' @seealso \code{\link[EnvStats]{gofTestCensored}}
#'
#' @export

gof.censored.zresidual <- function (censored.Zresidual)
{
  id.negtv.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual < 0)
  id.pos.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual > 0)
  censored.Zresidual[id.negtv.inf]<- -1e10
  censored.Zresidual[id.pos.inf]<- 1e10
  censored.status<-attr(censored.Zresidual, "censored.status")
  censored.Zresidual<- as.vector(censored.Zresidual)
  gofTestCensored(censored.Zresidual,censored=censored.status, test = "sf",
                  censoring.side = "right",
                  distribution = "norm")$p.value
}
