#' ANOVA Test for Z-Residuals
#'
#' Performs an ANOVA test to assess whether Z-residuals differ across
#' levels of a covariate. Numeric covariates with many distinct values
#' are discretized into bins. Small or empty bins are removed before testing.
#'
#' @param Zresidual Numeric vector of Z-residuals.
#' @param fitted.value Numeric or factor covariate to test against.
#' @param k.anova Integer; the maximum number of bins to discretize a numeric covariate (default 10).
#'
#' @details
#' The function handles covariates as follows:
#' \itemize{
#'   \item If \code{fitted.value} is a factor or has fewer than \code{k.anova} unique values, it is treated as categorical.
#'   \item Otherwise, numeric covariates are binned into \code{k.anova} bins using \code{\link{cut}}.
#'   \item Bins with fewer than 3 observations are removed.
#'   \item If insufficient bins remain, the covariate is log-transformed and binned again.
#' }
#' ANOVA is then performed with \code{lm(Zresidual ~ binned_covariate)} and the p-value for the first term is returned.
#'
#' @return Numeric p-value from the ANOVA F-test for the effect of the covariate on Z-residuals.
#'
#' @examples
#' \dontrun{
#' Zres <- rnorm(100)
#' x <- runif(100)
#' test.nl.aov(Zres, x, k.anova = 5)
#' }
#'
#' @seealso
#' \code{\link[stats]{anova}}, \code{\link[stats]{lm}}
#'
test.nl.aov <- function(Zresidual, fitted.value, k.anova=10)
{
  unique.vals <- unique(fitted.value)
  n.unique <- length(unique.vals)

  if (is.factor(fitted.value) || n.unique <= k.anova) {
    # treat as categorical / discrete variable
    if (!is.factor(fitted.value)) {
      fitted.value <- factor(fitted.value)
    }
    lpred.bin <- as.numeric(fitted.value) - 1
    anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
  } else if(!is.factor(fitted.value)){
    lpred.bin <- droplevels(cut(fitted.value, k.anova))
    less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
    is.bins2 <- (nlevels(lpred.bin) - length(less2_factor))>2

    if(!is.bins2) {
    fitted.value <- log(fitted.value)
    message("Contrasts can be applied only to factors with 2 or more levels. Fitted values converted to log.")
    lpred.bin <- droplevels(cut(fitted.value, k.anova))
    less2_factor<-which(tapply(lpred.bin,lpred.bin,length)<= 2)
    }

    if(rlang::is_empty(names(less2_factor))){
      anova(lm(Zresidual ~ lpred.bin))$`Pr(>F)`[1]
    }else {
      list_less2_factor<-list()
      for(j in 1:length(less2_factor)){
        list_less2_factor[[j]]<-which(lpred.bin==names(less2_factor[j]))
      }
      vector_less2_factor<-unlist(list_less2_factor, use.names = FALSE)
      new.lpred.bin<- lpred.bin[-vector_less2_factor]
      new.Zresidual<-Zresidual[-vector_less2_factor]
      anova(lm(new.Zresidual ~ new.lpred.bin))$`Pr(>F)`[1]
    }
  } else {
    warning("Problem in ANOVA.")
  }
}
