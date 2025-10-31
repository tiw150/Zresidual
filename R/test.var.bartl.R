#' Bartlett Test for Homogeneity of Variances of Z-Residual
#'
#' Performs Bartlett's test to assess whether the variance of Z-residuals
#' differs across levels of a covariate. Numeric covariates with many distinct
#' values are binned, and small or empty bins are removed before testing.
#'
#' @param Zresidual Numeric vector of Z-residuals.
#' @param fitted.value Numeric or factor covariate to test against.
#' @param k.bl Integer; the number of bins to discretize a numeric covariate (default 10).
#'
#' @details
#' The function handles covariates as follows:
#' \itemize{
#'   \item If \code{fitted.value} is a factor or has fewer than \code{k.bl} unique values, it is treated as categorical.
#'   \item Otherwise, numeric covariates are binned into \code{k.bl} bins.
#'   \item Bins with fewer than 3 observations are removed.
#'   \item If insufficient bins remain, the covariate is log-transformed and binned again.
#' }
#' Bartlett's test is then applied to the Z-residuals grouped by the factor or binned covariate.
#'
#' @return Numeric p-value from Bartlett's test for homogeneity of variances.
#'
#' @examples
#' \dontrun{
#' Zres <- rnorm(100)
#' x <- runif(100)
#' test.var.bartl(Zres, x, k.bl = 5)
#' }
#'
#' @seealso
#' \code{\link[stats]{bartlett.test}}, \code{\link[stats]{split}}
test.var.bartl <- function(Zresidual, fitted.value, k.bl=10)
{
  Zresidual <- Zresidual[which(!is.na(Zresidual))]
  fitted.value <- fitted.value[which(!is.na(Zresidual))]

  unique.vals <- unique(fitted.value)
  n.unique <- length(unique.vals)

  if (is.factor(fitted.value) || n.unique <= k.bl) {
    # treat as categorical / discrete variable
    if (!is.factor(fitted.value)) {
      fitted.value <- factor(fitted.value)
    }

    lpred.bin <- fitted.value
    Z_group<- split(Zresidual, lpred.bin)
    bl.test<-bartlett.test(Z_group)[["p.value"]]
  }
  if(!is.factor(fitted.value)){
    lpred.bin <- cut(fitted.value, k.bl)
    Z_group<- split(Zresidual, lpred.bin)
    check_Z_group<-rep(k.bl)

    less2_Zgroup<-which(lapply(Z_group,length)<= 2 | is.na(lapply(Z_group,length)))
    is.bins2 <- (k.bl - length(less2_Zgroup))>=2

    if(!is.bins2) {
      fitted.value <- log(fitted.value)
      message("'x' must be a list with at least 2 elements. Fitted values converted to log.")
      lpred.bin <- cut(fitted.value, breaks = quantile(fitted.value, probs = seq(0, 1, length.out = k.bl+1)), include.lowest = TRUE)
      less2_Zgroup<-which(lapply(Z_group,length)<= 2 | is.na(lapply(Z_group,length)))
      Z_group<- split(Zresidual, lpred.bin)
    }

    for(i in 1:k.bl)
    {
      fun<-function(x) x>2
      check_Z_group[i]<-fun(length(Z_group[[i]]))
    }
    if(all(check_Z_group!=0)){
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }else{
      Z_group<-Z_group[-which(check_Z_group==0)]
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }
  }
  bl.test
}
