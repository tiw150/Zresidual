#' A function to calculate ANOVA of Zresidual
#'
#' @param Zresidual A z-residual.
#' @param X Linear predictor or covariate
#' @param k.anova Number of bins if applicable
#' @export
#'

aov.test.zresid <- function (Zresidual, X = c("lp", "covariate"), k.anova = 10) {

  n <- NROW(Zresidual)
  p <- NCOL(Zresidual)
  choices <- c("lp", "covariate")

  is_user_vector <-
    (!missing(X)) &&
    (length(X) == n) &&
    (is.numeric(X) || is.factor(X) || is.logical(X) ||
       (is.character(X) && length(X) == n))

  if (!is_user_vector) {
    if (missing(X)) X <- "lp"
    if (is.character(X) && length(X) > 1 && all(X %in% choices)) {
      X <- match.arg(X, choices)
    }
  }

  if (is_user_vector) {
    xv <- X
    if (is.character(xv)) {
      xv_num <- suppressWarnings(as.numeric(xv))
      xv <- if (all(!is.na(xv_num))) xv_num else droplevels(factor(xv))
    }
    if (is.factor(xv)) xv <- droplevels(xv)
  } else if (is.character(X) && length(X) == 1) {
    if (X == "lp") {
      xv <- attr(Zresidual, "linear.pred")
      if (is.null(xv)) stop("attr(Zresidual, 'linear.pred') not found for X='lp'.")
    } else {
      fitted.value <- attr(Zresidual, "covariates")
      if (is.null(fitted.value)) stop("attr(Zresidual, 'covariates') is missing.")
      cov.name <- variable.names(fitted.value)

      if (X == "covariate") {
        i <- 1L
        cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
            cov.name, "\n")
      } else if (X %in% cov.name) {
        i <- which(cov.name == X)
      } else {
        stop(paste0("X must be one of covariate names: ", paste(cov.name, collapse = ", "), "."))
      }
      xv <- fitted.value[, i]
    }
  } else {
    stop("Invalid X. Provide a length-n vector, or 'lp', or 'covariate', or a covariate name contained in attr(Zresidual,'covariates').")
  }

  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf   <- which(is.infinite(Zresidual) & Zresidual > 0)
  if (length(id.negtv.inf)) Zresidual[id.negtv.inf] <- -1e10
  if (length(id.pos.inf))   Zresidual[id.pos.inf]   <-  1e10

  aov.pv <- rep(NA_real_, p)
  for (j in seq_len(p)) {
    idx <- which(!is.na(Zresidual[, j]) & !is.na(xv))
    if (length(idx) < 3L) { aov.pv[j] <- NA_real_; next }

    z_j  <- Zresidual[idx, j]
    x_jv <- if (is.factor(xv)) droplevels(xv[idx]) else xv[idx]

    aov.pv[j] <- tryCatch(
      test.nl.aov(z_j, x_jv, k.anova = k.anova),
      error = function(e) NA_real_
    )
  }

  return(aov.pv)
}


# aov.test.zresid <- function (Zresidual,X = c("lp", "covariate"), k.anova=10)
# {
#
#   if (missing(X)) X = "lp"
#   if (X == "lp") {
#     fitted.value <- attr(Zresidual, "linear.pred")
#     id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
#     id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
#     Zresidual[id.negtv.inf]<- -1e10
#     Zresidual[id.pos.inf]<- 1e10
#     aov.pv<-rep(0,ncol(Zresidual))
#     for(j in 1:ncol(Zresidual)){
#       aov.pv[j]<- test.nl.aov(Zresidual[,j], fitted.value, k.anova)
#     }
#     aov.pv
#   }
#   if (X != "lp") {
#     fitted.value <- attr(Zresidual, "covariates")
#     if(X == "covariate"){
#       i<-1
#       cat("To plot against other covariates, set X to be the covariate name. Please copy one of the covariate name:",
#           variable.names(fitted.value))
#
#     } else if(X %in% variable.names(fitted.value)){
#       cov.name<-variable.names(fitted.value)
#       i<- which(cov.name==X)
#     } else{stop(paste0("X must be the one of covariate name: ", variable.names(fitted.value),". "))}
#
#
#     id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
#     id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
#     Zresidual[id.negtv.inf]<- -1e10
#     Zresidual[id.pos.inf]<- 1e10
#
#     aov.pv<-rep(0,ncol(Zresidual))
#     for(j in 1:ncol(Zresidual)){
#       id.na <- which(is.na(Zresidual[,j]))
#       count.id <- which(!is.na(Zresidual[,j]))
#       new.Zresidual <- Zresidual[count.id, j]
#       #if(length(id.na) > 0) message("NAs omitted.")
#       aov.pv[j]<- test.nl.aov(new.Zresidual, fitted.value[,i][count.id], k.anova)
#     }
#     aov.pv
#   }
#   return(aov.pv)
# }
