#' Cross-Validated Z-Residuals for Cox Proportional Hazards Models
#'
#' @description
#' Computes cross-validated (CV) Z-residuals for a Cox proportional hazards
#' model (`coxph`), optionally including user-supplied data or performing
#' CV using the model's internal model frame. The function performs
#' K-fold cross-validation, fits the Cox model on each training subset, and
#' computes Z-residuals on the held-out fold.
#'
#' @param fit.coxph A fitted Cox proportional hazards model from
#'   \pkg{survival}, created using `coxph()`.
#'
#' @param data Optional data frame used for cross-validation.
#'   If `NULL`, the model frame of `fit.coxph` is used internally.
#'
#' @param nfolds Integer. Number of cross-validation folds.
#'   If `NULL`, folds are automatically constructed using a parallel
#'   heuristic: number of folds = `10 %/% ncores * ncores`.
#'
#' @param foldlist Optional list specifying fold indices.
#'   If `NULL`, folds are created using `make_fold()`, stratifying by
#'   survival time and censoring status.
#'
#' @param n.rep Integer. Number of repeated Z-residual samples per observation.
#'
#' @details
#' The function:
#'
#' 1. Extracts covariates and event information from either `data` or
#'    the model frame of `fit.coxph`.
#' 2. Creates or uses supplied folds.
#' 3. For each fold:
#'    * Fits the Cox model to the training subset.
#'    * Computes Z-residuals on the held-out subset using
#'      `Zresidual.coxph()`.
#'    * Handles failed model fits by filling the fold with `NA`.
#' 4. Combines results into a matrix of dimension
#'    **n × n.rep**, where *n* is the number of observations.
#'
#' The returned object also stores attributes:
#' * `"Survival.Prob"` – predicted survival probabilities
#' * `"linear.pred"` – Cox linear predictors
#' * `"censored.status"` – censoring indicator
#' * `"covariates"` – original covariate matrix
#' * `"object.model.frame"` – the reconstructed model frame
#'
#' @return
#' A matrix of CV Z-residuals with attributes (listed above). The matrix
#' is not given a specific S3 class by default, but is intended to be used
#' by downstream diagnostic visualization or analysis tools.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' fit <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'
#' # 5-fold CV Z-residuals
#' cvz <- CV.Zresidual.coxph(fit, data = lung, nfolds = 5, n.rep = 10)
#' }
#'
#' @seealso
#' `coxph()`,
#' `Zresidual.coxph()`,
#' `make_fold()`,
#' `CV.Zresidual()`
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel


CV.Zresidual.coxph <- function(fit.coxph, data, nfolds,foldlist,n.rep)
{
  if(!is.null(data)){
    mf <- model.frame(fit.coxph$formula, data)
    nc <- ncol (mf)
    fix_var<-mf[,-1,drop=FALSE]

    if(is.null(nfolds)){
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds<-10 %/% ncores*ncores
    }
    if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,1],
                                              k=nfolds,censor=mf[,1][,2])
    res_fold <- as.list(1:nfolds)
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      data.train <- data[-testcases, ]
      data.test <- data[testcases, ]

      fit_traindata <- tryCatch(
        coxph(fit.coxph$formula,data=data.train),
        error = function(e) NA)

      if(any(is.na(fit_traindata))) {
        res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      }

      if(any(!is.na(fit_traindata))){
        res_fold[[fid]]<-Zresidual.coxph(fit_coxph=fit_traindata,
                                         newdata=data.test,n.rep=n.rep)
      }
    }

    CV.Zresid<-matrix (0, nrow(data) , n.rep)

    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      CV.Zresid[testcases,]<-res_fold[[fid]]
    }
    col_name <- rep(0, n.rep)
    for (i in 1:n.rep) {
      col_name[i] <- paste("CV.Z-residual ", i, sep = "")
    }
    colnames(CV.Zresid) <- col_name
    attr(CV.Zresid,"type")<- "survival"
    attr(CV.Zresid, "Survival.Prob") <- rep(0, nrow(data))
    attr(CV.Zresid, "linear.pred")<- rep(0, nrow(data))
    attr(CV.Zresid,"censored.status")<- rep(0, nrow(data))
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      attr(CV.Zresid, "Survival.Prob")[testcases] <- attr(res_fold[[fid]],"Survival.Prob")
      attr(CV.Zresid, "linear.pred")[testcases] <- attr(res_fold[[fid]],"linear.pred")
      attr(CV.Zresid, "censored.status")[testcases] <- attr(res_fold[[fid]],"censored.status")
    }

    attr(CV.Zresid,"covariates")<- matrix(0,nrow(data),ncol(mf)-1)
    attr(CV.Zresid,"object.model.frame")<- matrix(0,nrow(mf),ncol(mf))
    for(fid in 1:length (foldlist)) {
      attr(CV.Zresid, "covariates")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"covariates"))
      attr(CV.Zresid, "object.model.frame")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"object.model.frame"))
    }
    attr(CV.Zresid,"covariates")<-as.data.frame(attr(CV.Zresid, "covariates"))
    colnames( attr(CV.Zresid,"covariates"))<- colnames(mf[,-1])
    attr(CV.Zresid,"object.model.frame")<-as.data.frame(attr(CV.Zresid, "object.model.frame"))
    colnames( attr(CV.Zresid,"object.model.frame"))<- colnames(mf)

  }

  if(is.null(data)){
    mf <- model.frame.coxph(fit.coxph)
    nc <- ncol (mf)
    fix_var<-mf[,-1,drop=FALSE]

    if(is.null(nfolds)){
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds<-10 %/% ncores*ncores
    }
    if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,1],
                                              k=nfolds,censor=mf[,1][,2])
    res_fold <- as.list(1:nfolds)
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      data.train <- mf[-testcases, ]
      data.test <- mf[testcases, ]

      data.train<-data.frame(t=data.train[[1]][,1], d=data.train[[1]][,2],
                             data.train[,-1])

      names(data.train)[1]<- as.character( (fit.coxph$formula[[2]])[[2]])
      names(data.train)[2]<- as.character( (fit.coxph$formula[[2]])[[3]])

      data.test<-data.frame(t=data.test[[1]][,1], d=data.test[[1]][,2],
                            data.test[,-1])
      rownames(data.test)<-row.names(mf[testcases,])
      names(data.test)[1]<- as.character( (fit.coxph$formula[[2]])[[2]])
      names(data.test)[2]<- as.character( (fit.coxph$formula[[2]])[[3]])

      fit_traindata <- tryCatch(
        coxph(fit.coxph$formula,data=data.train),
        error = function(e) NA)

      if(any(is.na(fit_traindata))) {
        res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      }

      if(any(!is.na(fit_traindata))){
        res_fold[[fid]]<-Zresidual.coxph(fit_coxph=fit_traindata,
                                         newdata=data.test,n.rep=n.rep)
      }
    }

    CV.Zresid<-matrix (0, nrow(mf) , n.rep)

    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      CV.Zresid[testcases,]<-res_fold[[fid]]
    }
    col_name <- rep(0, n.rep)
    for (i in 1:n.rep) {
      col_name[i] <- paste("CV.Z-residual ", i, sep = "")
    }
    colnames(CV.Zresid) <- col_name

    attr(CV.Zresid, "Survival.Prob") <- rep(0, nrow(mf))
    attr(CV.Zresid, "linear.pred")<- rep(0, nrow(mf))
    attr(CV.Zresid,"censored.status")<- rep(0, nrow(mf))
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      attr(CV.Zresid, "Survival.Prob")[testcases] <- attr(res_fold[[fid]],"Survival.Prob")
      attr(CV.Zresid, "linear.pred")[testcases] <- attr(res_fold[[fid]],"linear.pred")
      attr(CV.Zresid, "censored.status")[testcases] <- attr(res_fold[[fid]],"censored.status")
    }

    attr(CV.Zresid,"covariates")<- matrix(0,nrow(mf),ncol(mf)-1)
    attr(CV.Zresid,"object.model.frame")<- matrix(0,nrow(mf),ncol(mf))
    for(fid in 1:length (foldlist)) {
      attr(CV.Zresid, "covariates")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"covariates"))
      attr(CV.Zresid, "object.model.frame")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"object.model.frame"))
    }
    attr(CV.Zresid,"covariates")<-as.data.frame(attr(CV.Zresid, "covariates"))
    colnames( attr(CV.Zresid,"covariates"))<- colnames(mf[,-1])
    attr(CV.Zresid,"object.model.frame")<-as.data.frame(attr(CV.Zresid, "object.model.frame"))
    colnames( attr(CV.Zresid,"object.model.frame"))<- colnames(mf)

  }
  CV.Zresid
}
