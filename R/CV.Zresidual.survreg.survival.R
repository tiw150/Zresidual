#' Cross-validated Z-residuals for parametric survival regression models
#'
#' @description
#' S3 method for \code{CV.Zresidual()} applied to parametric survival
#' regression models fitted with \code{survival::survreg()}.
#'
#' @param object A fitted \code{survival::survreg} model object.
#' @param nfolds Integer. Number of folds for cross-validation.
#' @param foldlist Optional list specifying custom fold assignments. If
#'   \code{NULL}, folds are generated internally.
#' @param data Optional data frame used to refit the model during
#'   cross-validation. Required when \code{foldlist} is supplied or when
#'   the original model call does not contain the data explicitly.
#' @param nrep Integer. Number of repeated cross-validations to perform.
#'   Default is 1.
#' @param ... Further arguments passed to the internal worker function.
#'
#' @details
#' This method delegates the actual cross-validation work to
#' \code{CV_Zresidual_survreg_survival()}, which performs fold
#' construction, model refitting, and computation of cross-validated
#' Z-residuals.
#'
#' The returned object is tagged with class \code{"cvzresid"} in addition
#' to any classes returned by the internal worker.
#'
#' @return
#' An object of class \code{"cvzresid"} containing cross-validated
#' Z-residual diagnostics for the parametric survival model.
#'
#' @seealso
#' \code{CV_Zresidual_survreg_survival()} and the generic
#' \code{CV.Zresidual()}.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'   fit_weibull <- survreg(Surv(time, status) ~ age + sex,
#'                          data = lung, dist = "weibull")
#'   cv_out <- CV.Zresidual(fit_weibull, nfolds = 5, data = lung)
#' }
#'
#' @method CV.Zresidual survreg
#' @export
CV.Zresidual.survreg <- function(object,
                                 nfolds,
                                 foldlist = NULL,
                                 data     = NULL,
                                 nrep     = 1,
                                 ...) {

  cv_obj <- CV_Zresidual_survreg_survival(
    fit_survreg = object,
    data        = data,
    nfolds      = nfolds,
    foldlist    = foldlist,
    n.rep       = nrep,
    ...
  )

  class(cv_obj) <- c("cvzresid", class(cv_obj))
  cv_obj
}


#' Cross-validated Z-residuals for parametric survival regression models
#' (internal worker)
#'
#' @description
#' `CV_Zresidual_survreg_survival()` is the internal workhorse used by
#' [CV.Zresidual.survreg()] to compute cross-validated (CV) Z-residuals
#' for parametric survival regression models fitted by
#' \code{\link[survival]{survreg}}.
#'
#' It performs K-fold cross-validation, refits the \code{survreg} model on
#' each training subset, and computes Z-residuals on the held-out fold via
#' \code{Zresidual_survreg_survival()}.
#'
#' @param fit_survreg A fitted parametric survival regression model from
#'   \pkg{survival}, created using \code{\link[survival]{survreg}}.
#'
#' @param data Optional \code{data.frame} used for cross-validation. It must
#'   contain the survival response and all covariates appearing in
#'   \code{fit_survreg$terms}. If \code{NULL}, the model frame of
#'   \code{fit_survreg} (via \code{model.frame.survreg()}) is used internally.
#'
#' @param nfolds Integer. Number of cross-validation folds. If \code{NULL},
#'   the number of folds is chosen heuristically by the caller (e.g. the
#'   S3 method \code{CV.Zresidual.survreg()}).
#'
#' @param foldlist Optional list specifying fold indices. Each element should
#'   be an integer vector giving the row indices of the held-out (test) set
#'   for a given fold. If \code{NULL}, folds are created using
#'   \code{make_fold()}, stratifying by the survival response and censoring
#'   indicator.
#'
#' @param n.rep Integer. Number of repeated Z-residual samples per observation
#'   in each fold (i.e. number of Monte Carlo replications for censored
#'   observations in \code{Zresidual_survreg_survival()}).
#'
#' @param ... Further arguments passed to lower-level helper functions
#'   (if any).
#'
#' @details
#' The function works in two modes:
#' \itemize{
#'   \item If \code{data} is not \code{NULL}, folds are defined on the rows
#'         of \code{data}. For each fold, the model is refitted on
#'         \code{data[-test, ]} and Z-residuals are computed on
#'         \code{data[test, ]}.
#'   \item If \code{data} is \code{NULL}, the internal model frame of
#'         \code{fit_survreg} is used. The function reconstructs explicit
#'         time and status columns from the \code{Surv} response before
#'         refitting the \code{survreg} model within each fold.
#' }
#'
#' For each fold, `CV_Zresidual_survreg_survival()` attempts to refit the model
#' using the same formula and distribution as in \code{fit_survreg}. If the
#' model fit fails (due to convergence or other errors/warnings), the
#' corresponding fold residuals are filled with \code{NA}.
#'
#' The per-fold Z-residual matrices are then stacked into a single
#' \eqn{n \times} \code{n.rep} matrix, where \eqn{n} is the total number of
#' observations. Attributes such as survival probabilities, linear predictors,
#' censoring status, covariates, and the model frame are reassembled in the
#' original observation order.
#'
#' @return
#' A numeric matrix of dimension \eqn{n \times} \code{n.rep}, where \eqn{n}
#' is the number of rows in \code{data} (if supplied) or in the internal
#' model frame of \code{fit_survreg}. Columns are typically named
#' \code{"CV.Z-residual 1"}, \code{"CV.Z-residual 2"}, ..., up to
#' \code{n.rep}. The matrix usually carries attributes such as:
#' \itemize{
#'   \item \code{"type"}: character string \code{"survival"}.
#'   \item \code{"Survival.Prob"}: vector of predicted survival probabilities
#'         at the observed time for each observation.
#'   \item \code{"linear.pred"}: vector of linear predictors on the
#'         \code{survreg} scale.
#'   \item \code{"censored.status"}: event indicator (1 = event, 0 = censored).
#'   \item \code{"covariates"}: data frame of covariates used for the
#'         residual computation.
#'   \item \code{"object.model.frame"}: data frame representing the model
#'         frame underlying the CV residuals.
#' }
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   data(lung)
#'   fit_wb <- survreg(Surv(time, status) ~ age + sex,
#'                     data = lung, dist = "weibull")
#'
#'   cvz_wb <- CV_Zresidual_survreg_survival(
#'     fit_survreg = fit_wb,
#'     data        = lung,
#'     nfolds      = 5,
#'     foldlist    = NULL,
#'     n.rep       = 10
#'   )
#' }
#'
#' @seealso
#' \code{\link[survival]{survreg}},
#' \code{Zresidual_survreg_survival()},
#' \code{CV.Zresidual.survreg()},
#' \code{make_fold()}
#'
#' @keywords internal
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
CV_Zresidual_survreg_survival<- function( fit.survreg,  data, nfolds,foldlist,n.rep)
{
  if(!is.null(data)){
    mf <- model.frame(fit.survreg$terms, data)
    nc <- ncol (mf)
    fix_var<-mf[,-1,drop=FALSE]

    if(is.null(nfolds))
    {
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds<-10 %/% ncores*ncores
    }
    if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,1],
                                              k=nfolds,censor=mf[,1][,2])
    res_fold <- as.list(1:nfolds)
    for(fid in 1:length (foldlist))
    {
      testcases <- foldlist[[fid]]
      data.train <- data[-testcases, ]
      data.test <- data[testcases, ]

      fit_traindata <- tryCatch(
        survreg(formula=fit.survreg$terms,data=data.train,dist=fit.survreg$dist),
        error = function(e) NA,
        warning = function(w) NA)

      if(any(is.na(fit_traindata))) {
        res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      }

      if(any(!is.na(fit_traindata))){
        res_fold[[fid]]<-Zresidual.survreg(fit_survreg=fit_traindata,newdata = data.test,n.rep=n.rep)
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
    mf <- model.frame.survreg(fit_survreg)
    nc <- ncol (mf)
    fix_var<-mf[,-1,drop=FALSE]

    if(is.null(nfolds))
    {
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds<-10 %/% ncores*ncores
    }
    if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,1],
                                              k=nfolds,censor=mf[,1][,2])
    res_fold <- as.list(1:nfolds)
    for(fid in 1:length (foldlist))
    {
      testcases <- foldlist[[fid]]
      data.train <- mf[-testcases, ]
      data.test <- mf[testcases, ]

      data.train<-data.frame(t=data.train[[1]][,1], d=data.train[[1]][,2],
                             data.train[,-1])

      names(data.train)[1]<- as.character( (fit.survreg$terms[[2]])[[2]])
      names(data.train)[2]<- as.character( (fit.survreg$terms[[2]])[[3]])

      data.test<-data.frame(t=data.test[[1]][,1], d=data.test[[1]][,2],
                            data.test[,-1])
      rownames(data.test)<-row.names(mf[testcases,])
      names(data.test)[1]<- as.character( (fit.survreg$terms[[2]])[[2]])
      names(data.test)[2]<- as.character( (fit.survreg$terms[[2]])[[3]])

      fit_traindata <- tryCatch(
        survreg(formula=fit.survreg$terms,data=data.train,dist=fit.survreg$dist),
        error = function(e) NA,
        warning = function(w) NA)

      if(any(is.na(fit_traindata))) {
        res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      }

      if(any(!is.na(fit_traindata))){
        res_fold[[fid]]<-Zresidual_survreg_survival(fit_survreg=fit_traindata,newdata = data.test,n.rep=n.rep)
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
