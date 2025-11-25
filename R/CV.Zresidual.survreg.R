#' Cross-Validated Z-Residuals for Parametric Survival Regression Models
#'
#' @description
#' Computes cross-validated (CV) Z-residuals for a parametric survival
#' regression model fitted by \code{\link[survival]{survreg}}. The function
#' performs K-fold cross-validation, refits the \code{survreg} model on each
#' training subset, and computes Z-residuals on the held-out fold via
#' \code{Zresidual.survreg()}.
#'
#' @param fit.survreg A fitted parametric survival regression model from
#'   \pkg{survival}, created using \code{\link[survival]{survreg}}.
#'
#' @param data Optional \code{data.frame} used for cross-validation. It must
#'   contain the survival response and all covariates appearing in
#'   \code{fit.survreg$terms}. If \code{NULL}, the model frame of
#'   \code{fit.survreg} (via \code{model.frame.survreg()}) is used internally.
#'
#' @param nfolds Integer. Number of cross-validation folds. If \code{NULL},
#'   the number of folds is chosen heuristically based on the number of
#'   available cores, approximately as \code{10 %/% ncores * ncores}.
#'
#' @param foldlist Optional list specifying fold indices. Each element should
#'   be an integer vector giving the row indices of the held-out (test) set
#'   for a given fold. If \code{NULL}, folds are created using
#'   \code{make_fold()}, stratifying by the survival response and censoring
#'   indicator.
#'
#' @param n.rep Integer. Number of repeated Z-residual samples per observation
#'   in each fold (i.e. number of Monte Carlo replications for censored
#'   observations in \code{Zresidual.survreg}).
#'
#' @details
#' The function works in two modes:
#' \itemize{
#'   \item If \code{data} is not \code{NULL}, folds are defined on the rows
#'         of \code{data}. For each fold, the model is refitted on
#'         \code{data[-test, ]} and Z-residuals are computed on
#'         \code{data[test, ]}.
#'   \item If \code{data} is \code{NULL}, the internal model frame of
#'         \code{fit.survreg} is used. The function reconstructs explicit
#'         time and status columns from the \code{Surv} response before
#'         refitting the \code{survreg} model within each fold.
#' }
#'
#' For each fold, \code{CV.Zresidual.survreg()} attempts to refit the model
#' using the same formula and distribution as in \code{fit.survreg}. If the
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
#' model frame of \code{fit.survreg}. Columns are named
#' \code{"CV.Z-residual 1"}, \code{"CV.Z-residual 2"}, \dots, up to
#' \code{n.rep}. The matrix has the following attributes attached:
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
#'   ## 5-fold CV Z-residuals with 10 Monte Carlo replications
#'   cvz_wb <- CV.Zresidual.survreg(
#'     fit.survreg = fit_wb,
#'     data        = lung,
#'     nfolds      = 5,
#'     foldlist    = NULL,
#'     n.rep       = 10
#'   )
#' }
#'
#' @seealso
#' \code{\link[survival]{survreg}},
#' \code{Zresidual.survreg},
#' \code{CV.Zresidual.coxph},
#' \code{CV.Zresidual.coxph.frailty},
#' \code{make_fold}
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel


CV.Zresidual.survreg<- function( fit.survreg,  data, nfolds,foldlist,n.rep)
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
        res_fold[[fid]]<-Zresidual.survreg(fit_survreg=fit_traindata,newdata = data.test,n.rep=n.rep)
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
