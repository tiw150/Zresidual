#' Cross-validated Z-residuals for parametric survival regression models
#'
#' @description
#' S3 method for \code{\link{CV.Zresidual}()} applied to parametric survival
#' regression models fitted with \code{\link[survival]{survreg}}. This method
#' performs K-fold cross-validation to obtain external Z-residuals for model
#' diagnostics.
#'
#' @param object A fitted \code{\link[survival]{survreg}} model object.
#' @param nfolds Integer. Number of folds for cross-validation.
#' @param foldlist Optional list specifying custom fold assignments. If
#'   \code{NULL}, folds are generated internally, typically stratified by the
#'   survival response and censoring indicator.
#' @param data Optional data frame used to refit the model during
#'   cross-validation. It is **highly recommended** to supply the original data
#'   here to ensure correct model refitting in each fold, especially when the
#'   original call was complex.
#' @param nrep Integer. Number of repeated Z-residual samples per observation
#'   to generate. Defaults to \code{1}. Each replicate involves re-randomizing
#'   the imputed survival probability for censored observations.
#' @param ... Further arguments passed to the internal worker function
#'   \code{\link{CV_Zresidual_survreg_survival}}.
#' @details
#' This method delegates the actual cross-validation work to
#' \code{\link{CV_Zresidual_survreg_survival}}, which handles the iterative
#' refitting of the \code{survreg} model on $K-1$ folds and computes the
#' randomized Z-residuals on the held-out fold. 
#'
#' The randomized Z-residual, $Z_{ij}$, for the $j$-th observation in the $i$-th fold
#' is computed based on the predicted out-of-sample survival probability $\hat{S}_{\text{train}_i}(t_j)$.
#'
#' The returned object is tagged with class \code{"cvzresid"} in addition
#' to any classes returned by the internal worker.
#'
#' @return
#' An object of class \code{"cvzresid"} containing cross-validated
#' Z-residual diagnostics for the parametric survival model. It is a numeric
#' matrix with $N$ rows and \code{nrep} columns, accompanied by diagnostic attributes
#' (see \code{\link{CV_Zresidual_survreg_survival}} for details).
#'
#' @seealso
#' \code{\link{CV_Zresidual_survreg_survival}}, the generic
#' \code{\link{CV.Zresidual}}, and the `survival` fitting function
#' \code{\link[survival]{survreg}}.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'   # Fit a Weibull model
#'   fit_weibull <- survreg(Surv(time, status) ~ age + sex,
#'                          data = lung, dist = "weibull")
#'   # Compute 5-fold cross-validated Z-residuals
#'   cv_out <- CV.Zresidual(fit_weibull, nfolds = 5, data = lung, nrep = 10)
#'
#'   # Check the first few cross-validated residuals
#'   head(cv_out)
#' }
#'
#' @method CV.Zresidual survreg
#' @export
CV.Zresidual.survreg <- function(object,
                                 nfolds,
                                 foldlist = NULL,
                                 data     = NULL,
                                 nrep     = 1, ...) {

  cv_obj <- CV_Zresidual_survreg_survival(
    fit.survreg = object,
    data        = data,
    nfolds      = nfolds,
    foldlist    = foldlist,
    n.rep       = nrep, ...
  )

  class(cv_obj) <- c("cvzresid", class(cv_obj))
  cv_obj
}

#' @keywords internal
#' @description Internal function to compute cross-validated Z-residuals for
#' **parametric** survival regression models fitted with \code{\link[survival]{survreg}}.
#'
#' @param fit.survreg A fitted \code{\link[survival]{survreg}} model object.
#' @param data Optional \code{data.frame} used for cross-validation. Highly
#'   recommended if the original model was fit without specifying the \code{data}
#'   argument or if \code{foldlist} is supplied.
#' @param nfolds Integer. Number of folds for cross-validation ($K$ in K-fold CV).
#' @param foldlist Optional list specifying custom fold assignments. If \code{NULL},
#'   folds are generated internally, typically stratified by the survival response.
#' @param n.rep Integer. Number of repeated Z-residual samples to generate per
#'   observation (Monte Carlo replications for censored observations).
#' @param ... Additional arguments passed to the residual calculation function
#'   \code{Zresidual_survreg_survival()}.
#'
#' @details
#' This function implements the K-fold cross-validation procedure. In each fold,
#' the parametric survival model (preserving the original distribution, e.g., Weibull)
#' is refitted on the training data. The out-of-sample randomized Z-residuals
#' are then calculated on the held-out test data. 
#'
#' **Data Handling Note:** If \code{data} is \code{NULL}, the internal model frame
#' must be manually reconstructed into a standard data frame (un-packing the
#' \code{Surv} object) before the \code{survreg} model can be successfully refitted
#' on the training subset of each fold. Failed model fits during cross-validation
#' result in \code{NA} residuals for the corresponding test fold.
#'
#' @return A numeric matrix containing the cross-validated Z-residuals ($N \times n.rep$),
#'   where $N$ is the total number of observations. The matrix carries the following
#'   diagnostic attributes:
#' \itemize{
#'   \item \code{Survival.Prob}: Out-of-sample predicted survival probabilities.
#'   \item \code{linear.pred}: Out-of-sample linear predictors (on the \code{survreg} scale).
#'   \item \code{censored.status}: Event indicator (1 = event, 0 = censored).
#'   \item \code{covariates}: Data frame of covariates used.
#'   \item \code{object.model.frame}: The full model frame used for CV residual calculation.
#'   \item \code{type}: Character string, typically \code{"survival"}.
#' }
CV_Zresidual_survreg_survival<- function( fit.survreg, data ,
                                          nfolds ,
                                          foldlist ,
                                          n.rep, ...)
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
