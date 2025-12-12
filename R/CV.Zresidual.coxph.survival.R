#' Cross-validated Z-residuals for Cox proportional hazards models
#'
#' @description
#' S3 method for the generic function \code{CV.Zresidual()} applied to
#' Cox proportional hazards models fitted with \code{survival::coxph()}.
#' This method computes the cross-validated Z-residuals, automatically detecting
#' whether a frailty term is present in the model formula. It then dispatches
#' to the appropriate internal implementation for standard Cox models or shared
#' frailty Cox models. 
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @param object A fitted \code{survival::coxph} model object.
#' @param nfolds Integer. The number of folds for cross-validation (K in K-fold CV).
#' @param foldlist Optional list specifying custom fold assignments. If
#'   \code{NULL}, folds are generated internally, often ensuring a balanced
#'   distribution of events.
#' @param data Optional data frame used to refit the model during
#'   cross-validation. It is required if the original model call did not contain
#'   the data explicitly, or when \code{foldlist} is supplied.
#' @param nrep Integer. Number of repeated cross-validations to perform.
#'   Default is 1. Repeating CV provides a more stable estimate of the residuals.
#' @param ... Further arguments passed to the internal worker functions.
#'
#' @details
#' The method determines the correct worker function by examining
#' \code{attr(object$terms, "specials")$frailty}:
#' \itemize{
#'   \item **Standard Cox Models**: If no frailty term is found, it calls
#'     \code{\link{CV_Zresidual_coxph_survival}}.
#'   \item **Shared Frailty Cox Models**: If a frailty term (e.g., \code{frailty(group)})
#'     is present, it calls \code{\link{CV_Zresidual_coxph_frailty_survival}}.
#' }
#'
#' The returned object is a matrix of Z-residuals with crucial diagnostic
#' information stored as attributes.
#'
#' @return
#' A matrix of class \code{"cvzresid"} (and others inherited from the internal
#' worker) with dimension \eqn{N \times nrep}{N x nrep}
#' (where \eqn{N}{N} is the number of observations).
#' The matrix columns contain the cross-validated Z-residuals.
#'
#' The object also includes the following diagnostic attributes:
#' \itemize{
#'   \item \code{Survival.Prob}: Cross-validated predicted survival probabilities
#'     \eqn{\hat S_i(t_i)}{S_hat_i(t_i)}.
#'   \item \code{linear.pred}: Cross-validated linear predictors
#'     \eqn{\eta_i = \mathbf{x}_i^\top \hat{\boldsymbol{\beta}}}{eta_i = x_i^T beta_hat}.
#'   \item \code{censored.status}: The event indicator from the survival object.
#'   \item \code{covariates}: The covariates used in the model.
#'   \item \code{object.model.frame}: The model frame of the original data.
#' }
#'
#' @seealso
#' \code{\link[survival]{coxph}}, \code{\link{CV_Zresidual_coxph_survival}},
#' \code{\link{CV_Zresidual_coxph_frailty_survival}},
#' and the generic \code{\link{CV.Zresidual}}.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'   # Example 1: Standard Cox Model
#'   fit_std <- coxph(Surv(time, status) ~ age + sex, data = lung)
#'   cv_out_std <- CV.Zresidual(fit_std, nfolds = 5, data = lung)
#'
#'   # Example 2: Shared Frailty Cox Model (assuming the lung data has a group column 'inst')
#'   # lung$inst_factor <- as.factor(lung$inst)
#'   # fit_frty <- coxph(Surv(time, status) ~ age + sex + frailty(inst_factor), data = lung)
#'   # cv_out_frty <- CV.Zresidual(fit_frty, nfolds = 5, data = lung)
#' }
#'
#' @method CV.Zresidual coxph
#' @export
CV.Zresidual.coxph <- function(object,
                               nfolds,
                               foldlist = NULL,
                               data     = NULL,
                               nrep     = 1, ...) {

  frailty_terms <- attr(object$terms, "specials")$frailty

  if (!is.null(frailty_terms)) {
    ## shared frailty Cox model
    cv_obj <- CV_Zresidual_coxph_frailty_survival(
      fit_coxph = object,
      data      = data,
      nfolds    = nfolds,
      foldlist  = foldlist,
      n.rep     = nrep, ...
    )
  } else {
    ## standard Cox model (no frailty)
    cv_obj <- CV_Zresidual_coxph_survival(
      fit_coxph = object,
      data      = data,
      nfolds    = nfolds,
      foldlist  = foldlist,
      n.rep     = nrep, ...
    )
  }

  class(cv_obj) <- c("cvzresid", class(cv_obj))
  cv_obj
}

#' @title Cross-validated Z-residuals for Standard Cox Models
#' @name CV_Zresidual_coxph_survival
#' @keywords internal
#' @description Internal function to compute cross-validated Z-residuals for
#' **standard** Cox proportional hazards models (without frailty).
#'
#' @param fit.coxph A fitted \code{survival::coxph} model object.
#' @param data Optional data frame containing the data. Required if the original
#'   model was fit without specifying the \code{data} argument or if \code{foldlist} is supplied.
#' @param nfolds Integer. Number of folds for cross-validation (K in K-fold CV).
#' @param foldlist Optional list specifying custom fold assignments.
#' @param n.rep Integer. Number of repeated cross-validations.
#' @param ... Additional arguments.
#'
#' @details
#' This function implements the K-fold cross-validation procedure. In each fold,
#' the Cox model is refitted on the training data, and the Z-residuals are
#' calculated on the held-out test data. This approach is for models without
#' shared frailty terms.
#'
#' @return A matrix containing the cross-validated Z-residuals (\eqn{N \times nrep}{N x nrep}) with
#'  diagnostic attributes: \code{Survival.Prob}, \code{linear.pred},
#'  \code{censored.status}, \code{covariates}, and \code{object.model.frame}.
CV_Zresidual_coxph_survival <- function(fit.coxph, data, nfolds,foldlist,n.rep, ...)
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
        res_fold[[fid]]<-Zresidual_coxph_survival(fit_coxph=fit_traindata,
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
        res_fold[[fid]]<-Zresidual_coxph_survival(fit_coxph=fit_traindata,
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

#' @title Cross-validated Z-residuals for Shared Frailty Cox Models
#' @name CV_Zresidual_coxph_frailty_survival
#' @keywords internal
#' @description Internal function to compute cross-validated Z-residuals for
#' **shared frailty** Cox proportional hazards models.
#'
#' @param fit.coxph A fitted \code{survival::coxph} model object containing a frailty term.
#' @param data Optional data frame containing the data. Required if the original
#'   model was fit without specifying the \code{data} argument or if \code{foldlist} is supplied.
#' @param nfolds Integer. Number of folds for cross-validation (K in K-fold CV).
#' @param foldlist Optional list specifying custom fold assignments.
#' @param n.rep Integer. Number of repeated cross-validations.
#' @param ... Additional arguments.
#'
#' @details
#' This function implements the K-fold cross-validation procedure specifically
#' for shared frailty Cox models (i.e., those containing a \code{frailty()} term).
#' In each fold, the model is refitted on the training data, and the Z-residuals
#' are calculated on the held-out test data.
#'
#' @return A matrix containing the cross-validated Z-residuals (\eqn{N \times nrep}{N x nrep}) with
#'  diagnostic attributes: \code{Survival.Prob}, \code{linear.pred},
#'  \code{censored.status}, \code{covariates}, and \code{object.model.frame}.
CV_Zresidual_coxph_frailty_survival<- function( fit.coxph, data, nfolds,foldlist,n.rep, ...)
{
  if (is.null(data)) {
    mf <- model.frame.coxph(fit.coxph)
    mf_nc <- ncol (mf)
    groupid <- as.factor(mf[, ncol(mf)])
    if (!is.factor(groupid))
      stop("The group ID must be factor!")
    gpnumber <- length(levels(groupid))
    mm <- model.matrix.coxph(fit.coxph)
    mm_nc <- ncol (mm)
    fix_var <- mf[, -c(1, mf_nc), drop = FALSE]

    if (is.null(nfolds))
    {
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds <- 10 %/% ncores * ncores
    }
    if (is.null(foldlist))
      foldlist <- make_fold(
        fix_var = fix_var,
        y = mf[, mf_nc],
        k = nfolds,
        censor = mf[, 1][, 2]
      )
    res_fold <- as.list(1:nfolds)
    for (fid in 1:length (foldlist))
    {
      testcases <- foldlist[[fid]]
      data.train <- mf[-testcases,]
      data.test <- mf[testcases,]

      data.train<-data.frame(t=data.train[[1]][,1], d=data.train[[1]][,2],
                             data.train[,-c(1,ncol(data.train))],
                             id=data.train[,ncol(data.train)])

      names(data.train)[1]<- as.character( (fit.coxph$formula[[2]])[[2]])
      names(data.train)[2]<- as.character( (fit.coxph$formula[[2]])[[3]])
      names(data.train)[ncol(data.train)]<- group_id_name<-gsub(".*[(]([^.]+)[,].*", "\\1",
                                                                (fit.coxph$formula)[[3]])[3]
      data.train[,ncol(data.train)]<-as.factor(data.train[,ncol(data.train)])
      #frailty.form<-((fit.coxph$formula)[[3]])[[3]]
      #frailty.distr<-gsub(".*[(]([^.]+)[,].*", "\\1", frailty.form)[3]

      data.test<-data.frame(t=data.test[[1]][,1], d=data.test[[1]][,2],
                            data.test[,-c(1,ncol(data.test))],
                            id=data.test[,ncol(data.test)])
      rownames(data.test)<-row.names(mf[testcases,])
      names(data.test)[1]<- as.character( (fit.coxph$formula[[2]])[[2]])
      names(data.test)[2]<- as.character( (fit.coxph$formula[[2]])[[3]])
      names(data.test)[ncol(data.test)]<- group_id_name<-gsub(".*[(]([^.]+)[,].*", "\\1",
                                                              (fit.coxph$formula)[[3]])[3]
      data.test[,ncol(data.test)]<-as.factor(data.test[,ncol(data.test)])

      fit_traindata <- tryCatch(
        coxph(fit.coxph$formula,data=data.train),
        error = function(e)
          NA
      )

      if (any(is.na(fit_traindata))) {
        res_fold[[fid]] <- rep(NA, length(foldlist[[fid]]))
      }

      if (any(!is.na(fit_traindata))) {
        res_fold[[fid]] <- Zresidual.coxph.frailty(
          fit_coxph = fit_traindata,
          traindata = data.train,
          newdata = data.test,
          n.rep = n.rep
        )
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

    attr(CV.Zresid,"type")<- "survival"
    attr(CV.Zresid, "Survival.Prob") <- rep(0, nrow(mf))
    attr(CV.Zresid, "linear.pred")<- rep(0, nrow(mf))
    attr(CV.Zresid,"censored.status")<- rep(0, nrow(mf))
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      attr(CV.Zresid, "Survival.Prob")[testcases] <- attr(res_fold[[fid]],"Survival.Prob")
      attr(CV.Zresid, "linear.pred")[testcases] <- attr(res_fold[[fid]],"linear.pred")
      attr(CV.Zresid, "censored.status")[testcases] <- attr(res_fold[[fid]],"censored.status")
    }

    attr(CV.Zresid,"covariates")<- matrix(0,nrow(mf),ncol(mm)-1)
    attr(CV.Zresid,"object.model.frame")<- matrix(0,nrow(mf),ncol(mf))
    for(fid in 1:length (foldlist)) {
      attr(CV.Zresid, "covariates")[foldlist[[fid]],] <- attr(res_fold[[fid]],"covariates")
      attr(CV.Zresid, "object.model.frame")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"object.model.frame"))
    }
    attr(CV.Zresid,"object.model.frame")<-as.data.frame(attr(CV.Zresid, "object.model.frame"))
    colnames( attr(CV.Zresid,"object.model.frame"))<- colnames(mf)
  }


  if (!is.null(data)){
    mf <- model.frame(fit.coxph$formula, data)
    nc <- ncol (mf)
    fix_var<-mf[,-c(1,nc),drop=FALSE]

    if(is.null(nfolds))
    {
      ncores <- detectCores()
      cl <- makeForkCluster(ncores)
      registerDoParallel(cl)
      nfolds<-10 %/% ncores*ncores
    }
    if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,nc],
                                              k=nfolds,censor=mf[,1][,2])
    res_fold <- as.list(1:nfolds)
    for(fid in 1:length (foldlist))
    {
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
        res_fold[[fid]]<-Zresidual_coxph_frailty_survival(fit_coxph = fit_traindata,
                                                 traindata = data.train,
                                                 newdata = data.test,
                                                 n.rep = n.rep
        )
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

    attr(CV.Zresid, "Survival.Prob") <- rep(0, nrow(data))
    attr(CV.Zresid, "linear.pred")<- rep(0, nrow(data))
    attr(CV.Zresid,"censored.status")<- rep(0, nrow(data))
    for(fid in 1:length (foldlist)) {
      testcases <- foldlist[[fid]]
      attr(CV.Zresid, "Survival.Prob")[testcases] <- attr(res_fold[[fid]],"Survival.Prob")
      attr(CV.Zresid, "linear.pred")[testcases] <- attr(res_fold[[fid]],"linear.pred")
      attr(CV.Zresid, "censored.status")[testcases] <- attr(res_fold[[fid]],"censored.status")
    }

    attr(CV.Zresid,"covariates")<- matrix(0,nrow(mf),ncol(mf))
    attr(CV.Zresid,"object.model.frame")<- matrix(0,nrow(mf),ncol(mf))
    for(fid in 1:length (foldlist)) {
      attr(CV.Zresid, "covariates")[foldlist[[fid]],] <- attr(res_fold[[fid]],"covariates")
      attr(CV.Zresid, "object.model.frame")[foldlist[[fid]],] <- as.matrix(attr(res_fold[[fid]],"object.model.frame"))
    }
    attr(CV.Zresid,"object.model.frame")<-as.data.frame(attr(CV.Zresid, "object.model.frame"))
    colnames( attr(CV.Zresid,"object.model.frame"))<- colnames(mf)

  }
  CV.Zresid
}


