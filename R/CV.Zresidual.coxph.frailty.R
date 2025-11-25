#' Cross-Validated Z-Residuals for Shared Frailty Cox Models
#'
#' @description
#' Computes cross-validated (CV) Z-residuals for a shared frailty Cox
#' proportional hazards model fitted with \code{\link[survival]{coxph}}
#' using a term such as \code{frailty(group)}. The function performs
#' K-fold cross-validation, refits the frailty Cox model on each training
#' subset, and computes (possibly replicated) Z-residuals on the held-out
#' fold via \code{Zresidual.coxph.frailty()}.
#'
#' @param fit.coxph A fitted shared frailty Cox proportional hazards model
#'   from \pkg{survival}, created using \code{\link[survival]{coxph}} with
#'   a frailty term (e.g. \code{Surv(time, status) ~ x1 + x2 + frailty(group)}).
#'
#' @param data Optional \code{data.frame} used for cross-validation. It must
#'   contain the survival response, all fixed-effect covariates, and the
#'   frailty grouping factor used in \code{fit.coxph$formula}. If \code{NULL},
#'   the model frame of \code{fit.coxph} (via \code{model.frame.coxph()}) is
#'   used internally.
#'
#' @param nfolds Integer. Number of cross-validation folds. If \code{NULL},
#'   the number of folds is chosen heuristically based on the number of
#'   available cores, approximately as \code{10 \%/\%  ncores * ncores}.
#'
#' @param foldlist Optional list specifying fold indices. Each element should
#'   be an integer vector giving the row indices of the held-out (test) set
#'   for a given fold. If \code{NULL}, folds are created using
#'   \code{make_fold()}, stratifying primarily by the frailty grouping factor
#'   and ensuring that factor-by-censoring combinations are represented in
#'   each training set.
#'
#' @param n.rep Integer. Number of repeated Z-residual samples per observation
#'   in each fold (i.e. number of Monte Carlo replications for censored
#'   observations in \code{Zresidual.coxph.frailty}).
#'
#' @details
#' The function operates in two modes:
#' \itemize{
#'   \item If \code{data} is not \code{NULL}, the folds are constructed (or
#'     interpreted) in terms of the rows of \code{data}, and each fold fit
#'     uses \code{data[-test, ]} as the training set and \code{data[test, ]}
#'     as the test set.
#'   \item If \code{data} is \code{NULL}, the internal model frame of
#'     \code{fit.coxph} is used. In this case the function reconstructs a
#'     data frame with explicit time and status columns plus covariates and
#'     the frailty grouping factor, consistent with \code{fit.coxph$formula},
#'     before refitting the model in each fold.
#' }
#'
#' In both cases, the workflow is:
#' \enumerate{
#'   \item Extract the model frame and identify the frailty grouping factor
#'         (the grouping variable in \code{frailty()}).
#'   \item Create or use supplied folds via \code{foldlist}. By default,
#'         \code{make_fold()} is called with the covariates, the grouping
#'         factor, and the event indicator to enforce reasonable balance.
#'   \item For each fold:
#'     \itemize{
#'       \item Fit the shared frailty Cox model to the training subset.
#'       \item Compute Z-residuals on the held-out subset using
#'             \code{Zresidual.coxph.frailty()}.
#'       \item If the model fit fails, fill the corresponding fold residuals
#'             with \code{NA}.
#'     }
#'   \item Stack the per-fold residual matrices into a single
#'         \eqn{n \times} \code{n.rep} matrix, where \eqn{n} is the total
#'         number of observations.
#'   \item Collect and re-assemble attributes (survival probabilities, linear
#'         predictors, censoring status, covariates, and model frame) in the
#'         original observation order.
#' }
#'
#' @return
#' A numeric matrix of dimension \eqn{n \times} \code{n.rep}, where \eqn{n}
#' is the number of rows in \code{data} (if supplied) or in the internal
#' model frame of \code{fit.coxph}. Columns are named
#' \code{"CV.Z-residual 1"}, \code{"CV.Z-residual 2"}, \dots, up to
#' \code{n.rep}. The matrix has the following attributes attached:
#' \itemize{
#'   \item \code{"type"}: character string \code{"survival"}.
#'   \item \code{"Survival.Prob"}: vector of predicted survival probabilities
#'         at the observed time for each observation, assembled from all folds.
#'   \item \code{"linear.pred"}: vector of fixed-effect linear predictors.
#'   \item \code{"censored.status"}: event indicator (1 = event, 0 = censored).
#'   \item \code{"covariates"}: matrix or data frame of covariates used in
#'         the Z-residual computation (structure depends on whether
#'         \code{data} is supplied).
#'   \item \code{"object.model.frame"}: data frame representing the model
#'         frame underlying the CV residuals (assembled from all folds).
#' }
#'
#' This object is intended for downstream diagnostic and visualization
#' tools (e.g. QQ-plots, residual-vs-covariate plots) tailored to shared
#' frailty Cox models.
#'
#' @examples
#' \dontrun{
#'   library(survival)
#'
#'   ## Example data with cluster/grouping variable
#'   data(lung)
#'   set.seed(1)
#'   lung$cluster_id <- factor(sample(1:10, size = nrow(lung), replace = TRUE))
#'
#'   ## Shared frailty Cox model
#'   fit_frail <- coxph(Surv(time, status) ~ age + sex + frailty(cluster_id),
#'                      data = lung)
#'
#'   ## 5-fold CV Z-residuals with 10 Monte Carlo replications
#'   cvz_frail <- CV.Zresidual.coxph.frailty(
#'     fit.coxph = fit_frail,
#'     data      = lung,
#'     nfolds    = 5,
#'     foldlist  = NULL,
#'     n.rep     = 10
#'   )
#' }
#'
#' @seealso
#' \code{\link[survival]{coxph}},
#' \code{Zresidual.coxph.frailty},
#' \code{CV.Zresidual.coxph},
#' \code{make_fold}
#'
#' @export
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel


CV.Zresidual.coxph.frailty<- function( fit.coxph, data, nfolds,foldlist,n.rep)
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
        res_fold[[fid]]<-Zresidual.coxph.frailty(fit_coxph = fit_traindata,
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
