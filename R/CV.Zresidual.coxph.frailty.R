#' CV Zresidual for shared fraitly model using coxph function
#'
#' @importFrom parallel detectCores
#' @importFrom foreach registerDoParallel
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
