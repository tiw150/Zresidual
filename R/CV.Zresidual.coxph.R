cv.zresidual.coxph <- function(fit, data, nfolds=NULL,foldlist=NULL)
{
  # Required packages:
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "foreach",
    "parallel",
    "doParallel",
    "stringr"
  )

  mf <- model.frame(fit$formula, data)
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
      coxph(fit$formula,data=data.train),
      error = function(e) NA,
      warning = function(w) NA)

    if(any(is.na(fit_traindata))) {
      res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
    }

    if(any(!is.na(fit_traindata))){
      res_fold[[fid]]<-Zresidual.coxph(fit_coxph=fit_traindata,newdata = data.test)
    }
  }
  cv_Zresid<-cv_censored.Zresid<-cv_SP<-cv_ucs<-cv_mcs<-cv_martg<-cv_dev<-rep(0, nrow(data))

  for(fid in 1:length (foldlist)) {
    testcases <- foldlist[[fid]]
    cv_Zresid[testcases]<-res_fold[[fid]]$Zresid
    cv_censored.Zresid[testcases]<-res_fold[[fid]]$censored.Zresid
    cv_SP[testcases]<-res_fold[[fid]]$SP
    cv_ucs[testcases]<-res_fold[[fid]]$ucs
    cv_mcs[testcases]<-res_fold[[fid]]$mcs
    cv_martg[testcases]<-res_fold[[fid]]$martg
    cv_dev[testcases]<-res_fold[[fid]]$dev
  }

  list(cv_Zresid=cv_Zresid,cv_censored.Zresid=cv_censored.Zresid,cv_SP=cv_SP,
       cv_ucs=cv_ucs,cv_mcs=cv_mcs,cv_martg=cv_martg,cv_dev=cv_dev)
}
