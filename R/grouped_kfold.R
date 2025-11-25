#' Create k-fold indices (internal helper)
#'
#' Generate k-fold cross-validation indices with simple stratification
#' based on \code{y}.
#'
#' @param y Outcome used for stratification. Can be a \code{Surv} object,
#'   factor, character, or numeric vector. Numeric vectors are coerced to
#'   factor.
#' @param k Integer; number of folds.
#' @param list Logical; if \code{TRUE}, return a list of index vectors
#'   for each fold, otherwise return an integer vector of fold labels.
#' @param returnTrain Logical; when \code{list = TRUE}, return training
#'   indices instead of test indices.
#'
#' @return A list of index vectors (one per fold) or an integer vector of
#'   fold labels, depending on \code{list}.
#'
#' @keywords internal
#' @noRd


kfold_fn<-function (y, k, list = TRUE, returnTrain = FALSE)
{
  if (class(y)[1] == "Surv")
    y <- y[, "time"]
  if (is.numeric(y)) y<-as.factor(y)
  # {
  #   cuts <- floor(length(y)/k)
  #   if (cuts < 2)
  #     cuts <- 2
  #   if (cuts > 5)
  #     cuts <- 5
  #   breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
  #   y <- cut(y, breaks, include.lowest = TRUE)
  # }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      min_reps <- numInClass[i]%/%k
      if (min_reps > 0) {
        spares <- numInClass[i]%%k
        seqVector <- rep(1:k, min_reps)
        if (spares > 0)
          seqVector <- c(seqVector, sample(1:k, spares))
        foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
      }
      else {
          foldVector[which(y == names(numInClass)[i])] <-
            sample(1:k,size = numInClass[i])
      }
    }
    fold_num<-length(table(foldVector))
    if (fold_num !=k){
      subsets <- replicate(1,sample(length(y)))
      fold <- rep(seq_len(k), length.out = length(y))
      foldVector<-fold[order(subsets)]
    }
  }
  else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  }
  else out <- foldVector
  out
}

#' Build k-fold splits with factor/censoring checks (internal)
#'
#' Generate a k-fold split that preserves factor levels in training sets
#' and avoids empty factor-by-censoring cells in each training fold.
#'
#' @param fix_var Matrix or data frame of covariates.
#' @param y Outcome used to seed the fold assignment (passed to
#'   \code{kfold_fn}).
#' @param k Integer; number of folds.
#' @param censor Vector of censoring/status indicators, same length as \code{y}.
#'
#' @return A list of length \code{k}, each element containing the test indices
#'   for that fold.
#'
#' @keywords internal
#' @noRd

make_fold<-function(fix_var,y,k,censor)
{
  ########Find out which column is factor data
  col<-as.logical(rep(0,ncol(fix_var)))
  for(i in 1:ncol(fix_var)){
    col[i]<-is.factor(fix_var[,i])
  }
  categ_col<-which(col==TRUE)

  #######combine fix variable and censor column#####
  fix_censor<-data.frame(fix_var,censor)
  fix_censor_col<-ncol(fix_censor)

  #####make a fold list
  if(length(categ_col)==0) fold_list<-kfold_fn(y,k)
  if(length(categ_col)!=0){
    counter <- 0
    allfoldisgood<-FALSE
    while(isFALSE(allfoldisgood)){
      fold_list<-kfold_fn(y,k)
      k2<-length(fold_list)
      if(k2 != k) warning(paste("It cannot get", k, "folds,the fold number is",k2))
      test_data <- as.list(1:k2)
      train_data<- as.list(1:k2)
      foldisgood<-matrix(FALSE,nrow =length(categ_col),ncol =k2)

      train_data_censor<-as.list(1:k2)
      check_table<-matrix(FALSE,nrow =length(categ_col),ncol =k2)

      for(i in 1:length(categ_col)){
        for(fid in 1:length (fold_list)) {
          #check factor column in test data within train data
          test_data[[fid]] <- fix_var[fold_list[[fid]], ,drop=FALSE]
          train_data[[fid]]<- fix_var[-fold_list[[fid]], ,drop=FALSE]
          foldisgood[i,fid]<- all(test_data[[fid]][,categ_col[i]]%in%
                                    train_data[[fid]][,categ_col[i]])

          #check cross table contain zero or not in r (factor cross with censor)
          train_data_censor[[fid]] <- fix_censor[-fold_list[[fid]], ,drop=FALSE]
          check_table[i,fid]<-all(apply(table(train_data_censor[[fid]][,categ_col[i]],
                                           train_data_censor[[fid]][,fix_censor_col]),
                                    2,function(x) x!=0))
        }
      }

      allfoldisgood<-all(foldisgood,check_table)
      counter <- counter+1
      if (counter>500) stop("It cannot get the fold list ")
    }
  }
  fold_list
}
