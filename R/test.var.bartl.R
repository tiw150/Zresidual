#' A function to calculate Bartlett test of Zresidual
#'
#'
#'
#' @param Zresidual A Z-residual.
#' @param fitted.value Fitted values
#' @param k.bl Number of bins if applicable
#'
test.var.bartl <- function(Zresidual, fitted.value, k.bl=10)
{
  Zresidual <- Zresidual[which(!is.na(Zresidual))]
  fitted.value <- fitted.value[which(!is.na(Zresidual))]

  if(is.factor(fitted.value)){
    lpred.bin <- fitted.value
    Z_group<- split(Zresidual, lpred.bin)
    bl.test<-bartlett.test(Z_group)[["p.value"]]
  }
  if(!is.factor(fitted.value)){
    lpred.bin <- cut(fitted.value, k.bl)
    Z_group<- split(Zresidual, lpred.bin)
    check_Z_group<-rep(k.bl)

    less2_Zgroup<-which(lapply(Z_group,length)<= 2 | is.na(lapply(Z_group,length)))
    is.bins2 <- (k.bl - length(less2_Zgroup))>=2

    if(!is.bins2) {
      fitted.value <- log(fitted.value)
      message("'x' must be a list with at least 2 elements. Fitted values converted to log.")
      lpred.bin <- cut(fitted.value, breaks = quantile(fitted.value, probs = seq(0, 1, length.out = k.bl+1)), include.lowest = TRUE)
      less2_Zgroup<-which(lapply(Z_group,length)<= 2 | is.na(lapply(Z_group,length)))
      Z_group<- split(Zresidual, lpred.bin)
    }

    for(i in 1:k.bl)
    {
      fun<-function(x) x>2
      check_Z_group[i]<-fun(length(Z_group[[i]]))
    }
    if(all(check_Z_group!=0)){
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }else{
      Z_group<-Z_group[-which(check_Z_group==0)]
      bl.test<-bartlett.test(Z_group)[["p.value"]]
    }
  }
  bl.test
}
