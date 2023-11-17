test.var.bartl <- function(Zresidual, fitted.values, k.bl=10)
{
  if(is.factor(fitted.values)){
    lpred.bin <- fitted.values
    Z_group<- split(Zresidual, lpred.bin)
    bl.test<-bartlett.test(Z_group)[["p.value"]]
  }
  if(!is.factor(fitted.values)){
    lpred.bin <- cut(fitted.values, k.bl)
    Z_group<- split(Zresidual, lpred.bin)
    check_Z_group<-rep(k.bl)
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

