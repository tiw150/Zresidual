Zresidual.hurdle<-function(model_count, model_zero, data)
{
  name.y <- names(model_count$frame)[[1]]
  y <- data[,name.y]
  m<-model_zero$modelInfo$nobs
  family<-family(model_count)$family
  mu<- rep(0, dim(data)[1])
  mu2<-predict(model_count,newdata=data[data[,name.y]>0,],type="conditional")
  mu[which(data[,name.y]>0)] <- mu2
  size<-sigma(model_count)
  prob_nb<-size/(size+mu)
  pi<- 1-predict(model_zero,newdata = data, type="zprob")

  if(family[1]=="truncated_poisson")
  {
    pvalue=pzmpois(y-1,mu,pi)+dzmpois(y,mu,pi)* runif(m)
  }
  else if(family[1]=="truncated_nbinom2")
  {
    pvalue=pzmnbinom(y-1,size,prob_nb,pi)+dzmnbinom(y,size,prob_nb,pi)*runif(m)
  }
  pvalue<-pmin(pmax(pvalue,10^{-10}),1-10^{-10})
  Zresid=qnorm(pvalue)
  list(pvalue=pvalue,Zresid=Zresid)
}


