Zresidual.ZI<-function(fit_ZI)
{
  family<-family(fit_ZI)$family
  mu<-predict(fit_ZI, type="conditional")
  size<-sigma(fit_ZI)
  p<-predict(fit_ZI, type="zprob")
  n<-fit_ZI$modelInfo$nobs
  y<-fit_ZI$frame[,1]
  dzpois <- function(x,lambda,p) {
    return((1-p)*dpois(x,lambda)+p*(x==0))
    }
  pzpois <- function(x,lambda,p) {
    return((1-p)*ppois(x,lambda)+p*(x>=0))
    }
  dznbinom <- function(x,size,mu,p) {
    return((1-p)*dnbinom(x,size = size, mu = mu)+p*(x==0))
    }
  pznbinom <- function(x,size,mu,p) {
    return((1-p)*pnbinom(x,size = size, mu = mu)+p*(x>=0))
    }
  if(1*(fit_ZI$modelInfo$allForm$ziformula==~0)==1)##no zero-inflation
  {
    if(family[1]=="gaussian")
    {
      pvalue=pnorm(y,mu,size)
    }
    else if(family[1]=="poisson")
    {
      pvalue=ppois (y-1,mu) + dpois (y,mu) * runif(n)
    }
    else if(family[1]=="nbinom2")
    {
      pvalue=pnbinom (y-1,size=size, mu=mu) + dnbinom (y,size = size, mu = mu) * runif(n)
    }
  }
  else if(1*(fit_ZI$modelInfo$allForm$ziformula==~0)==0)##zero-inflation
  {
    if(family[1]=="poisson")
    {
      pvalue=pzpois(y-1,mu,p) + dzpois (y, mu,p) * runif(n)
    }
    else if(family[1]=="nbinom2")
    {
      pvalue=pznbinom (y-1,size,mu,p) + dznbinom (y,size,mu,p) * runif(n)
    }
  }
  pvalue<-pmin(pmax(pvalue,10^{-10}),1-10^{-10})
  Zresid=qnorm(pvalue)
  list(pvalue=pvalue,Zresid=Zresid)
}




