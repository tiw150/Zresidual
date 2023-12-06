Zresidual.hurdle<-function(fit_hd, n.rep=nrep)
{
  y <- fit_hd$frame[,1]
  covs<- as.data.frame(fit_hd$frame[,-1])
  colnames(covs)<-colnames(fit_hd$frame)[-1]
  family<-family(fit_hd)$family
  m<-fit_hd$modelInfo$nobs
  fixed.terms<-attr(fit_hd$modelInfo$terms$cond$fixed,"term.labels")
  fixed.pars<- fit_hd$fit$par[c(1:(length(fixed.terms)+1))]
  zi.terms<- attr(fit_hd$modelInfo$reTrms$zi$terms$fixed,"term.labels")
  zi.pars<- fit_hd$fit$par[-c(1:(length(fixed.terms)+1))]

  dhurdle.poisson <- function(x,lambda,pi, k=0,log=FALSE) {
    y <- ifelse(x<=k, log(pi),
                log(1 - pi)+dpois(x,lambda,log=TRUE)-
                  log(1- dpois(k,lambda)))
    if (!log) {y <- exp(y) }
  }

  phurdle.poisson <- function(q,lambda,pi, k=0,lower.tail = TRUE,
                              log.p = FALSE) {
    y <- log(1 - pi) + ppois(q,lambda, lower.tail = FALSE, log.p = TRUE)
    y <- y - log(1 - dpois(x=k,lambda))
    if (lower.tail) {
      y <- 1 - exp(y)
      if (log.p) { y <- log(y) }
    } else {
      if (!log.p) { y <- exp(y)}
    }
    y
  }

  dhurdle.nb<- function(x,mu, size, pi, k=0,log=FALSE) {
    y <- ifelse(x<=k, log(pi),
                log(1 - pi) + dnbinom(x=x,mu=mu,size=size, log=TRUE)-
                  log(1- dnbinom(x=k,mu=mu,size=size)))
    if (!log) {y <- exp(y) }
  }

  phurdle.nb <- function(q, mu, size, pi, k=0,lower.tail = TRUE,
                         log.p = FALSE) {
    y <- log(1 - pi) + pnbinom(q=q,mu=mu,size=size,lower.tail = FALSE, log.p = TRUE)
    y <- y - log(1- dnbinom(x=k,mu=mu,size=size))
    if (lower.tail) {
      y <- 1 - exp(y)
      if (log.p) { y <- log(y) }
    } else {
      if (!log.p) { y <- exp(y)}
    }
    y
  }

  Zresid <- matrix(0, n, n.rep)
  col_name <- rep(0, n.rep)
  for (i in 1:n.rep) {
  if(family[1]=="truncated_poisson"){

    lp.fixed<- as.numeric(fixed.pars[1])+ (as.matrix(covs[fixed.terms]) %*% as.numeric(fixed.pars[-1]))
    mu <- as.vector(exp(lp.fixed))
    lp.zi<- as.numeric(zi.pars[1])+ (as.matrix(covs[zi.terms]) %*% as.numeric(zi.pars[-1]))
    p.zi<- exp(lp.zi)/(1+exp(lp.zi))
    pvalue=phurdle.poisson(y,mu,p.zi,lower.tail=FALSE) + dhurdle.poisson(y,mu,p.zi)* runif(m)

  }else if(family[1]=="truncated_nbinom2"){
    lp.fixed<- as.numeric(fixed.pars[1])+ (as.matrix(covs[fixed.terms]) %*% as.numeric(fixed.pars[-1]))
    mu <- as.vector(exp(lp.fixed))
    lp.zi<- as.numeric(zi.pars[1])+ (as.matrix(covs[zi.terms]) %*% as.numeric(zi.pars[-c(1,length(zi.pars))]))
    p.zi<- exp(lp.zi)/(1+exp(lp.zi))
    size<-sigma(fit_hd)
    # prob_nb<-size/(size+mu)
    pvalue=phurdle.nb(q=y,mu=mu,size=size,pi=p.zi,lower.tail=FALSE) + dhurdle.nb(x=y,mu=mu,size=size,pi=p.zi)*runif(m)
  }
    pvalue <- pmin(pmax(pvalue, 10 ^ {-10 }), 1 - 10 ^ { -10})
    Zresid[,i] = -qnorm(pvalue)
    col_name[i]<-paste("Z-residual ",i,sep = "")
  }
  colnames(Zresid)<- col_name

  #######
  zero.ind<- ifelse(y == 0, 1, 0)

  Zresid.value<-as.matrix(Zresid)

  attributes(Zresid.value) <- c(attributes(Zresid.value), list(
    zero.indicator= zero.ind,
    linear.pred = log(mu),
    covariates = covs
  ))
  return(Zresid.value)
}




