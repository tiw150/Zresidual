#' @importFrom nortest sf.test
#' @export sf.test.zresid
sf.test.zresid <- function (Zresidual)
{
  id.negtv.inf <- which(is.infinite(Zresidual) & Zresidual < 0)
  id.pos.inf <- which(is.infinite(Zresidual) & Zresidual > 0)
  Zresidual[id.negtv.inf]<- -1e10
  Zresidual[id.pos.inf]<- 1e10
  sf.pv<-rep(0,ncol(Zresidual))
  for(i in 1:ncol(Zresidual)){
    sf.pv[i]<-sf.test(Zresidual[,i])$p.value
  }
  sf.pv
}
