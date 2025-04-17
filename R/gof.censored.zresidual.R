#' @export gof.censored.zresidual
gof.censored.zresidual <- function (censored.Zresidual)
{
  id.negtv.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual < 0)
  id.pos.inf <- which(is.infinite(censored.Zresidual) & censored.Zresidual > 0)
  censored.Zresidual[id.negtv.inf]<- -1e10
  censored.Zresidual[id.pos.inf]<- 1e10
  censored.status<-attr(censored.Zresidual, "censored.status")
  censored.Zresidual<- as.vector(censored.Zresidual)
  gofTestCensored(censored.Zresidual,censored=censored.status, test = "sf",
                  censoring.side = "right",
                  distribution = "norm")$p.value
}
