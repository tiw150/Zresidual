#' A function to calculate pdf of Truncated Negative Binomial
#'
#' @param y y values.
#' @param mu mu parameter of TNB distribution.
#' @param size size parameter of TNB distribution.
#' @export

pdf.tnb <- function(y, mu, size, log.p = FALSE)
{

  n.y <- max(length(y), length(mu),length(size),length(pi))
  y <- rep(y,length=n.y)
  mu <- rep(mu, length=n.y)
  size <- rep(size, length=n.y)

  log_prob <- rep (-Inf, n.y)
  #if(count_only) i.py <- which(y>0) else i.py <- 1:n
  i.py <- which (y>0) ## for y > 0

  ## computing log probabilities for positive y only (y>0)
  log_prob [i.py] <- dnbinom(x=y[i.py],mu=mu[i.py],size=size[i.py], log = TRUE) -
    pnbinom(q=0,mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE)

  prob <- exp(log_prob)

  if(log.p) log_prob
  else prob

}
