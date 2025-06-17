#' A function to calculate pmf of Truncated Poisson
#'
#' @param y y values.
#' @param lambda lambda parameter of TP distribution.
#'

pdf.tp <- function(y, lambda, log.p = FALSE)
{

  n.y <- max(length(y), length(lambda))
  y <- rep(y,length=n.y)
  lambda <- rep(lambda, length=n.y)

  log_prob <- rep (-Inf, n.y)
  i.py <- which (y>0) ## for y > 0

  ## computing log probabilities for positive y only (y>0)
  log_prob [i.py] <- dpois(x=y[i.py],lambda=lambda[i.py], log = TRUE) -
    ppois(q=0,lambda=lambda[i.py], lower.tail = FALSE, log.p = TRUE)

  prob <- exp(log_prob)

  if(log.p) log_prob
  else prob

}
