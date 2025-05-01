#' A function to calculate cdf of Truncated Poisson
#'
#' @param y y values.
#' @param lambda lambda parameter of TP distribution.
#'

cdf.tp<- function(y, lambda, lower.tail = FALSE, log.p = FALSE) {

  n.y <- max(length(y), length(lambda))
  y <- rep(y,length=n.y)
  lambda <- rep(lambda, length=n.y)

  log_upper_prob <- rep (0, n.y)
  i.py <- which (y>=1) ## for y > 0

  ## computing upper tail probabilities for positive y only (y>0)
  log_upper_prob [i.py] <- ppois(q=y[i.py],lambda = lambda[i.py], lower.tail = FALSE, log.p = TRUE) -
    ppois(q=0,lambda = lambda[i.py], lower.tail = FALSE, log.p = TRUE)

  upper_prob <- exp(log_upper_prob)
  log_lower_prob <- log_diff_exp(0, log_upper_prob)
  lower_prob <- exp (log_lower_prob)

  if (lower.tail==FALSE)
  {
    if(log.p) log_upper_prob
    else upper_prob
  } else
  {
    if(log.p) log_lower_prob
    else lower_prob
  }
}
