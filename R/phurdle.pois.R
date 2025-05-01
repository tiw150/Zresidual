#' A function to calculate predictive probability of hurdle poisson.
#'
#' @param y y values.
#' @param lambda lambda parameter.
#' @param pi hu parameter.
#'

phurdle.pois <- function(y, lambda, pi, lower.tail = FALSE, log.p = FALSE)
{
  log_diff_exp <- function (la, lb)
  {
    ifelse(la >= lb, la + log1mexp(la-lb), NaN)
  }

  # compute the logarithm of the sum of exponential values.
  log_sum_exp <- function (lx)
  {
    mlx <- max (lx)
    log (sum (exp (lx - mlx))) + mlx
  }

  n.y <- max(length(y), length(lambda),length(pi))
  y <- rep(y,length=n.y)
  lambda <- rep(lambda, length=n.y)
  pi <- rep(pi, length=n.y)

  # this code computes log upper tail probability
  log_upper_prob <- rep (0, n.y)
  i.py <- which (y>0) ## for y < 0, upper tail probability == 1 (log=0)
  i.zy <- which (y==0)

  ## computing upper tail probabilities for zero y only
  log_upper_prob [i.zy] <- log_diff_exp(0, log(pi[i.zy]))

  ## computing upper tail probabilities for positive y only
  log_upper_prob [i.py] <-
    log_diff_exp(0, log(pi[i.py])) +
    ppois(q=y[i.py], lambda=lambda[i.py], lower.tail = FALSE, log.p = TRUE) +
    - ppois(q=0, lambda=lambda[i.py], lower.tail = FALSE, log.p = TRUE)

  #- log_diff_exp(0, dnbinom(x=0,mu=mu[i.py],log=TRUE))
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
