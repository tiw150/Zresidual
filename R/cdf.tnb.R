#' A function to calculate cdf of Truncated Negative Binomial
#'
#' @param y y values.
#' @param mu mu parameter of TNB distribution.
#' @param size size parameter of TNB distribution.
#'

cdf.tnb <- function(y, mu, size, lower.tail = FALSE, log.p = FALSE) {

  n.y <- max(length(y), length(mu),length(size),length(pi))
  y <- rep(y,length=n.y)
  mu <- rep(mu, length=n.y)
  size <- rep(size, length=n.y)
  # pi <- rep(pi, length=n.y)

  log_upper_prob <- rep (0, n.y)
  #if(count_only) i.py <- which(y>=1) else i.py <- 1:n

  i.py <- which (y>=1) ## for y > 0

  ## computing upper tail probabilities for positive y only (y>0)
  log_upper_prob [i.py] <- pnbinom(q=y[i.py],mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE) -
    pnbinom(q=0,mu=mu[i.py],size=size[i.py], lower.tail = FALSE, log.p = TRUE)

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
