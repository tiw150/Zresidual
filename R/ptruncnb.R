#' @rdname truncnb
#' @export
ptruncnb <- function(y, mu, size, lower.tail = FALSE, log.p = FALSE) {
  
  n.y <- max(length(y), length(mu), length(size))
  y <- rep(y, length.out = n.y)
  mu <- rep(mu, length.out = n.y)
  size <- rep(size, length.out = n.y)
  
  log_upper_prob <- rep(0, n.y)
  i.py <- which(y >= 1)
  
  log_upper_prob[i.py] <- pnbinom(
    q = y[i.py],
    mu = mu[i.py],
    size = size[i.py],
    lower.tail = FALSE,
    log.p = TRUE
  ) - pnbinom(
    q = 0,
    mu = mu[i.py],
    size = size[i.py],
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  upper_prob <- exp(log_upper_prob)
  log_lower_prob <- log_diff_exp(0, log_upper_prob)
  lower_prob <- exp(log_lower_prob)
  
  if (!lower.tail) {
    if (log.p) log_upper_prob else upper_prob
  } else {
    if (log.p) log_lower_prob else lower_prob
  }
}


