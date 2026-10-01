ptweedie_saddle_tail <- function(q, mu, phi, power, lower.tail = TRUE) {
  lambda <- mu^(2 - power) / (phi * (2 - power))
  alpha  <- (2 - power) / (power - 1)
  gam    <- phi * (power - 1) * mu^(power - 1)
  
  that  <- (1 - (lambda * alpha * gam / q)^(1/(alpha + 1))) / gam
  Kt    <- lambda * ((1 - gam * that)^(-alpha) - 1)
  Kppt  <- lambda * alpha * (alpha + 1) * gam^2 * (1 - gam * that)^(-alpha - 2)
  
  w <- sign(that) * sqrt(2 * (that * q - Kt))
  u <- that * sqrt(Kppt)
  
  leading   <- stats::pnorm(w)
  corr_term <- stats::dnorm(w) * (1/w - 1/u)
  Fy <- leading + corr_term
  
  ratio <- abs(corr_term / leading)   # diagnostic: should be << 1 for a trustworthy result
  
  if (!lower.tail) Fy <- 1 - Fy
  
  list(cdf = Fy, ratio = ratio)
}