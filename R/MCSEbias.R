#' Monte Carlo SE of bias
#'
#' Finds Monte Carlo standard error of bias for a vector of estimates.
#'
#' @param theta.hat vector estimated response rates
#'
#' @return Monte Carlo standard error of bias
#'
#' @export
MCSEbias <- function(theta.hat){
  nsim <- length(theta.hat)
  theta.bar <- mean(theta.hat)
  est <- sqrt(sum((theta.hat-theta.bar)^2)/((nsim-1)*nsim))
  est
}

