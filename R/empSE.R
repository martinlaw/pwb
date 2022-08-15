#' Empirical standard error
#'
#' Finds empirical standard error for a vector of estimates.
#'
#' @param theta.hat vector estimated response rates
#' @return Estimated empirical standard error
#' @export
#'
empSE <- function(theta.hat){
  nsim <- length(theta.hat)
  theta.bar <- mean(theta.hat)
  est <- sqrt(sum((theta.hat-theta.bar)^2)/(nsim-1))
  est
}

