#' Meta-analysis simulation (binary outcome)
#'
#' Simulates binary data for a meta-analysis
#'
#' @param N sample size (shared for all trials)
#' @param theta response rate (shared for all trials)
#' @param n.studies number of trials to simulate
#' @param nsims number of meta-analyses to simulate
#'
#' @export
#'
simMeta <- function(N, theta, n.studies, nsims){

  responses <- vector("list", nsims)

  for(i in 1:nsims){
    responses[[i]] <- rbinom(n=n.studies, size=N, prob=theta)
  }
  responses <- do.call(rbind, responses)

  # Correction in case there are 0 responses or N responses:
  responses[responses==0] <- 0.5
  responses[responses==N] <- N-0.5

  theta.hat <- responses/N
  se <- sqrt((theta.hat)*(1-theta.hat)/N)


  colnames(theta.hat) <- paste("Study ", 1:n.studies, sep="")
  colnames(se) <- paste("Study ", 1:n.studies, sep="")

  MA.output <- list(theta.hat=theta.hat,
                    se=se)
  MA.output
}


