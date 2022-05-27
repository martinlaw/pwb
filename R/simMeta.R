#' Meta-analysis simulation (binary outcome)
#'
#' Simulates binary data for a meta-analysis.
#'
#' @param N sample size (shared for all trials if scalar)
#' @param theta response rate (shared for all trials)
#' @param n.studies number of trials to simulate
#' @param nsims number of meta-analyses to simulate
#' @return A list of two vectors of length N: response rate theta hat and SE
#' @export
#'
simMeta <- function(N, theta, n.studies, nsims){

  if(length(N)!=1 & length(N)!=n.studies) stop("N must be length 1 (shared sample size) or length equal to n.studies")

  if(length(N)==1) N <- rep(N, n.studies) # Convert N into a vector if a scalar

  responses <- vector("list", nsims)
  for(i in 1:nsims){
      responses[[i]] <- rbinom(n=n.studies, size=N, prob=theta)
  }

  responses <- do.call(rbind, responses)

  # Correction in case there are 0 responses or N responses:
  responses[responses==0] <- 0.5
  for(i in 1:n.studies){
    responses[responses==N[i]] <- N[i]-0.5
  }

  theta.hat <- t(apply(responses, 1, function(x) x/N))

  se <- matrix(NA, nrow=nsims, ncol=n.studies)
  for(i in 1:n.studies){
     se[, i] <- sqrt(theta.hat[, i]*(1-theta.hat[, i])/N[i])
  }

  colnames(theta.hat) <- paste("Study ", 1:n.studies, sep="")
  colnames(se) <- paste("Study ", 1:n.studies, sep="")

  MA.output <- list(theta.hat=theta.hat,
                    se=se)
  MA.output
}


