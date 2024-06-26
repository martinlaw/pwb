#' pwbSimon
#'
#' Examines bias for a single Simon design
#'
#' @importFrom clinfun ph2simon
#' @importFrom clinfun ph2single
#' @param theta true response probablity
#' @param des a realisation of a Simon design, formatted as per output of clinfun::ph2simon
#' @param nsims number of simulations (default 1e5)
#' @return list of length 3:
#' results: individual sims
#' ests: estimates
#' mc.error
#' @export

pwbSimon <- function(theta, des, nsims=1e5){
  ##### Setup #####
  n1 <- des[names(des)=="n1"]
  r1 <- des[names(des)=="r1"]
  N <- des[names(des)=="n"] # max sample size
  r <- des[names(des)=="r"]
  reject <- s1 <- s <- numeric(nsims)


  ##### Simulate many trials of this design #####
  for(i in 1:nsims){
    onetrial <- rbinom(n=N, size=1, prob=theta) # Simulate one trial of N results
    s1[i] <- cumsum(onetrial)[n1] # Number of responses at n1
    s[i] <- sum(onetrial) # Number of responses at N
  }

  ##### Combine results into data frame #####
  results <- data.frame(early.stop=s1<=r1)
  results$responses <- s
  results$responses[results$early.stop] <- s1[results$early.stop]
  results$reject <- s>r & results$early.stop==FALSE
  results$n <- rep(N, nsims)
  results$n[results$early.stop] <- n1
  results$complete <- results$n==N

  results$theta.hat <- results$responses/results$n
  results$se <- sqrt(results$theta.hat * (1-results$theta.hat) / results$n)

  # Add correction for (rare) instances of zero responses and all responses:
  results$responses.cor <- results$responses
  results$responses.cor[results$responses.cor==0] <- 0.5 # Add 0.5 when zero responses.
  # Subtract 0.5 when all responses:
  results$responses.cor[results$responses.cor==results$n] <- results$responses.cor[results$responses.cor==results$n]-0.5
  # Calculate SE using corrected responses. Note: generally identical to uncorrected.
  results$theta.hat.cor <- results$responses.cor/results$n
  results$se.cor <- sqrt(results$theta.hat.cor * (1-results$theta.hat.cor) / results$n)

  early.stop <- results[results$complete==FALSE, ] # Trials stopped early
  stop.at.N <- results[results$complete==TRUE, ] # Trials not stopped early


  ##### Analyse results #####
  early.stop.theta.bar <- mean(early.stop$theta.hat.cor)
  stop.at.N.theta.bar <- mean(stop.at.N$theta.hat.cor)
  all.theta.bar <- mean(results$theta.hat.cor)

  # Empirical SE:
  early.stop.emp.SE <- sqrt((1/(nrow(early.stop)-1)) * sum((early.stop$theta.hat.cor-early.stop.theta.bar)^2))
  stop.at.N.emp.SE <- sqrt((1/(nrow(stop.at.N)-1)) * sum((stop.at.N$theta.hat.cor-stop.at.N.theta.bar)^2))
  all.emp.SE <- sqrt((1/(nsims-1)) * sum((results$theta.hat.cor-all.theta.bar)^2))
  emp.SE <- c(early.stop.emp.SE,
              stop.at.N.emp.SE,
              all.emp.SE,
              NA)

  # SE:
  mean.SE <- c(mean(early.stop$se.cor),
               mean(stop.at.N$se.cor),
               mean(results$se.cor),
               weighted.mean(results$se.cor, w=1/(results$se.cor)^2))

  # Estimate of response rate:
  theta.bar.vec <- c(early.stop.theta.bar,
                     stop.at.N.theta.bar,
                     all.theta.bar,
                     weighted.mean(x=results$theta.hat.cor, w=1/(results$se.cor)^2))
  bias <- theta.bar.vec-theta

  trials <- c(nrow(early.stop),
              nrow(stop.at.N),
              nsims,
              nsims)

  ests <- data.frame(nsims=trials,
                     bias,
                     mean.SE,
                     emp.SE,
                     theta)
  row.names(ests) <- c("Stopped early", "Complete", "All", "All (precision-weighted)")
  ests$type <- c("Stopped early (unadjusted)", "Complete (unadjusted)", "All (unadjusted)", "All (precision-weighted)")

  mc.error.bias.all <- sqrt(sum((results$theta.hat-all.theta.bar)^2) / (nsims*(nsims-1)) )
  mc.error.se.all <- all.emp.SE/sqrt(2*(nsims-1))

  all.data <- list(results=results,
                   ests=ests,
                   mc.error=c(mc.error.bias.all, mc.error.se.all))
  return(all.data)
}
