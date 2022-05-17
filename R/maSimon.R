#' maSimon
#'
#' Examines bias for a single Simon design WRT a meta-analysis of non-adaptive designs
#'
#' @param theta true response probablity
#' @param des a realisation of a Simon design, formatted as per output of clinfun::ph2simon
#' @param nsims number of simulations (default 1e5)
#' @param n.studies number of non-adaptive studies in meta-analysis
#' @return Data frame of results
#' @export
#'
maSimon <- function(theta=0.5, des, nsims=1e5, n.studies=4){

  # Obtain bias for single true response probability:
  simon.data <- pwbSimon(theta=theta, des=des, nsims=nsims)$results

  #### Meta-analysis: Simulate further trials  ####
  N <- des[names(des)=="n"] # max sample size of Simon design
  responses <- vector("list", nsims)

  for(i in 1:nsims){
    responses[[i]] <- rbinom(n=n.studies, size=N, prob=theta)
  }
  responses <- do.call(rbind, responses)
  theta.hat <- responses/N
  se <- sqrt((theta.hat)*(1-theta.hat)/N)

  # Combine results from Simon simulation and MA simulation:
  theta.hat <- cbind(simon.data$theta.hat.cor, theta.hat)
  se <- cbind(simon.data$se.cor, se)


  # Summary estimate and bias for all 5 trials:
  wtdmeans.all <- rep(NA, nsims)
  for(i in 1:nrow(theta.hat)){
    wtdmeans.all[i] <- weighted.mean(x=theta.hat[i, ], w=1/se[i, ]^2)
  }
  bias.all <- wtdmeans.all-theta

  bias.all.unwtd <- apply(theta.hat, 1, mean) - theta # Bias using simple mean of estimates

  # Summary estimate and bias, excluding adaptive design if stopped early:
  wtdmeans.exclude <- wtdmeans.all
  for(i in (1:nsims)[simon.data$early.stop==TRUE]){
    wtdmeans.exclude[i] <- weighted.mean(x=theta.hat[i, 2:(n.studies+1)], w=1/se[i, 2:(n.studies+1)]^2)
  }
  bias.exclude <- wtdmeans.exclude-theta

  mean.trials <- n.studies+mean(simon.data$early.stop==FALSE)

  # Present results #
  ma.df <- data.frame(reps=rep(nsims, 3),
                      mean.studies=c(n.studies+1,n.studies+1, mean.trials),
                      bias=c(mean(bias.all.unwtd), mean(bias.all), mean(bias.exclude)),
                      mean.se=c(NA, sd(wtdmeans.all), sd(wtdmeans.exclude)),
                      theta=rep(theta, 3),
                      type=c("All trials (unwtd)", "All trials (wtd)", "Exclude early stopped")
                      )
  row.names(ma.df) <- c("All trials (unwtd)", "All trials (wtd)", "Exclude early stopped")
  return(ma.df)
}

