#' ma2stage
#'
#' Examines bias for a 2-stage binary outcome design WRT a meta-analysis of non-adaptive designs
#'
#' @param nsims number of simulations (default 1e5)
#' @param n.studies number of non-adaptive studies in meta-analysis
#' @return Data frame of results
#' @export
#'
ma2stage <- function(N=2428, theta0=0.2, theta1=0.15, nsims=1e5, n.studies=4){

  # Obtain bias for single true response probability:
  ad.data <- pwb2stage(N=N, theta0=theta0, theta1=theta1, nsims=nsims)$results

  #### Meta-analysis: Simulate further trials  ####
  #N <- des[names(des)=="n"] # max sample size of Simon design
  C.responses <- vector("list", nsims)
  E.responses <- vector("list", nsims)

  for(i in 1:nsims){
    C.responses[[i]] <- rbinom(n=n.studies, size=N/2, prob=theta0)
    E.responses[[i]] <- rbinom(n=n.studies, size=N/2, prob=theta1)
  }
  C.responses <- do.call(rbind, C.responses)
  E.responses <- do.call(rbind, E.responses)

  pC <- C.responses/(N/2)
  pE <- E.responses/(N/2)

  RDs <- pE - pC
  SEs <- sqrt(pC*(1-pC)/(N/2) + pE*(1-pE)/(N/2))

  #se <- sqrt((theta.hat)*(1-theta.hat)/N)

  # Combine results from Simon simulation and MA simulation:
  RDs <- cbind(ad.data$RD, RDs)
  se <- cbind(ad.data$SE, SEs)


  # Summary estimate and bias for all 5 trials:
  wtdmeans.all <- rep(NA, nsims)
  for(i in 1:nsims){
    wtdmeans.all[i] <- weighted.mean(x=RDs[i, ], w=1/se[i, ]^2)
  }
  bias.all <- wtdmeans.all-(theta1-theta0)


  # Summary estimate and bias, excluding adaptive design if stopped early:
  wtdmeans.exclude <- wtdmeans.all
  for(i in (1:nsims)[ad.data$early.stop==TRUE]){
    wtdmeans.exclude[i] <- weighted.mean(x=RDs[i, 2:(n.studies+1)], w=1/se[i, 2:(n.studies+1)]^2)
  }
  bias.exclude <- wtdmeans.exclude-(theta1-theta0)

  mean.trials <- n.studies+mean(ad.data$early.stop==FALSE)

  # Present results #
  ma.df <- data.frame(reps=rep(nsims, 2),
                      mean.studies=c(n.studies+1, mean.trials),
                      bias=c(mean(bias.all), mean(bias.exclude)),
                      mean.se=c(sd(wtdmeans.all), sd(wtdmeans.exclude)),
                      theta0=rep(theta0, 2),
                      theta1=rep(theta1, 2),
                      RD=rep(theta1-theta0, 2),
                      type=c("All trials", "Exclude early stopped")
                      )
  row.names(ma.df) <- c("All trials", "Exclude early stopped")
  return(ma.df)
}

