#### Find Simon designs and select optimal designs ####
all.designs <- clinfun::ph2simon(pu=0.5, pa=0.6, ep1=0.05, ep2=0.1, nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
des <- all.designs$out[opt.index, ]


maSimon <- function(theta=0.5, nsims=10000, n.studies=4, des){

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
  theta.hat <- cbind(simon.data$theta.hat, theta.hat)
  se <- cbind(simon.data$se, se)


  # Summary estimate and bias for all 5 trials:
  wtdmeans.all <- rep(NA, nsims)
  for(i in 1:nrow(theta.hat)){
    wtdmeans.all[i] <- weighted.mean(x=theta.hat[i, ], w=1/se[i, ]^2)
  }
  bias.all <- wtdmeans.all-theta


  # Summary estimate and bias, excluding adaptive design if stopped early:
  wtdmeans.exclude <- wtdmeans.all
  for(i in (1:nsims)[simon.data$early.stop==TRUE]){
    wtdmeans.exclude[i] <- weighted.mean(x=theta.hat[i, 2:(n.studies+1)], w=1/se[i, 2:(n.studies+1)]^2)
  }
  bias.exclude <- wtdmeans.exclude-theta

  mean.trials <- n.studies+mean(simon.data$early.stop==FALSE)

  # Present results #
  ma.df <- data.frame(reps=rep(nsims, 2),
                      mean.studies=c(n.studies+1, mean.trials),
                      bias=round(c(mean(bias.all), mean(bias.exclude)), 5),
                      mean.se=round(c(sd(wtdmeans.all), sd(wtdmeans.exclude)), 3)
                      )
  row.names(ma.df) <- c("All trials", "Exclude early stopped")
  return(ma.df)
}

maSimon(des=des)
