#' maGS
#'
#' Examines bias for a single group sequential design WRT a meta-analysis of non-adaptive designs
#'
#' @param theta true response probablity
#' @param des a realisation of a group sequential design, formatted as per
#'  output of curtailment::singlearmDesign
#' @param bounds interim analysis bounds, also from curtailment::singlearmDesign
#' @param nsims number of simulations (default 1e5)
#' @param n.studies number of non-adaptive studies in meta-analysis
#' @return Data frame of results
#' @export
#'


maGS <- function(theta=0.5, des, bounds, nsims=1e5, n.studies=4){



  # Obtain bias for single true response probability:
  gs.data <- pwbGS(theta=theta,
                   des=des,
                   bounds=bounds$bounds.mat,
                   nsims=nsims)$results

  #### Meta-analysis: Simulate further trials  ####
  N <- des$all.des[colnames(des$all.des)=="n"] # max sample size of Simon design
  responses <- vector("list", nsims)

  for(i in 1:nsims){
    responses[[i]] <- rbinom(n=n.studies, size=N, prob=theta)
  }
  responses <- do.call(rbind, responses)
  theta.hat <- responses/N
  se <- sqrt((theta.hat)*(1-theta.hat)/N)

  # Combine results from GS simulation and MA simulation:
  theta.hat <- cbind(gs.data$theta.hat, theta.hat)
  se <- cbind(gs.data$se, se)


  # Summary estimate and bias for all trials:
  wtdmeans.all <- rep(NA, nsims)
  for(i in 1:nrow(theta.hat)){
    wtdmeans.all[i] <- weighted.mean(x=theta.hat[i, ], w=1/se[i, ]^2)
  }
  bias.all <- wtdmeans.all-theta


  # Summary estimate and bias, excluding adaptive design if stopped early:
  wtdmeans.exclude <- wtdmeans.all
  for(i in (1:nsims)[gs.data$complete==FALSE]){
    wtdmeans.exclude[i] <- weighted.mean(x=theta.hat[i, 2:(n.studies+1)], w=1/se[i, 2:(n.studies+1)]^2)
  }
  bias.exclude <- wtdmeans.exclude-theta

  mean.trials <- n.studies+mean(gs.data$complete==TRUE)

  # Present results #
  ma.df <- data.frame(reps=rep(nsims, 2),
                      mean.studies=c(n.studies+1, mean.trials),
                      bias=c(mean(bias.all), mean(bias.exclude)),
                      mean.se=c(sd(wtdmeans.all), sd(wtdmeans.exclude)),
                      theta=rep(theta, 2),
                      type=c("All trials", "Exclude early stopped")
  )
  row.names(ma.df) <- c("All trials", "Exclude early stopped")
  return(ma.df)
}

