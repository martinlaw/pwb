# usethis::use_package("ggplot2")
# usethis::use_package("gridExtra")
# usethis::use_package("clinfun")
# usethis::use_package("xtable")
# usethis::use_package("DHARMa")
# usethis::use_package("librarian")
# usethis::use_dev_package("curtailment", remote="martinlaw/curtailment")
# usethis::use_dev_package("singlearm", remote="mjg211/singlearm")

# What happens when you weight the estimates only?
nsims <- 1e5
N <- 1000
theta <- 0.3

responses <- rbinom(nsims, N, theta)
theta.hat <- responses/N
naive <- mean(theta.hat-theta)
se <- sqrt(theta.hat*(1-theta.hat)/N)
pwbias <- weighted.mean(x=theta.hat-theta, w=1/se^2)
pw.est <- weighted.mean(x=theta.hat, w=1/se^2)-theta

######## Ian's design but intervention response rate is allowed to vary:

########### Single trial with early stopping


twoArmTwoStage <- function(N=2428, nsim=1e4, p1=0.15){

  RD_interim = RD_final = SE_interim = SE_final = pval = rep(NA, N)

  for(i in 1:nsim){
    # Stage 1:
    X0 = rbinom(N/4, 1, 0.2)
    X1 = rbinom(N/4, 1, p1)

    EE = sum(X1) # Events on Experimental arm at interim
    EN = N/4 - EE  # Non-events on Experimental arm at interim

    CE = sum(X0) # Events on Control arm at interim
    CN = N/4 - CE # Non-events on Control arm at interim

    RD_interim[i] = (EE/(EE+EN)) - (CE/(CE+CN))

    SE_interim[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)

    pval[i] = pnorm(RD_interim[i]/SE_interim[i])

    # Stage 2:
    Y0 = rbinom(N/4, 1, 0.2)
    Y1 = rbinom(N/4, 1, p1)

    EE = sum(X1+Y1) # Events on Experimental arm at end
    EN = N/2 - EE # Non-events on Experimental arm at end

    CE = sum(X0+Y0) # Events on Control arm at end
    CN = N/2 - CE # Non-events on Control arm at end

    RD_final[i] = (EE/(EE+EN)) - (CE/(CE+CN))

    SE_final[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)

  }

  stop_early = (pval <= 0.001)
  continue = (pval > 0.001)

  trueRD <- 0.2-p1 # True RD no longer fixed at 0.05

  # Table (slide 18)
  one.trial <- rbind(
    c(sum(stop_early), mean(RD_interim[stop_early]) + trueRD,
      mean(SE_interim[stop_early])),

    c(sum(continue), mean(RD_final[continue]) + trueRD,
      mean(SE_final[continue])),

    c(nsim, mean(c(RD_interim[stop_early], RD_final[continue])) + trueRD,
      mean(c(SE_interim[stop_early],SE_final[continue]))),

    c(nsim, weighted.mean(c(RD_interim[stop_early], RD_final[continue]),
                          c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)) + trueRD,
      weighted.mean(c(SE_interim[stop_early], SE_final[continue]),
                    c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)))
  )
  colnames(one.trial) <- c("nsims", "Bias in RD", "Mean SE of RD")
  rownames(one.trial) <- c("Stopped early", "Not stopped early", "All", "All, precision weighted")
  signif(one.trial, 3)
}
twoArmTwoStage(p1=0.1)

