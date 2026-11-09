#' pwb2stage
#'
#' For two-arm two-stage design, obtains bias for a single true response rate.
#'
#' @param N Maximum sample size
#' @param theta0 true response probability for control arm
#' @param theta1 true response probability for treatment arm
#' @param nsims number of simulations (default 1e4)
#' @return list of length 2:
#' ests: estimates
#' mc.error
#' @export



pwb2stage <- function(N, theta0, theta1, nsims){

  RD_interim = RD_final = SE_interim = SE_final = z = rep(NA, nsims)

  results <- vector("list", nsims)
  original = gsDesign(k=2, test.type = 1, beta = 0.2, sfu='OF', n.fix = 400)
  e1 = original$upper$bound[1]
  #f1 = -Inf

  for(i in 1:nsims){
    # Stage 1:
    X0 = rbinom(N/4, 1, theta0)
    X1 = rbinom(N/4, 1, theta1)

    EE = sum(X1) # Events on Experimental arm at interim
    EN = N/4 - EE  # Non-events on Experimental arm at interim

    CE = sum(X0) # Events on Control arm at interim
    CN = N/4 - CE # Non-events on Control arm at interim

    RD_interim[i] = (EE/(EE+EN)) - (CE/(CE+CN))

    SE_interim[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)

    #pval[i] = pnorm(RD_interim[i]/SE_interim[i])
    z[i] <- RD_interim[i]/SE_interim[i]

    # Stage 2:
    Y0 = rbinom(N/4, 1, theta0)
    Y1 = rbinom(N/4, 1, theta1)

    EE = sum(X1+Y1) # Events on Experimental arm at end
    EN = N/2 - EE # Non-events on Experimental arm at end

    CE = sum(X0+Y0) # Events on Control arm at end
    CN = N/2 - CE # Non-events on Control arm at end

    RD_final[i] = (EE/(EE+EN)) - (CE/(CE+CN))

    SE_final[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)

    results[[i]] <- c(sum(X0), sum(X1), sum(X0+Y0), sum(X1+Y1))

  }
  #stop_early = (pval <= 0.001)
  #continue = (pval > 0.001)
  stop_early <- z>e1
  continue <- z<=e1 # Note: No futility stopping

  results <- do.call(rbind, results)
  results <- data.frame(early.stop=stop_early,
                        RD=RD_final,
                        SE=SE_final)
  results$RD[stop_early] <- RD_interim[stop_early]
  results$SE[stop_early] <- SE_interim[stop_early]
  #names(results) <- c("C.int", "E.int", "C.final", "E.final")

  RD <- stop_early*RD_interim + continue*RD_final
  SE <- stop_early*SE_interim + continue*SE_final
  RD.stop_early <- RD[stop_early]
  RD.continue <- RD[continue]
  trueRD <- theta1-theta0

  # Table
  ntrials <- c(sum(stop_early), sum(continue), nsims, nsims)
  RD.bias <- c(mean(RD.stop_early), mean(RD.continue), mean(RD), weighted.mean(x=RD, w=1/(SE)^2)) - trueRD
# SE.RD <- c(sd(RD.stop_early), sd(RD.continue), sd(RD), sd(RD)) # By defn, same for unweighted and PW?
#  mean.SE <- c(SE[stop_early], SE[continue], mean(SE), mean(RD)) # By defn, same for unweighted and PW?

  # Table (slide 18)

  # ests <- rbind(
  #   c(sum(stop_early), mean(RD_interim[stop_early]) - trueRD,
  #     mean(SE_interim[stop_early])),
  #
  #   c(sum(continue), mean(RD_final[continue]) - trueRD,
  #     mean(SE_final[continue])),
  #
  #   c(nsims, mean(c(RD_interim[stop_early], RD_final[continue])) - trueRD,
  #     mean(c(SE_interim[stop_early],SE_final[continue]))),
  #
  #   c(nsims, weighted.mean(c(RD_interim[stop_early], RD_final[continue]),
  #                         c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)) - trueRD,
  #     weighted.mean(c(SE_interim[stop_early], SE_final[continue]),
  #                   c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)))
  # )
  all.theta.bar <- mean(RD)
  mc.error.bias <- sqrt(1/(nsims*(nsims-1)) * sum((RD-all.theta.bar)^2))
  ests <- data.frame(type=c("Stopped early", "Not stopped early", "All (unweighted)", "All (precision-weighted)"),
                     ntrials,
                     theta0,
                     theta1,
                     theta=trueRD,
                     "Bias in RD"=RD.bias,
                     bias.100=RD.bias*100,
                     mc.error.bias)
  rownames(ests) <- c("Stopped early", "Not stopped early", "All (unweighted)", "All (precision-weighted)")
  return(list(results=results,
              ests=ests,
              mc.error.bias=mc.error.bias))
}
