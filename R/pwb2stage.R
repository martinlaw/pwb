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



pwb2stage <- function(N=2428, theta0=0.2, theta1=0.15, nsims=1e4){

  RD_interim = RD_final = SE_interim = SE_final = pval = rep(NA, nsims)

  #results <- vector("list", nsims)


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

    pval[i] = pnorm(RD_interim[i]/SE_interim[i])

    # Stage 2:
    Y0 = rbinom(N/4, 1, theta0)
    Y1 = rbinom(N/4, 1, theta1)

    EE = sum(X1+Y1) # Events on Experimental arm at end
    EN = N/2 - EE # Non-events on Experimental arm at end

    CE = sum(X0+Y0) # Events on Control arm at end
    CN = N/2 - CE # Non-events on Control arm at end

    RD_final[i] = (EE/(EE+EN)) - (CE/(CE+CN))

    SE_final[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)

    #results[[i]] <- c(sum(X0), sum(X1), sum(X0+Y0), sum(X1+Y1))

  }

  stop_early = (pval <= 0.001)
  continue = (pval > 0.001)

  #results <- do.call(rbind, results)
  results <- data.frame(early.stop=stop_early,
                        RD=RD_final,
                        SE=SE_final)

  results$RD[stop_early] <- RD_interim[stop_early]
  results$SE[stop_early] <- SE_interim[stop_early]

  #names(results) <- c("C.int", "E.int", "C.final", "E.final")



  # results$early.stop <- stop_early
  # results$SE.int <- SE_interim
  # results$SE.final <- SE_final
  # results$RD.int <- RD_interim
  # results$RD.final <- RD_final
  # results$n <- rep(N/2, nsims)
  # results$n[results$early.stop] <- N/4
  # results$RD <- RD_final
  # results$RD[stop_early] <- RD_interim
  # results$SE <- SE_final
  # results$SE[stop_early] <- SE_interim

  trueRD <- theta1-theta0 # True RD no longer fixed at -0.05

  # Table (slide 18)
  ests <- rbind(
    c(sum(stop_early), mean(RD_interim[stop_early]) - trueRD,
      mean(SE_interim[stop_early])),

    c(sum(continue), mean(RD_final[continue]) - trueRD,
      mean(SE_final[continue])),

    c(nsims, mean(c(RD_interim[stop_early], RD_final[continue])) - trueRD,
      mean(c(SE_interim[stop_early],SE_final[continue]))),

    c(nsims, weighted.mean(c(RD_interim[stop_early], RD_final[continue]),
                          c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)) - trueRD,
      weighted.mean(c(SE_interim[stop_early], SE_final[continue]),
                    c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)))
  )
  ests <- data.frame(ests,
                     rep(theta0, 4),
                     rep(theta1, 4),
                     c("Stopped early", "Not stopped early", "All (unadjusted)", "All (precision-weighted)"))
  names(ests) <- c("nsims", "Bias in RD", "Mean SE of RD", "theta0", "theta1", "type")
  rownames(ests) <- c("Stopped early", "Not stopped early", "All (unadjusted)", "All (precision-weighted)")
  #signif(ests, 3)
  theta.hat <- c(RD_interim[stop_early], RD_final[continue])
  all.theta.bar <- mean(theta.hat)
  mc.error.bias <- sqrt( 1/(nsims*(nsims-1)) * sum((theta.hat-all.theta.bar)^2) )
  return(list(results=results,
              ests=ests,
              mc.error.bias=mc.error.bias))
}
