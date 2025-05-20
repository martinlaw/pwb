pwbSimonPlusNonAD <- function (theta, des, nsims = 1e3)
{
  n1 <- des[names(des) == "n1"]
  r1 <- des[names(des) == "r1"]
  N <- des[names(des) == "n"]
  r <- des[names(des) == "r"]
  reject <- s1 <- s <- numeric(nsims)
  for (i in 1:nsims) {
    onetrial <- rbinom(n = N, size = 1, prob = theta)
    s1[i] <- cumsum(onetrial)[n1]
    s[i] <- sum(onetrial)
  }
  results <- data.frame(early.stop = s1 <= r1)
  results$responses <- s
  results$responses.all <- s
  results$responses[results$early.stop] <- s1[results$early.stop]
  results$reject <- s > r & results$early.stop == FALSE
  results$n <- rep(N, nsims)
  results$n[results$early.stop] <- n1
  results$complete <- results$n == N
  results$theta.hat <- results$responses/results$n
  results$se <- sqrt(results$theta.hat * (1 - results$theta.hat)/results$n)
  results$responses.cor <- results$responses
  results$responses.cor[results$responses.cor == 0] <- 0.5
  results$responses.cor[results$responses.cor == results$n] <- results$responses.cor[results$responses.cor == results$n]-0.5
  results$theta.hat.cor <- results$responses.cor/results$n
  results$se.cor <- sqrt(results$theta.hat.cor * (1 - results$theta.hat.cor)/results$n)
  early.stop <- results[results$complete == FALSE, ]
  stop.at.N <- results[results$complete == TRUE, ]
  early.stop.theta.bar <- mean(early.stop$theta.hat.cor)
  stop.at.N.theta.bar <- mean(stop.at.N$theta.hat.cor)
  all.theta.bar <- mean(results$theta.hat.cor)
  early.stop.emp.SE <- sqrt((1/(nrow(early.stop) - 1)) * sum((early.stop$theta.hat.cor -
                                                                early.stop.theta.bar)^2))
  stop.at.N.emp.SE <- sqrt((1/(nrow(stop.at.N) - 1)) * sum((stop.at.N$theta.hat.cor -
                                                              stop.at.N.theta.bar)^2))
  all.emp.SE <- sqrt((1/(nsims - 1)) * sum((results$theta.hat.cor - all.theta.bar)^2))
  emp.SE <- c(early.stop.emp.SE, stop.at.N.emp.SE, all.emp.SE,
              NA)
  mean.SE <- c(mean(early.stop$se.cor), mean(stop.at.N$se.cor),
               mean(results$se.cor), weighted.mean(results$se.cor, w = 1/(results$se.cor)^2))
  theta.bar.vec <- c(early.stop.theta.bar, stop.at.N.theta.bar,
                     all.theta.bar, weighted.mean(x = results$theta.hat.cor,
                                                  w = 1/(results$se.cor)^2))
  bias <- theta.bar.vec - theta
  trials <- c(nrow(early.stop), nrow(stop.at.N), nsims, nsims)



# Add fixed/non-AD measures:
results$theta.hat.nonAD <- results$responses.all/N
# In case of zero successes or N successes:
results$responses.all.cor <- results$responses.all
results$responses.all.cor[results$responses.all.cor == 0] <- 0.5
results$responses.all.cor[results$responses.all.cor == N] <- results$responses.all.cor[results$responses.all.cor == N] - 0.5
results$theta.hat.nonAD.cor <- results$responses.all.cor/N
results$se.nonAD.cor <- sqrt(results$theta.hat.nonAD.cor * (1 - results$theta.hat.nonAD.cor)/N)

theta.bar.vec.nonAD <- c(mean(results$theta.hat.nonAD.cor),
                         weighted.mean(x = results$theta.hat.nonAD.cor,
                                       w = 1/(results$se.nonAD.cor)^2))
bias.nonAD <- theta.bar.vec.nonAD - theta
all.theta.bar.nonAD <-  mean(results$theta.hat.nonAD.cor)
mean.SE.nonAD <- c(mean(results$se.nonAD.cor), weighted.mean(results$se.nonAD.cor, w = 1/(results$se.nonAD.cor)^2))
emp.SE.nonAD <- sqrt((1/(nsims - 1)) * sum((results$theta.hat.nonAD.cor - all.theta.bar.nonAD)^2))

 ests <- data.frame(nsims,
                    bias=c(bias[3:4], bias.nonAD),
                    mean.SE=c(mean.SE[3:4], mean.SE.nonAD),
                    emp.SE=c(emp.SE[3:4], emp.SE.nonAD, NA),
                    theta,
                    Bias=c("Unconditional", "Precision-weighted"),
                    Design=c("Simon", "Simon", "Fixed", "Fixed"))
  row.names(ests) <- c("Simon uncond", "Simon PW", "Fixed uncond", "Fixed PW")
  ests$Bias <- factor(ests$Bias)
  ests$Design <- factor(ests$Design)
  ests$early.stop <- sum(results$early.stop)/nsims

mc.error.bias.all <- sqrt(sum((results$theta.hat.cor - all.theta.bar)^2)/(nsims*(nsims - 1)))
mc.error.se.all <- all.emp.SE/sqrt(2 * (nsims - 1))
mc.error.bias.all.nonAD <- sqrt(sum((results$theta.hat.nonAD.cor - all.theta.bar.nonAD)^2)/(nsims*(nsims - 1)))
mc.error.se.all.nonAD <- emp.SE.nonAD/sqrt(2 * (nsims - 1))


# Summary:
  all.data <- list(results = results,
                   ests = ests,
                   mc.error = c(mc.error.bias.all,
                                mc.error.se.all,
                                mc.error.bias.all.nonAD,
                                mc.error.se.all.nonAD))
  return(all.data)
}

# single.simon <- pwbSimonPlusNonAD(theta=0.5,
#                          des=simon.des,
#                          nsims=1e3)
# head(single.simon$results)
# single.simon$ests
# single.simon$mc.error

