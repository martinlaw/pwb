# This script simulates the results of single-arm multi-stage with binary
# outcome, for a single value of response rate theta.
#
# Find a multi-stage design for p0=0.5, p1=0.6, alpha=0.1, power=0.8 and
# satisfying various other parameters:

ad <- curtailment::singlearmDesign(nmin=50,
                                   nmax=140,
                                   C=20,
                                   minstop=40,
                                   p0=0.5,
                                   p1=0.6,
                                   alpha=0.1,
                                   power=0.8,
                                   minthetaE=1,
                                   maxthetaF=0,
                                   max.combns=1e3)
diag <- drawDiagram(ad) # Diagram and stopping bounds


##### Set up and run simulation of this design for single value of response rate theta #####

N <- ad$all.des[, "n"]
theta <- 0.7 # equal to p1
nsims <- 1e5
interims <- seq(from=80, to=140, by=20)
finite.bounds <- diag$bounds.mat[!is.infinite(diag$bounds.mat$success) | !is.infinite(diag$bounds.mat$fail), ]
results <- vector("list", nsims)
states.data <- vector("list", nsims+1)
set.seed(53)

for(i in 1:nsims){
    states.data[[i]] <- DHARMa::getRandomState(NULL)
    onetrial <- rbinom(n=N, size=1, prob=theta) # Simulate one trial of N results
    s.list <- cumsum(onetrial)[interims] # Number of responses at each interim
    fail.stage <- match(TRUE, s.list <= finite.bounds$fail, nomatch=100) # Earliest crossing of futility bound
    success.stage <- match(TRUE, s.list >= finite.bounds$success, nomatch=100)  # Earliest crossing of efficacy bound
    reject <- success.stage<fail.stage
    stop.stage <- min(fail.stage, success.stage)
    final.n <- interims[stop.stage] # Trial sample size
    final.s <- s.list[stop.stage] # Trial number of responses
    results[[i]] <- c(final.s, final.n, reject)
}
states.data[[i+1]] <- DHARMa::getRandomState(NULL)


##### Combine simulation results into data frame: #####

results <- do.call(rbind, results)
colnames(results) <-  c("s", "n", "reject")
results <- data.frame(results)
results$complete <- results$n==140 # Trial continued to final stage
results$theta.hat <- results$s/results$n
results$se <- sqrt(results$theta.hat * (1-results$theta.hat) / results$n)


##### Analyse results #####

early.stop <- results[results$complete==FALSE, ]
stop.at.N <- results[results$complete==TRUE, ]

# Estimates of response rate:
early.stop.theta.bar <- mean(early.stop$theta.hat)
stop.at.N.theta.bar <- mean(stop.at.N$theta.hat)
all.theta.bar <- mean(results$theta.hat)

# Empirical SEs:
early.stop.emp.SE <- sqrt((1/(nrow(early.stop)-1)) * sum((early.stop$theta.hat-early.stop.theta.bar)^2))
stop.at.N.emp.SE <- sqrt((1/(nrow(stop.at.N)-1)) * sum((stop.at.N$theta.hat-stop.at.N.theta.bar)^2))
all.emp.SE <- sqrt((1/(nsims-1)) * sum((results$theta.hat-all.theta.bar)^2))
emp.SE <- c(early.stop.emp.SE, stop.at.N.emp.SE, all.emp.SE, NA)

mean.SE <- c(mean(results$se[results$complete==FALSE]),
             mean(results$se[results$complete==TRUE]),
             mean(results$se),
             weighted.mean(results$se, w=1/(results$se)^2)
             )

theta.bar.vec <- c(early.stop.theta.bar,
                   stop.at.N.theta.bar,
                   all.theta.bar,
                   weighted.mean(x=results$theta.hat, w=1/(results$se)^2)
                   )
bias <- theta.bar.vec-theta

trials <- c(nrow(early.stop),
            nrow(stop.at.N),
            nsims,
            nsims
            )

ests <- data.frame(nsims=trials,
                   bias=round(bias, 4),
                   mean.SE=round(mean.SE, 4),
                   emp.SE=round(emp.SE, 4),
                   theta)
row.names(ests) <- c("Stopped early", "Complete", "All", "All (precision-weghted)")
ests$type <- c("Stopped early (naive)", "Complete (naive)", "All (naive)", "All (precision-weghted)")
ests

mc.error.bias.all <- sqrt( 1/(nsims*(nsims-1)) * sum((results$theta.hat-all.theta.bar)^2) )
mc.error.bias.all
mc.error.se.all <- all.emp.SE/sqrt(2*(nsims-1))
mc.error.se.all


