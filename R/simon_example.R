# This script simulates the results of a Simon design
# for a single value of response rate theta.

# Find possible Simon designs for p0=0.5, p1=0.6, alpha=0.05, power=0.9:
all.designs <- clinfun::ph2simon(pu=0.5, pa=0.6, ep1=0.05, ep2=0.1, nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
des <- all.designs$out[opt.index, ]

# Parameters of Simon design:
n1 <- des[names(des)=="n1"]
r1 <- des[names(des)=="r1"]
N <- des[names(des)=="n"] # max sample size
r <- des[names(des)=="r"]


##### Set up and run simulation of this design for single value of response rate theta #####
theta <- 0.6
nsims <- 1e5
s1 <- s <- numeric(nsims)
states.data <- vector("list", nsims+1)

set.seed(9)
for(i in 1:nsims){
  states.data[[i]] <- DHARMa::getRandomState(NULL)
  onetrial <- rbinom(n=N, size=1, prob=theta) # Simulate one trial of N results
  s1[i] <- cumsum(onetrial)[n1] # Number of responses at n1
  s[i] <- sum(onetrial) # Number of responses at N
}
states.data[[i+1]] <- DHARMa::getRandomState(NULL)


##### Combine simulation results into data frame: #####

results <- data.frame(early.stop=s1<=r1)

results$responses <- s
results$responses[results$early.stop] <- s1[results$early.stop]

results$reject <- s>r & results$early.stop==FALSE

results$n <- rep(N, nsims)
results$n[results$early.stop] <- n1

results$complete <- results$n==N

results$theta.hat <- results$responses/results$n # naive estimate of response rate
results$se <- sqrt(results$theta.hat * (1-results$theta.hat) / results$n)


#### Naive and precision weighted bias: ####
mean(results$theta.hat)-theta
weighted.mean(x=results$theta.hat, w=1/results$se^2)-theta

