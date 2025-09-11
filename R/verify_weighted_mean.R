N=233
n1=104
theta=0.5
nsims=20
s1 <- s <- numeric(nsims)
for(i in 1:nsims){
  onetrial <- rbinom(n=N, size=1, prob=theta) # Simulate one trial of N results
  s1[i] <- cumsum(onetrial)[n1] # Number of responses at n1
  s[i] <- sum(onetrial) # Number of responses at N
}

theta.hat <- s/N
se <- sqrt(theta.hat * (1-theta.hat) / N)

# Approach 1:
theta.bar <- weighted.mean(x=theta.hat, w=1/(se)^2)
bias <- theta.bar-theta

# Approach 2:
weighted.mean(x=theta.hat-theta, w=1/(se)^2)
