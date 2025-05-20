sims <- 1000
prop <- seq(0,1,0.01)
n1 <- 104
n <- 233
V <- NULL

for(i in 1:101){
  sim1 <- rep(c(n1,n), sims*c(prop[i], 1-prop[i]))
  V[i] <- var(sim1)
}

plot(prop, V)
prop[which.max(V)]
