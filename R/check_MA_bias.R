library(pwb)
nsims <- 100000
p05 <- simMeta(N=200, theta=0.5, n.studies=4, nsims=nsims)

wtdmeans.05 <- rep(NA, nsims)
for(i in 1:nsims)    wtdmeans.05[i] <- weighted.mean(x=p05$theta.hat[i, ], w=1/p05$se[i, ]^2)
mean(wtdmeans.05 - 0.5)

p03 <- simMeta(N=200, theta=0.3, n.studies=4, nsims=nsims)
wtdmeans.03 <- rep(NA, nsims)
for(i in 1:nsims)    wtdmeans.03[i] <- weighted.mean(x=p03$theta.hat[i, ], w=1/p03$se[i, ]^2)
mean(wtdmeans.03 - 0.3)

p07 <- simMeta(N=200, theta=0.7, n.studies=4, nsims=nsims)
wtdmeans.07 <- rep(NA, nsims)
for(i in 1:nsims)    wtdmeans.07[i] <- weighted.mean(x=p07$theta.hat[i, ], w=1/p07$se[i, ]^2)
mean(wtdmeans.07 - 0.7)
