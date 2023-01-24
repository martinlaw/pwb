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



all.designs <- clinfun::ph2simon(pu=0.5,
                                 pa=0.6,
                                 ep1=0.05,
                                 ep2=0.1,
                                 nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
simon.des <- all.designs$out[opt.index, ]

theta <- 0.3
set.seed(16)
single.simon2 <- pwbSimon(theta=0.3,
                          des=simon.des,
                          nsims=nsims)
round(single.simon2$ests[3:4, 2:4], 5)
#>                              bias mean.SE  emp.SE
#> All                       0.00014 0.04468 0.04491
#> All (precision-weighted) -0.00390 0.04450      NA
mcse.simon2 <- MCSEbias(single.simon2$results$theta.hat)
mcse.simon2 * 1e5
MCSE(single.simon2$results$theta.hat - theta)
mean(single.simon2$results$se)

head(single.simon2)



