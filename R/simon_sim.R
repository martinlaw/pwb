# This script simulates the results of a Simon design for a vector of response rates.


# We vary response rate theta from 0.1 to 0.9

#### Find Simon designs and select optimal designs ####
all.designs <- clinfun::ph2simon(pu=0.5, pa=0.6, ep1=0.05, ep2=0.1, nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
des <- all.designs$out[opt.index, ]


##### Obtain bias for range of true response probabilities ####
theta.vec <- seq(0.1, 0.9, 0.1)
summary.data <- vector("list", length(theta.vec))
raw.data <- vector("list", length(theta.vec))
nsims <- 10000
for(i in 1:length(theta.vec)){
  one.run <- pwb::pwbSimon(theta=theta.vec[i], des=des, nsims=nsims)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
}

all.simon <- do.call(rbind, summary.data)

stop.early.count.index <- grepl("Stopped early", rownames(all.simon))
stop.early.count <- round(all.simon$nsims[stop.early.count.index]/nsims, 3)


##### UMVUE for Simon design #####
d2 <- des_gs(J=2,
             pi0=0.5,
             pi1=0.6,
             alpha=0.05,
             beta=0.1,
             Nmin=233,
             Nmax=233,
             futility=TRUE,
             efficacy=FALSE,
             optimality="null_ess",
             summary=TRUE)
est <- est_gs(des=d2, pi=theta.vec, method=c("umvue"))

umvue <- data.frame(matrix(ncol = length(names(all.simon)),
                           nrow = nrow(est$perf)))
names(umvue) <- names(all.simon)
umvue$theta <- est$perf$pi
umvue$bias <- est$perf$`Bias(hat(pi)|pi)`
umvue$type <- "UMVUE"

all.simon <- rbind(all.simon, umvue)



#### Plot bias for all theta values ####
simon.plot <- ggplot(data=all.simon, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.15, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)
simon.plot
#ggsave(filename="simon_bias.png", device="png", plot=simon.plot, width=9, height=6)


#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
simon.plot2 <- ggplot(data=all.simon[all.simon$type %in% c("All (naive)", "All (precision-weghted)"), ], mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.01, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)
simon.plot2

gridExtra::grid.arrange(simon.plot, simon.plot2)

# See results for individual theta:
showTable(bias.df=all.simon, theta=0.4)
