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
# Store interim points and corresponding stopping boundaries
diag <- curtailment::drawDiagram(ad) # diagram of design and stopping bounds
finite.bounds <- diag$bounds.mat[!is.infinite(diag$bounds.mat$success) | !is.infinite(diag$bounds.mat$fail), ]
interims <- finite.bounds$m


##### Obtain bias for range of true response probabilities ####
theta.vec <- seq(0.2, 0.8, 0.1)
summary.data <- vector("list", length(theta.vec))
raw.data <- vector("list", length(theta.vec))
set.seed(53)

for(i in 1:length(theta.vec)){
  one.run <- pwbGS(theta=theta.vec[i], des=ad, interims=interims, nsims=1e5)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
}

all.gs <- do.call(rbind, summary.data)

stop.early.count.index <- grepl("Stopped early", rownames(all.gs))
stop.early.count <- round(all.gs$nsims[stop.early.count.index]/1e5, 3)


#### Plot results ####

#### Plot bias for all theta values ####
gs.plot <- ggplot(data=all.gs, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in GS design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.15, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)
gs.plot

#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
gs.plot2 <- ggplot(data=all.gs[all.gs$type %in% c("All (naive)", "All (precision-weghted)"), ], mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in GS design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.01, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)
gs.plot2

gridExtra::grid.arrange(gs.plot, gs.plot2)


# See results for individual theta:
pwb::showTable(bias.df=all.gs, theta=0.4)
