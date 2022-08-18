## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  dev='png',
  fig.path="figures/")

## ----setup--------------------------------------------------------------------
library(pwb)
library(ggplot2)

## ----"GS_bounds", fig.height=4, fig.width=6, echo=FALSE, warning=FALSE--------
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
# Obtain stopping bounds:
  bounds <- curtailment::drawDiagram(ad)
  bounds$diagram

## ----GS_single----------------------------------------------------------------
# Find bias for this design, for a single true response probability theta:
  nsims <- 1e5
  single.GS <- pwbGS(theta=0.5,
                     des=ad,
                     bounds=bounds$bounds.mat,
                     nsims=nsims)
  single.GS$ests[, 1:5]
  round(single.GS$mc.error, 5)

## ---- multiple GS-------------------------------------------------------------
theta.vec <- seq(0.2, 0.8, 0.1) # Vector of true response probabilities
summary.data <- raw.data <- vector("list", length(theta.vec))
stop.early.count <- rep(NA, length(theta.vec))
set.seed(53)

for(i in 1:length(theta.vec)){
  one.run <- pwbGS(theta=theta.vec[i], des=ad, bounds=bounds$bounds.mat, nsims=nsims)
  summary.data[[i]] <- one.run$ests
  stop.early.count[i] <- one.run$ests["Stopped early", "nsims"]/nsims
}

all.gs <- do.call(rbind, summary.data)

## ----"GS", fig.height=4, fig.width=4, echo=FALSE, warning=FALSE---------------
#### Plot bias for all theta values ####
type.colours <- setNames(c("#F8766D", "#00BFC4", "#7CAE00", "#C77CFF"),
                         c("All (unadjusted)", "All (precision-weghted)", "Complete (unadjusted)", "Stopped early (unadjusted)"))
gs.plot <- ggplot(data=all.gs, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in GS design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.15, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)+
  scale_color_manual(values=type.colours)+
  theme(legend.position = 'bottom')
gs.plot

## ----"GS_subset", fig.height=4, fig.width=4, echo=FALSE-----------------------
#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
gs.plot2 <- ggplot(data=all.gs[all.gs$type %in% c("All (unadjusted)", "All (precision-weghted)"), ],
                   mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in GS design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.01, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  #scale_color_manual(values=type.colours)+
  scale_x_continuous(breaks=theta.vec)+
  theme(legend.position = 'bottom')
gs.plot2

