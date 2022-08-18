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

## ----find_design--------------------------------------------------------------
all.designs <- clinfun::ph2simon(pu=0.5,
                                 pa=0.6,
                                 ep1=0.05,
                                 ep2=0.1,
                                 nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
simon.des <- all.designs$out[opt.index, ]

## ---- find bias single--------------------------------------------------------
# Find bias for this design, for a single true response probability theta:
  nsims <- 1e5
  set.seed(16)
  single.simon <- pwbSimon(theta=0.5,
                                des=simon.des,
                                nsims=nsims)
  round(single.simon$ests[, 1:5], 5)
  round(single.simon$mc.error, 5)

## ---- find bias vector--------------------------------------------------------
##### Obtain bias for range of true response probabilities ####
theta.vec <- seq(0.1, 0.9, 0.1) # Vector of true response probabilities
summary.data <- raw.data <- mc.error <- vector("list", length(theta.vec))
stop.early.count <- rep(NA, length(theta.vec))

for(i in 1:length(theta.vec)){
  one.run <- pwbSimon(theta=theta.vec[i],
                           des=simon.des,
                           nsims=nsims)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
  mc.error[[i]] <- one.run$mc.error
  stop.early.count[i] <- one.run$ests["Stopped early", "nsims"]/nsims
}

all.simon <- do.call(rbind, summary.data)

## ---- umvue, results="hide", message=F, warning=F-----------------------------

##### UMVUE for Simon design #####
d2 <- singlearm::des_gs(J=2,
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
est <- singlearm::est_gs(des=d2, pi=theta.vec, method=c("umvue"))

umvue <- data.frame(matrix(ncol = length(names(all.simon)),
                           nrow = nrow(est$perf)))
names(umvue) <- names(all.simon)
umvue$theta <- est$perf$pi
umvue$bias <- est$perf$`Bias(hat(pi)|pi)`
umvue$type <- "UMVUE"

all.simon <- rbind(all.simon, umvue)

## ---- showtable---------------------------------------------------------------
showTable(bias.df=all.simon, theta=0.4)

## ----"simon", echo=FALSE, fig.height=4, fig.width=4, warning=TRUE-------------

#### Plot bias for all theta values ####
simon.colours <- setNames(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
                         c("All (unadjusted)", "All (precision-weighted)", "Complete (unadjusted)", "Stopped early (unadjusted)", "UMVUE"))

simon.plot <- ggplot(data=all.simon, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.11, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  ylim(max(abs(all.simon$bias)) * c(-1.2, 1.2) )+
  scale_x_continuous(breaks=theta.vec)+
  scale_color_manual(values=simon.colours)+
  theme(legend.position="bottom")
simon.plot

# MC error:
mc.error.df <- do.call(rbind, mc.error)
mc.error.df <- data.frame(theta.vec, mc.error.df)
names(mc.error.df) <- c("theta", "MC error: bias", "MC error: SE")
round(mc.error.df, 5)

## ----"simon_subset", echo=FALSE, fig.height=4, fig.width=4, warning=TRUE------
#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
simon.sub.df <- all.simon[all.simon$type %in% c("All (unadjusted)", "All (precision-weighted)"), ]
simon.plot.sub <- ggplot(data=simon.sub.df, mapping=aes(x=theta, y=bias, col=type))+  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.011, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  ylim(max(abs(simon.sub.df$bias)) * c(-1.2, 1.2) )+
  scale_x_continuous(breaks=theta.vec)+
 # scale_color_manual(values=simon.colours)+
  theme(legend.position="bottom")
simon.plot.sub

## ----"simon_high_p"-----------------------------------------------------------
all.designs2 <- clinfun::ph2simon(pu=0.8,
                                 pa=0.9,
                                 ep1=0.05,
                                 ep2=0.1,
                                 nmax=150)
# Choose optimal Simon design:
opt.index <- which.min(all.designs2$out[, "EN(p0)"])
simon.des2 <- all.designs2$out[opt.index, ]

## ---- echo=FALSE--------------------------------------------------------------
##### Obtain bias for range of true response probabilities ####
theta.vec2 <- seq(0.1, 0.9, 0.05) # Vector of true response probabilities
summary.data2 <- raw.data2 <- mc.error2 <- vector("list", length(theta.vec2))
stop.early.count2 <- rep(NA, length(theta.vec2))

for(i in 1:length(theta.vec2)){
  one.run <- pwbSimon(theta=theta.vec2[i],
                           des=simon.des2,
                           nsims=nsims)
  raw.data2[[i]] <- one.run$results
  summary.data2[[i]] <- one.run$ests
  mc.error2[[i]] <- one.run$mc.error
  stop.early.count2[i] <- one.run$ests["Stopped early", "nsims"]/nsims
}

all.simon2 <- do.call(rbind, summary.data2)

## ----"simon2_subset", fig.height=4, fig.width=4, echo=FALSE, warning=T--------
#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
simon.sub.df2 <- all.simon2[all.simon2$type %in% c("All (unadjusted)", "All (precision-weighted)"), ]
#simon.sub.df2$type <- factor(simon.sub.df2$type)
#levels(simon.sub.df2$type)
simon.plot.sub2 <- ggplot(data=simon.sub.df2, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design (p0=0.8, p1=0.9)",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec2, y=-0.015, label=round(stop.early.count2,2), size=2)+
  geom_vline(aes(xintercept=0.8), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.9), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  ylim(max(abs(simon.sub.df2$bias)) * c(-1.1, 1.1) )+
  scale_x_continuous(breaks=seq(0.1, 0.9, 0.1))+
  #scale_color_manual(values=simon.colours)+
  theme(legend.position="bottom")
simon.plot.sub2

# MC error:
mc.error.df2 <- do.call(rbind, mc.error2)
mc.error.df2 <- data.frame(theta.vec2, mc.error.df2)
names(mc.error.df2) <- c("theta", "MC error: bias", "MC error: SE")
round(mc.error.df2, 5)

