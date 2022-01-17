## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(pwb)
library(ggplot2)

## ---- simon--------------------------------------------------------------
all.designs <- clinfun::ph2simon(pu=0.5,
                                 pa=0.6,
                                 ep1=0.05,
                                 ep2=0.1,
                                 nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
simon.des <- all.designs$out[opt.index, ]

## ---- find bias single---------------------------------------------------
# Find bias for this design, for a single true response probability theta:
  nsims <- 1e4
  single.simon <- pwbSimon(theta=0.5,
                                des=simon.des,
                                nsims=nsims)
  round(single.simon$ests[, 1:5], 5)
  round(single.simon$mc.error, 5)

## ---- find bias vector---------------------------------------------------
##### Obtain bias for range of true response probabilities ####
theta.vec <- seq(0.1, 0.9, 0.1) # Vector of true response probabilities
summary.data <- raw.data <- vector("list", length(theta.vec))
stop.early.count <- rep(NA, length(theta.vec))

for(i in 1:length(theta.vec)){
  one.run <- pwbSimon(theta=theta.vec[i],
                           des=simon.des,
                           nsims=nsims)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
  stop.early.count[i] <- one.run$ests["Stopped early", "nsims"]/nsims
}

all.simon <- do.call(rbind, summary.data)

## ---- umvue, results="hide", message=FALSE, warning=FALSE----------------

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

## ---- showtable----------------------------------------------------------
showTable(bias.df=all.simon, theta=0.4)

## ----simon plot, fig.height=4, fig.width=6, echo=FALSE, warning=FALSE----

#### Plot bias for all theta values ####
simon.colours <- setNames(c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"),
                         c("All (naive)", "All (precision-weghted)", "Complete (naive)", "Stopped early (naive)", "UMVUE"))

simon.plot <- ggplot(data=all.simon, mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.15, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)+
  scale_color_manual(values=simon.colours)
simon.plot

## ----simon plot subset, fig.height=4, fig.width=6, echo=FALSE, warning=FALSE----
#### Plot bias -- no subsetting of "stopped early" or "stopped at N" ####
simon.plot2 <- ggplot(data=all.simon[all.simon$type %in% c("All (naive)", "All (precision-weghted)"), ], mapping=aes(x=theta, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability theta",
       y="Bias (estimate - theta)",
       title="Bias in Simon design",
       subtitle="Includes proportion of trials stopped early")+
  annotate("text", x=theta.vec, y=-0.01, label=round(stop.early.count,2))+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  geom_hline(aes(yintercept=0), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=theta.vec)+
  scale_color_manual(values=simon.colours)
simon.plot2

