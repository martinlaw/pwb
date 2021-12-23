# This script simulates the results of a Simon design for a vector of response rates.
# Uses function pwbSimon for finding precision weighted bias:
source("pwbSimon.R")

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
for(i in 1:length(theta.vec)){
  one.run <- pwbSimon(theta=theta.vec[i], des=des, nsims=1e5)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
}

all.simon <- do.call(rbind, summary.data)

stop.early.count.index <- grepl("Stopped early", rownames(all.simon))
stop.early.count <- round(all.simon$nsims[stop.early.count.index]/1e5, 2)


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
       title="Bias in Simon design -- precision-weighted vs naive",
       subtitle="Includes naive estimates for early stopped vs complete trials and also prop'n stopped early")+
  annotate("text", x=theta.vec, y=-0.15, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=seq(0,1,0.1))
simon.plot
#ggsave(filename="simon_bias.png", device="png", plot=simon.plot, width=9, height=6)


# showTable: shows results for single value of theta
#
# Arguments:
#
# theta (required)
# latex: show Latex code for table (TRUE/FALSE)
# digit: number of rounding digits for bias and SE
#
showTable <- function(theta, latex=FALSE, digit=4){
  single.theta.df <- all.simon[abs(all.simon$theta-theta)<0.001, ]
  rownames(single.theta.df) <- single.theta.df$type
  if(latex==TRUE){
    print(xtable(single.theta.df[, 1:4],
               caption=paste("Bias, mean SE and emp SE for Simon design when theta=", theta, sep=""),
               align=rep("r", 5),
               digits=c(0, 0, digit, digit, digit),
               label="tab:simon"),
        include.rownames=TRUE,
        booktabs=TRUE,
        sanitize.text.function=function(x){x})
  }
  single.theta.df$bias <- round(single.theta.df$bias, digit)
  single.theta.df$mean.SE <- round(single.theta.df$mean.SE, digit)
  single.theta.df$emp.SE <- round(single.theta.df$emp.SE, digit)
  single.theta.df
}
showTable(theta=0.4)
