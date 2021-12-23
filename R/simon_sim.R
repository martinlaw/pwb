# This script contains a function that simulates the results of a Simon design
# for a vector of response rates


# Parameters that may be changed: alpha, power, p0, p1, p
# Fix alpha, power, p1-p0
# Vary p0, p1: p0=0.1, 0.2, ..., 0.7. p1=p0+0.2
# Vary p from 0.1 to 0.9

x <- clinfun::ph2simon(pu=0.5, pa=0.6, ep1=0.05, ep2=0.1, nmax=500)

opt.index <- which.min(x$out[,"EN(p0)"])
opt <- x$out[opt.index, ]

# function to examine bias for a single Simon design
pwbSimon <- function(theta, des, nsims=1e5){
  # p: true response probablity
  # des: a realisation of a Simon design, formatted as per output of clinfun::ph2simon
  
   ##### Simulate many trials of this design: #####
  n1 <- des[names(des)=="n1"]
  r1 <- des[names(des)=="r1"]
  N <- des[names(des)=="n"] # max sample size
  r <- des[names(des)=="r"]
  reject <- s1 <- s <- numeric(nsims)
  states.data <- vector("list", nsims+1)
  
  for(i in 1:nsims){
    states.data[[i]] <- DHARMa::getRandomState(NULL)
    onetrial <- rbinom(n=N, size=1, prob=theta) # Simulate one trial of N results
    s1[i] <- cumsum(onetrial)[n1] # Number of responses at n1
    s[i] <- sum(onetrial) # Number of responses at N
  }
states.data[[i+1]] <- DHARMa::getRandomState(NULL)
  
  ##### Combine results into data frame: #####
  results <- data.frame(early.stop=s1<=r1)
  results$responses <- s
  results$responses[results$early.stop] <- s1[results$early.stop]
  results$reject <- s>r & results$early.stop==FALSE
  results$n <- rep(N, nsims)
  results$n[results$early.stop] <- n1
  results$complete <- results$n==N
  
  results$theta.hat <- results$responses/results$n
  results$se <- sqrt(results$theta.hat * (1-results$theta.hat) / results$n)
  
  # print(paste("p=", p, ": number of zero responses and all responses", sep=""))
  # print(sum(results$responses==0))
  # print(sum(results$responses==results$n))

  # Add correction for instances of zero responses and all responses:
  results$responses.cor <- results$responses
  results$responses.cor[results$responses.cor==0] <- 0.5 # Add 0.5 when zero responses.
  # Subtract 0.5 when p-hat=1.0:
  results$responses.cor[results$responses.cor==results$n] <- results$responses.cor[results$responses.cor==results$n]-0.5 
  # Calculate theta hat then SE using corrected responses:
  results$theta.hat.cor <- results$responses.cor/results$n
  results$se.cor <- sqrt(results$theta.hat.cor * (1-results$theta.hat.cor) / results$n)

  # Unweighted estimates of response rate:
  #summary(results$theta.hat)
  early.stop <- results[results$complete==FALSE, ]
  #summary(early.stop$theta.hat)
  stop.at.N <- results[results$complete==TRUE, ]
  #summary(stop.at.N$theta.hat)
  #early.stop.reject <- early.stop[early.stop$reject==1, ]
  #early.stop.accept <- early.stop[early.stop$reject==0, ]
  #summary(early.stop.reject$theta.hat)
  #summary(early.stop.accept$theta.hat)
  
  # Precision-weighted estimates of response rate:
  # weighted.mean(x=results$theta.hat, w=1/(results$se.cor)^2)
  # weighted.mean(x=early.stop$theta.hat, w=1/(early.stop$se.cor)^2)
  # weighted.mean(x=stop.at.N$theta.hat, w=1/(stop.at.N$se.cor)^2)
  #weighted.mean(x=early.stop.reject$theta.hat, w=1/(early.stop.reject$se.cor)^2)
  #weighted.mean(x=early.stop.accept$theta.hat, w=1/(early.stop.accept$se.cor)^2)
  
  early.stop.theta.bar <- mean(early.stop$theta.hat)
  stop.at.N.theta.bar <- mean(stop.at.N$theta.hat)
  all.theta.bar <- mean(results$theta.hat)
  
  early.stop.emp.SE <- sqrt((1/(nrow(early.stop)-1)) * sum((early.stop$theta.hat-early.stop.theta.bar)^2))
  stop.at.N.emp.SE <- sqrt((1/(nrow(stop.at.N)-1)) * sum((stop.at.N$theta.hat-stop.at.N.theta.bar)^2))
  all.emp.SE <- sqrt((1/(nsims-1)) * sum((results$theta.hat-all.theta.bar)^2))
  
  emp.SE <- c(early.stop.emp.SE, stop.at.N.emp.SE, all.emp.SE, NA)
  
  mean.SE <- c(mean(early.stop$se.cor),
               mean(stop.at.N$se.cor),
               mean(results$se.cor),
               weighted.mean(results$se.cor, w=1/(results$se.cor)^2)
  )
  
  
  theta.bar.vec <- c(early.stop.theta.bar,
                     stop.at.N.theta.bar,
                     all.theta.bar,
                     weighted.mean(x=results$theta.hat, w=1/(results$se.cor)^2))
  bias <- theta.bar.vec-theta
  
  trials <- c(nrow(early.stop),
              nrow(stop.at.N),
              nsims,
              nsims)
  
  ests <- data.frame(nsims=trials,
                     bias,
                     mean.SE,
                     emp.SE,
                     theta)
  row.names(ests) <- c("Stopped early", "Complete", "All", "All (precision-weghted)")
  ests$type <- c("Stopped early (naive)", "Complete (naive)", "All (naive)", "All (precision-weghted)")
  
  mc.error.bias.all <- sqrt( 1/(nsims*(nsims-1)) * sum((results$theta.hat-all.theta.bar)^2) )
  mc.error.se.all <- all.emp.SE/sqrt(2*(nsims-1))
  
  all.data <- list(results=results,
                   ests=ests)
  return(all.data)
} # END OF FN

##### Obtain bias for range of true response probabilities: ####
p.vec <- seq(0.1, 0.9, 0.1)
summary.data <- vector("list", length(p.vec))
raw.data <- vector("list", length(p.vec))
for(i in 1:length(p.vec)){
  one.run <- pwbSimon(p=p.vec[i], des=opt, nsims=1e5)
  raw.data[[i]] <- one.run$results
  summary.data[[i]] <- one.run$ests
}

all.simon <- do.call(rbind, summary.data)

stop.early.count.index <- grepl("Stopped early", rownames(all.simon))
stop.early.count <- round(all.simon$nsims[stop.early.count.index]/1e5, 2)


##### UMVUE for Simon design:
library(singlearm)
d2 <- des_gs(J=2, pi0=0.5, pi1=0.6, alpha=0.05, beta=0.1, Nmin=233, Nmax=233, futility=TRUE, efficacy=F, optimality="null_ess", summary=T)
est <- est_gs(des=d2, pi=p.vec, method=c("umvue"))
est.test <- est_gs(des=d2, pi=0.6, method=c("naive"))
est.test$perf
est$perf$`Var(hat(pi)|pi)`

umvue <- data.frame(matrix(ncol = length(names(all.simon)), nrow = nrow(est$perf)))
names(umvue) <- names(all.simon)
umvue$p <- est$perf$pi
umvue$bias <- est$perf$`Bias(hat(pi)|pi)`  
umvue$type <- "UMVUE"

all.simon <- rbind(all.simon, umvue)

simon.plot <- ggplot(data=all.simon, mapping=aes(x=p, y=bias, col=type))+
  geom_line(alpha=0.4, size=1)+
  geom_point()+
  labs(x="True response probability p",
       y="Bias (estimate - p)",
       title="Bias in Simon design -- precision-weighted vs naive",
       subtitle="Includes naive estimates for early stopped vs complete trials and also prop'n stopped early")+
  annotate("text", x=p.vec, y=-0.15, label=stop.early.count)+
  geom_vline(aes(xintercept=0.5), col="grey", linetype="dashed")+
  geom_vline(aes(xintercept=0.6), col="grey", linetype="dashed")+
  scale_x_continuous(breaks=seq(0,1,0.1))
simon.plot
#ggsave(filename="simon_bias.png", device="png", plot=simon.plot, width=9, height=6)

# Why is PWB worse? Possible reasons:
# (1) Problem is caused by corrections when there are zero responses or all responses.
# This would happen more at the extreme values of p, and this is indeed where
# we see the biggest difference in bias.
# (2) Mistake in code: when all trials end early and when all trials stop at N,
# precision-weighted bias should be identical to naive estimate. Check this,
# bearing in mind the fact that (1) above could be causing the difference in bias.



# ggplot(results, aes(x=final.n, y=responses/final.n)) +
#   geom_jitter(width=0.45, height=0.015, shape=".", alpha=0.1) +
#   #geom_point(size=0.1, alpha=0.01)+
#   xlab("Sample size") + ylab("Naive estimate of response rate") +
#   geom_hline(yintercept=0.3, color="red") +
#   theme_bw()
# 
# ggplot(results, aes(x=responses/final.n)) +
#   geom_density(adjust=1.25) +
#   xlab("Naive estimate of response rate") +
#   geom_vline(xintercept=0.3, color="red") +
#   theme_bw()




library(xtable)

point <- all.simon[all.simon$p==0.6, ]
rownames(point) <- point$type
print(xtable(point[, 1:4],
             caption="Bias, mean SE and emp SE for Simon design when p=0.6", 
             align=rep("r", 5),
             digits=c(0, 0, 3, 3, 3),
             label="tab:simon"),
      include.rownames=TRUE, 
      booktabs=TRUE, 
      file="simon_tab.tex",
      sanitize.text.function=function(x){x})



sum(pmf_i$`f(s,m|pi)`*(est_i$`hat(pi)(s,m)`)^2) - perf$`E(hat(pi)|pi)`[i]^2

responses <- rbinom(1e6, 104, 0.5)
PET <- sum(responses <= 54)/1e6
earlystop <- responses[responses<=54]
earlystop.p <- earlystop/104
cont <- responses[responses>54]
stage2s <- rbinom(length(cont), 233-104, 0.5)
length(cont)
length(stage2s)
final.s <- cont + stage2s
final.p <- final.s/233
all.p.hat <- c(earlystop.p, final.p)
summary(all.p.hat)
