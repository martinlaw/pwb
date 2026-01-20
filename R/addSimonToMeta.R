# Add Simon design to MA and analyse results
addSimonToMeta <- function(nsims=1e3,
                           K=5,
                           theta,
                           theta0,
                           theta1,
                           alpha=0.05,
                           power=0.90){

# Find Simon design:
  all.designs <- clinfun::ph2simon(pu=theta0,
                                   pa=theta1,
                                   ep1=alpha,
                                   ep2=1-power,
                                   nmax=500)

# Choose optimal Simon design:
  opt.index <- which.min(all.designs$out[, "EN(p0)"])
  simon.des <- all.designs$out[opt.index, ]
  N <- simon.des["n"]
  n1 <- simon.des["n1"]

# Find bias of Simon design:
  single.simon <- pwbSimon(theta=theta,
                           des=simon.des,
                           nsims=nsims)
  simon.unwtd.bias <- single.simon$ests["All", "bias"]
  simon.pw.bias <- single.simon$ests["All (precision-weighted)", "bias"]
  print(simon.pw.bias)
  simon.mcse <- MCSEbias(single.simon$results$theta.hat)

# Simulate all other studies, find summary estimate for each replication. These studies all have sample size equal to the maximum sample size of the Simon design.
  ma.results <- simMeta(N=N,
                        theta=theta,
                        n.studies=K-1,
                        nsims=nsims)
  wtdmeans.ex <- rep(NA, nsims)
  for(i in 1:nsims) wtdmeans.ex[i] <- weighted.mean(x=ma.results$theta.hat[i, ], w=1/ma.results$se[i, ]^2)

# Add Simon design to the meta-analysis and find summary estimates again:
  combined.results <- ma.results
  combined.results$theta.hat <- cbind(single.simon$results$theta.hat, ma.results$theta.hat)
  combined.results$se <- cbind(single.simon$results$se, ma.results$se)
  wtdmeans.inc <- rep(NA, nsims)
  for(i in 1:nsims) wtdmeans.inc[i] <- weighted.mean(x=combined.results$theta.hat[i, ], w=1/combined.results$se[i, ]^2)

# Find bias in the meta-analysis, exc. and inc. Simon design:
  bias.ex <- mean(wtdmeans.ex-theta)
  unwtd.bias.ex <- mean(rowMeans(ma.results$theta.hat)) - theta # Bias using simple mean of estimates (i.e. without correct weighting)
  bias.inc <- mean(wtdmeans.inc-theta)
  unwtd.bias.inc <- mean(rowMeans(combined.results$theta.hat)) - theta

  single.simon$results$test.wt <- 1
  early.stop.wt <- (((K-1)*N)+n1)/(K*N) # Relative total sample size of MA when early stop / no early stop
  single.simon$results$test.wt[single.simon$results$early.stop==T] <- early.stop.wt
  bias.inc.test <- weighted.mean(x=wtdmeans.inc, w=single.simon$results$test.wt)-theta

# MCSE:
  mcse.ex <- MCSEbias(wtdmeans.ex)
  mcse.inc <- MCSEbias(wtdmeans.inc)

# Bias, PWB, MCSE of each non-AD study:
  bias.nonAD <- colMeans(ma.results$theta.hat)-theta
  pwb.nonAD <- rep(NA, ncol(ma.results$theta.hat))
  for(i in 1:ncol(ma.results$theta.hat)){
    pwb.nonAD[i] <- weighted.mean(x=ma.results$theta.hat[, i], w=1/ma.results$se[,i]^2) - theta
  }
  mean.pwb.nonAD <- mean(pwb.nonAD)
  mcse.nonAD <- apply(ma.results$theta.hat, 2, MCSEbias)

# Precision
  emp.SE <- apply(combined.results$theta.hat, 2, empSE) # emp SE of each column
  precision <- 1/(emp.SE^2)
  prec.simon <- precision[1]  # Precision of Simon
  # Precision of Simon design as fraction of precision of all designs:
  simon.frac <- prec.simon/sum(precision)
  ave.prec.non.ad <- mean(precision[-1]) # Mean precision of non-ADs
# Fraction of Simon vs mean of non-adaptive:
  prec.frac <- prec.simon/ave.prec.non.ad

# Results: Long table
bias.studies.unw <- c(bias.nonAD, mean(bias.nonAD), simon.unwtd.bias)
bias.studies.pwb <- c(pwb.nonAD, mean(pwb.nonAD), simon.pw.bias)
studies.MCSE <- c(mcse.nonAD, NA, simon.mcse)
studies.rown <- rep(c(paste("Study", 1:(K-1)), paste0("Mean 1-", K-1), "Simon"))
sample.size <- c(rep(N, K-1), NA, mean(single.simon$results$n))
#browser()
smaller.table <- data.frame(unw=bias.studies.unw,
                            pwb=bias.studies.pwb,
                            mcse=studies.MCSE,
                            Nbar=sample.size,
                            row.names=studies.rown)
wide.results <- data.frame(theta=theta,
                           mean.pwb=mean.pwb.nonAD,
                           bias.ex=bias.ex,
                           mcse.ex=mcse.ex,
                           bias.inc=bias.inc,
                           mcse.inc=mcse.inc,
                           change=bias.inc-bias.ex,
                           W=simon.frac,
                           simon.pw.bias=simon.pw.bias,
                           pred.change=simon.frac*(simon.pw.bias-bias.ex))
# Results: output
results <- list(full.ma=smaller.table,
                wide=wide.results,
                theta=theta,
                mean.pwb=mean.pwb.nonAD,
                bias.ex=bias.ex,
                mcse.ex=mcse.ex,
                bias.inc=bias.inc,
                mcse.inc=mcse.inc,
                change=bias.inc-bias.ex,
                W=simon.frac,
                simon.pw.bias=simon.pw.bias,
                pred.change=simon.frac*(simon.pw.bias-bias.ex),
                prec.frac=prec.frac)
return(results)
}
