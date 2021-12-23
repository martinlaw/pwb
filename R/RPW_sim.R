########### Single trial with early stopping

# source('RPW.R')

nsim = 10^6

N = 1000

p0 = 0.25
p1 = 0.5

theta = log(1-p1)-log(1-p0)

# set.seed(7)
# 
# theta_U = theta_D = theta_C = ID = IU = N1 = rep(NA, nsim)
# 
# Sys.time()
# 
# for(i in 1:nsim){
#   
#   trial = randomiseRPW(N, p0, p1)
#   
#   N0 = trial[1]
#   X0 = trial[2]
#   
#   N1[i] = trial[3]
#   X1 = trial[4]
#   
#   theta_U[i] = log(1 - X1/N1[i]) - log(1 - X0/N0)
#   IU[i] = 1/(1/(N1[i]-X1) + 1/(N0-X0) - 1/N1[i] - 1/N0)
#   
#   phat = (X1/N1[i] + X0/N0)/2
#   
#   v = (3 + 4*(phat-0.5))/(1 - 4*(phat-0.5))
#   
#   vT = (1-X0/N0)/(2*(1-phat))
#   
#   ID[i] = 1/(v/(N*vT*(1-vT)))
#   
#   Fhat = ID[i]/IU[i]
#   
#   theta_D[i] = log(N0/N1[i])
#   
#   theta_C[i] = (theta_U[i] - Fhat*theta_D[i])/(1-Fhat)
#   
# }
# 
# Sys.time()


load("RPWsim.Rdata")


IC = IU - ID

varC = 1/IC
varU = 1/IU


### Plot

cond_bias = data.frame(allocation = N1/N, uncond = (theta_U - theta)/theta,
                       cond = (theta_C - theta)/theta)

library(ggplot2)

ggplot(cond_bias, aes(x=allocation, y=100*uncond)) +
  stat_summary_bin(fun='mean', bins=100, col = 'blue',
                   size=2, geom='point') +
  xlab('Allocation proportion') + ylab('Conditional bias (%)')

ggplot(cond_bias, aes(x=allocation, y=100*cond)) +
  stat_summary_bin(fun='mean', bins=100, col = 'red',
                   size=2, geom='point') + ylim(c(-40, 40)) + 
  xlab('Allocation proportion') + ylab('Conditional bias (%)')



### Tables

rbind(
  c(mean(theta_C) - theta,
    mean(sqrt(varC))),
  
  c(mean(theta_U) - theta,
    mean(sqrt(varU))),
  
  c(weighted.mean(theta_C, 1/varC) - theta,
    weighted.mean(sqrt(varC), 1/varC)),
  
  c(weighted.mean(theta_U, 1/varU) - theta,
    weighted.mean(sqrt(varU), 1/varU))
)


under60 = (N1/N <= 0.60)
over60 = (N1/N > 0.60)

c(sum(under60), sum(over60))

rbind(
  c(mean(theta_C[under60]) - theta,
    mean(sqrt(varC[under60]))),
  
  c(mean(theta_U[under60]) - theta,
    mean(sqrt(varU[under60]))),
  
  c(mean(theta_C[over60]) - theta,
    mean(sqrt(varC[over60]))),
  
  c(mean(theta_U[over60]) - theta,
    mean(sqrt(varU[over60]))),
  
  c(mean(theta_C) - theta,
    mean(sqrt(varC))),
  
  c(mean(theta_U) - theta,
    mean(sqrt(varU))),
  
  c(weighted.mean(theta_C, 1/varC) - theta,
    weighted.mean(sqrt(varC), 1/varC)),
  
  c(weighted.mean(theta_U, 1/varU) - theta,
    weighted.mean(sqrt(varU), 1/varU))
)

