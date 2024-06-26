## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(pwb)
library(ggplot2)

## ---- set up sim--------------------------------------------------------------
# Setup
nsim = 10^6
N = 1000
p0 = 0.25
p1 = 0.5
theta = log(1-p1)-log(1-p0)
set.seed(7)
theta_U = theta_D = theta_C = ID = IU = N1 = rep(NA, nsim)

## ---- simulate, eval=FALSE----------------------------------------------------
#  # Simulate trials:
#  for(i in 1:nsim){
#  
#    trial = randomiseRPW(N, p0, p1)
#  
#    N0 = trial[1]
#    X0 = trial[2]
#  
#    N1[i] = trial[3]
#    X1 = trial[4]
#  
#    theta_U[i] = log(1 - X1/N1[i]) - log(1 - X0/N0)
#    IU[i] = 1/(1/(N1[i]-X1) + 1/(N0-X0) - 1/N1[i] - 1/N0)
#  
#    phat = (X1/N1[i] + X0/N0)/2
#  
#    v = (3 + 4*(phat-0.5))/(1 - 4*(phat-0.5))
#  
#    vT = (1-X0/N0)/(2*(1-phat))
#  
#    ID[i] = 1/(v/(N*vT*(1-vT)))
#  
#    Fhat = ID[i]/IU[i]
#  
#    theta_D[i] = log(N0/N1[i])
#  
#    theta_C[i] = (theta_U[i] - Fhat*theta_D[i])/(1-Fhat)
#  
#  }
#  
#  
#  IC = IU - ID
#  varC = 1/IC
#  varU = 1/IU
#  
#  RPWsim <- data.frame(varC=varC,
#                       varU=varU,
#                       theta_C=theta_C,
#                       theta_U=theta_U,
#                       N1=N1
#                       )

## ---- fig.height=4, fig.width=6, echo=FALSE-----------------------------------
# Load simulated data (to save time):
data(RPWsim, package = "pwb")

cond_bias = data.frame(allocation = RPWsim$N1/N,
                       uncond = (RPWsim$theta_U - theta)/theta,
                       cond = (RPWsim$theta_C - theta)/theta)


ggplot(cond_bias, aes(x=allocation, y=100*uncond)) +
  stat_summary_bin(fun='mean', bins=100, col = 'blue',
                   size=2, geom='point') +
  xlab('Allocation proportion') + ylab('Unconditional bias (%)')

ggplot(cond_bias, aes(x=allocation, y=100*cond)) +
  stat_summary_bin(fun='mean', bins=100, col = 'red',
                   size=2, geom='point') + ylim(c(-40, 40)) +
  xlab('Allocation proportion') + ylab('Conditional bias (%)')


## ---- tables, echo=FALSE------------------------------------------------------

### Tables

full.results <- rbind(
  c(mean(RPWsim$theta_C) - theta,
    mean(sqrt(RPWsim$varC))),

  c(mean(RPWsim$theta_U) - theta,
    mean(sqrt(RPWsim$varU))),

  c(weighted.mean(RPWsim$theta_C, 1/RPWsim$varC) - theta,
    weighted.mean(sqrt(RPWsim$varC), 1/RPWsim$varC)),

  c(weighted.mean(RPWsim$theta_U, 1/RPWsim$varU) - theta,
    weighted.mean(sqrt(RPWsim$varU), 1/RPWsim$varU))
)
colnames(full.results) <- c("bias", "SE")
rownames(full.results) <- c("cond",
                            "uncond",
                            "precision weighted cond",
                            "precision weighted uncond")
full.results

## ---- allocation, echo=FALSE--------------------------------------------------
under60 = (RPWsim$N1/N <= 0.60)
over60 = (RPWsim$N1/N > 0.60)
c(sum(under60), sum(over60))

## ---- tables 2, echo=FALSE----------------------------------------------------

split.results <- rbind(
  c(mean(RPWsim$theta_C[under60]) - theta,
    mean(sqrt(RPWsim$varC[under60]))),

  c(mean(RPWsim$theta_U[under60]) - theta,
    mean(sqrt(RPWsim$varU[under60]))),

  c(mean(RPWsim$theta_C[over60]) - theta,
    mean(sqrt(RPWsim$varC[over60]))),

  c(mean(RPWsim$theta_U[over60]) - theta,
    mean(sqrt(RPWsim$varU[over60])))
)
colnames(split.results) <- c("bias", "SE")
rownames(split.results) <- c("cond   (<=0.60)",
                             "uncond (<=0.60)",
                             "cond   (>0.60)",
                             "uncond (>0.60)")
split.results



