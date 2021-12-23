########### Single trial with early stopping

nsim = 10^4

N = 2428

RD_interim = RD_final = SE_interim = SE_final = pval = rep(NA, N)

set.seed(7)

for(i in 1:nsim){
  # Stage 1:
  X0 = rbinom(N/4, 1, 0.2)
  X1 = rbinom(N/4, 1, 0.15)
  
  EE = sum(X1) # Events on Experimental arm at interim
  EN = N/4 - EE  # Non-events on Experimental arm at interim
  
  CE = sum(X0) # Events on Control arm at interim
  CN = N/4 - CE # Non-events on Control arm at interim
  
  RD_interim[i] = (EE/(EE+EN)) - (CE/(CE+CN))
  
  SE_interim[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)
  
  pval[i] = pnorm(RD_interim[i]/SE_interim[i])
  
  # Stage 2:
  Y0 = rbinom(N/4, 1, 0.2)
  Y1 = rbinom(N/4, 1, 0.15)
  
  EE = sum(X1+Y1) # Events on Experimental arm at end
  EN = N/2 - EE # Non-events on Experimental arm at end
  
  CE = sum(X0+Y0) # Events on Control arm at end
  CN = N/2 - CE # Non-events on Control arm at end
  
  RD_final[i] = (EE/(EE+EN)) - (CE/(CE+CN))
  
  SE_final[i] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)
  
}

stop_early = (pval <= 0.001)
continue = (pval > 0.001) 

# Plot

plot(RD_interim[continue], RD_final[continue], col = 'blue',
     xlim = c(-0.15, 0.05), ylim = c(-0.1, 0),
     xlab = 'Risk difference interim',
     ylab = 'Risk difference final')

points(RD_interim[stop_early], RD_final[stop_early], col = 'red')


# Table (slide 18)

rbind(
c(sum(stop_early), mean(RD_interim[stop_early]) + 0.05,
  mean(SE_interim[stop_early])),

c(sum(continue), mean(RD_final[continue]) + 0.05,
  mean(SE_final[continue])),

c(nsim, mean(c(RD_interim[stop_early], RD_final[continue])) + 0.05,
  mean(c(SE_interim[stop_early],SE_final[continue]))),

c(nsim, weighted.mean(c(RD_interim[stop_early], RD_final[continue]),
                      c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)) + 0.05,
  weighted.mean(c(SE_interim[stop_early], SE_final[continue]),
                c(1/SE_interim[stop_early]^2, 1/SE_final[continue]^2)))
)




########### 4 more studies with fixed sample size


RD_fixed = SE_fixed = matrix(nrow = nsim, ncol = 4)

for(j in 1:4){
  
  for(i in 1:nsim){
    
    X0 = rbinom(N/2, 1, 0.2)
    X1 = rbinom(N/2, 1, 0.15)
    
    EE = sum(X1)
    EN = N/2 - EE
    
    CE = sum(X0)
    CN = N/2 - CE
    
    RD_fixed[i,j] = (EE/(EE+EN)) - (CE/(CE+CN))
    
    SE_fixed[i,j] = sqrt(EE*EN/(EE+EN)^3 + CE*CN/(CE+CN)^3)
  }
}


# Meta-analysis

MA_all = MA_exclude = rep(NA, nsim)


for(i in 1:nsim){
  
  if(pval[i] <= 0.001){
    
    MA_all[i] = weighted.mean(c(RD_interim[i], RD_fixed[i,]),
                           c(1/SE_interim[i]^2, 1/SE_fixed[i,]^2))
    
    MA_exclude[i] = weighted.mean(c(RD_fixed[i,]), c(1/SE_fixed[i,]^2))
    
  } else {
    
    MA_all[i] = MA_exclude[i] = weighted.mean(c(RD_final[i], RD_fixed[i,]),
                           c(1/SE_final[i]^2, 1/SE_fixed[i,]^2))
    
  }
  
}


# Table (slide 19)

rbind(
c(nsim, 5, mean(MA_all)+0.05, sd(MA_all)),

c(nsim, 4+mean(continue), mean(MA_exclude)+0.05, sd(MA_exclude))
)
