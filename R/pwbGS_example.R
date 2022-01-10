\dontrun{
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
  interims <- c(80, 100, 120, 140)
  pwbGS(theta=0.5, des=ad, interims=interims)
}
