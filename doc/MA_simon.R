## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(pwb)

## ---- simon--------------------------------------------------------------
all.designs <- clinfun::ph2simon(pu=0.5,
                                 pa=0.6,
                                 ep1=0.05,
                                 ep2=0.1,
                                 nmax=500)
# Choose optimal Simon design:
opt.index <- which.min(all.designs$out[, "EN(p0)"])
simon.des <- all.designs$out[opt.index, ]

## ------------------------------------------------------------------------
ma.simon.results <- maSimon(theta=0.5,
                            des=simon.des,
                            nsims=1e4,
                            n.studies=4)
ma.simon.results

