
# rm(list=ls())
#
# library(ineqx)
# source("./create/createExampleData.R")
#
# dat1 <-
#   crDat(N.T=4,
#         N.G=3,
#         LP_n=c("1000", "1000", "1000"),
#         LP_x = c("0.1*i", "0.5", "0.5"),
#         LP_mu = "1000",
#         LP_sigma = "0.1+1000*(group==1)*x*t*year+5*(group==2)*x*t+0*(group==3)*x*t",
#         seed=T,
#         sav=F) %>%
#   dplyr::mutate(xt=x*t)
#
#
# ineqx(y="c.inc", ystat = "CV2", groupvar = "i.group", cf=1, dat = dat1)
#
# ineqx(x="x", y="c.inc", ystat = "CV2", groupvar = "i.group", cf=1, dat = dat1)
#
# ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", cf=1, dat = dat1)
#
# ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", timevar="i.year", cf=1, dat = dat1)
#
# ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", timevar = "i.year", cf=list(n=c(1000,1000,100),
#                                                                                                 mu=c(10,10,10),
#                                                                                                 sigma=c(0.1,0.1,0.1),
#                                                                                                 beta=c(0,0,0),
#                                                                                                 lambda=c(0,0,0)), dat = dat1)












