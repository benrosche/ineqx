
# rm(list=ls())

library(ineqx)
library(stringr)

dat1 <-
  crDat(2, 3, c("1000", "1000", "1000"), "10+10*(group==1)*x+250*(group==2)*x+300*(group==3)*x", "0.1", "between", X_const = F, seed=T, sav=F)

dat2 <-
  crDat(2, 3, c("10", "10", "10"), "10+10*(group==1)*x+250*(group==2)*x+300*(group==3)*x*year", "0.1", "within", X_const = T, seed=T, sav=F)



wibe(y="inc", groupvar = "group", timevar = "year", dat=dat1)[[1]]

wibe(y="inc", groupvar = "group", timevar = "year", dat=dat1)[[2]] %>% dplyr::mutate(dCV2W=CV2W-first(CV2W),
                                                                                     dCV2B=CV2B-first(CV2B),
                                                                                     dCV2T=CV2T-first(CV2T),
                                                                                     dsum1=dCV2B+dCV2W,
                                                                                     dVarB=VarB-first(VarB),
                                                                                     dVarW=VarW-first(VarW),
                                                                                     dsum2=dVarB+dVarW,
                                                                                     dVarT=VarT-first(VarT))

# variance function regression will find 0.1 as SD but wibe shows that within variation is bigger because of x being continuous






