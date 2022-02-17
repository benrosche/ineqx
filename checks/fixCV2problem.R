rm(list=ls())

library(ineqx)

source("./create/createExampleData.R")

dat.WX2 <-
  crDat(
    N.T=2,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c("0.9", "0.5", "0.1*i"), #
    LP_mu="100*(group==1)+100*(group==2)+100*(group==3)*year",
    LP_sigma="0.1", #
    x_const = F,
    sav = F,
    seed = 1
  ) %>% group_by(id, year, group) %>% dplyr::filter(row_number()==1) %>% ungroup()

# Variance decomposition
ineqx.WX21 <- ineqx(treat="i.x", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.WX2)
plot(ineqx.WX21, type="dPA")
ineqx.WX21$dT[[1]]
# correct

# CV2 decomposition
ineqx.WX22 <- ineqx(treat="i.x", y="c.inc", ystat="CV2", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.WX2)
plot(ineqx.WX22, type="dPA")
ineqx.WX22$dT[[1]]
# incorrect



