
rm(list=ls())

library(ineqx)
source("./create/createExampleData.R")

dat1 <-
  crDat(N.T=4,
        N.G=3,
        LP_n=c("1000", "1000", "1000"),
        LP_x = c("0.1", "0.1", "0.1"),
        LP_mu = "1000",
        LP_sigma = "0.1+1000*(group==1)*x*t*year+5*(group==2)*x*t+0*(group==3)*x*t",
        seed=T,
        sav=F) %>%
  dplyr::mutate(xt=x*t)



f1 <- ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", timevar="i.year", ref=2, dat = dat1)



ineqx(y="c.inc", ystat = "CV2", groupvar = "i.group", ref=1, dat = dat1)

ineqx(x="x", y="c.inc", ystat = "CV2", groupvar = "i.group", ref=1, dat = dat1)

ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", ref=1, dat = dat1)


ineqx(x="x", t="t", y="c.inc", ystat = "CV2", groupvar = "i.group", timevar = "i.year", ref=list(n=c(1000,1000,100),
                                                                                                 mu=c(10,10,10),
                                                                                                 sigma=c(0.1,0.1,0.1),
                                                                                                 beta=c(0,0,0),
                                                                                                 lambda=c(0,0,0)), dat = dat1)


plot(f1, type="dPA")




# crDat(5, 3, c("100", "100", "100"), c("0.1*i", "0.5", "0.5"), "10+10*(group==1)*x*t + 20*(group==2)*x*t + 10*(group==3)*x*t*year", "1", X_const=T, seed=T, sav="R")
#
# incdat %>%
#   dplyr::filter(x>0,t>0) %>%
#   group_by(year) %>%
#   dplyr::summarise(CV2T=ineq::var.coeff(inc, square = T)) %>%
#   dplyr::mutate(d=CV2T-first(CV2T))
#
# incdat %>% dplyr::select(-x) %>% group_by(year, group, t) %>% dplyr::summarise(n=n())
#
# incdat %>%
#   dplyr::select(-x, -lninc) %>%
#   bind_rows(incdat %>% dplyr::select(-x) %>% dplyr::filter(t>0) %>% group_by(year, group, id) %>% dplyr::summarise(inc=mean(inc), t=1)) %>%
#   arrange(year, group, id, t) %>%
#   dplyr::filter(t<=1) %>%
#   dplyr::filter(t>0) %>%
#   group_by(year) %>%
#   dplyr::summarise(CV2T=ineq::var.coeff(inc, square = T)) %>%
#   dplyr::mutate(d=CV2T-first(CV2T))
#
# set.seed(1)
# incdat2 <-
#   incdat %>%
#   dplyr::mutate(z=runif(5*3*3*100), z=case_when(t==0 ~ 0, TRUE ~ z)) %>%
#   group_by(year, group, id) %>%
#   dplyr::arrange(z, .by_group = T) %>%
#   dplyr::filter(row_number()<=2) %>%
#   ungroup() %>%
#   dplyr::select(-z)
#
# incdat2 %>% group_by(year, group, x, t) %>% dplyr::summarise(n=n())
#
# incdat2 %>%
#   dplyr::filter(x>0,t>0) %>%
#   group_by(year) %>%
#   dplyr::summarise(CV2T=ineq::var.coeff(inc, square = T)) %>%
#   dplyr::mutate(d=CV2T-first(CV2T))
#
# # How does this change?
#
# library(ineqx)
#
# t1 <- ineqx(x="i.x", t="i.t", y="inc", ystat="CV2", groupvar = "i.group", timevar = "i.year", ref=1, dat=incdat)
# t2 <- ineqx(x="i.x", t=NULL, y="inc", ystat="CV2", groupvar = "i.group", timevar = "i.year", ref=1, dat=incdat)
#
# t1$dT[[1]][,c("dT", "CV2T", "dCV2T")]
#
# plot(t1, type = "dPA")
# plot(t2, type = "dPA")





