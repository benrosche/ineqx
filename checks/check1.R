# Check for the ineqx package
# Benjamin Rosche
# July 2021 (last updated: Feb 2022)

rm(list=ls())

library(ineqx)

source("./create/createExampleData.R")

# ------------------------------------------------------------------------------------------------ #
# 1 Violation of the PD principle of transfers ----
# ------------------------------------------------------------------------------------------------ #

## Situation 1 =====================================================================================

dat.PD1 <-
  crDat(
    N.T=2,
    N.G=3,
    LP_n     = c("300", "600", "100"),
    LP_x     = c("0", "0", "0"),
    LP_mu    = "1000*(group==1) + 41000*(group==2) + 41000*(group==3) - 40000 * (group==2) * (year>1) + 240000 * (group==3) * (year>1)",
    LP_sigma = "1",
    sav = F
  )

library(gglorenz)
ggplot(aes(x=inc, group=as.factor(year), color=as.factor(year)), data = dat.PD1) + stat_lorenz(desc = T) + geom_abline(linetype = "dashed")

dat.PD1 %>% group_by(year) %>% dplyr::summarise(minc=mean(inc), mlninc=mean(lninc))
# -> The income change from year 1 to 2 is overall mean-preserving on the $ scale but mean changing on the log-$ scale

wibe(y="inc", group="group", time="year", long=F, dat=dat.PD1)[[1]] %>% dplyr::select(year, group, mu) %>% pivot_wider(names_from=year, names_prefix = "y", values_from=mu)
wibe(y="lninc", group="group", time="year", long=F, dat=dat.PD1)[[1]] %>% dplyr::select(year, group, mu) %>% pivot_wider(names_from=year, names_prefix = "y", values_from=mu)
# -> We can see that both variance and log variance record a change that makes the poor poorer and the rich richer

wibe(y="inc", group="group", time="year", long=F, dat=dat.PD1)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
wibe(y="lninc", group="group", time="year", long=F, dat=dat.PD1)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
# -> The log variance, however, records a decrease of inequality
# -> violation of the PD principle of transfers

## Situation 2 =====================================================================================

dat.PD2 <-
  tibble(year=1, group=1, inc=c(0.05,0.05,0.05)) %>%
  add_row(data.frame(year=1, group=2, inc=c(2,2.8284,11.3137))) %>%
  add_row(data.frame(year=2, group=1, inc=c(0.0499,0.05,0.0501))) %>%
  add_row(data.frame(year=2, group=2, inc=c(2.1,2.6284,11.4137))) %>%
  dplyr::mutate(inc=inc*1000, lninc=log(inc))

wibe(y="inc", group="group", time="year", long=F, dat=dat.PD2)[[1]] %>% dplyr::select(year, group, mu, sigma2)
wibe(y="lninc", group="group", time="year", long=F, dat=dat.PD2)[[1]] %>% dplyr::select(year, group, mu, sigma2)
# -> Here we see that the change does not affect the between-group distances; mu's are unchanged both with variance and log-variance
# -> Moreover, the within-group variance and the within-group log-variance is larger at time point 2

wibe(y="inc", group="group", time="year", long=F, dat=dat.PD2)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
wibe(y="lninc", group="group", time="year", long=F, dat=dat.PD2)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
# -> Accordingly, the total variance goes up
# -> The log-variance, however, goes down!
# -> This is because the between-group variance is affected by the within-group changes in the case of the log-variance. As we can see, log-VarB goes down!
# -> Therefore, while the within/between-group decomposition seemingly works, it does not properly differentiate the two (for the arithmetic mean)

# More info:
# For how Gini violates decomposability: Cowell 2011, p.66
# For how log-variance violates decomposability: Cowell 1988

# ------------------------------------------------------------------------------------------------ #
# 2 Check dW + dB ----
# ------------------------------------------------------------------------------------------------ #

dat.dWB1 <-
  crDat(
    N.T=10,
    N.G=3,
    LP_n = c("500", "500", "500"),
    LP_x = c("0.1", "0.1", "0.1"),
    LP_mu = "10+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma = "1+10*(group==1)*x*t+10*(group==2)*x*t+10*(group==3)*x*t+10*(group==3)*x*t*year",
    x_const = T
  )

ggplot(data=
         wibe(dat = dat.dWB1 %>% dplyr::filter(x==1), y="inc", group="group", time="year", long=T)[[2]] %>% filter(variable %in% c("CV2B", "CV2W", "CV2T")),
       aes(x=year, y=value, group=variable, color=factor(variable))) +
  geom_line() +
  labs(x="", color="")

ineqx.dWB1 <- ineqx(treat="x", post="t", y="inc", ystat="CV2", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.dWB1)
ineqx.dWB2 <- ineqx(treat="x", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.dWB1)
ineqx.dWB3 <- ineqx(treat="x", y="c.lninc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.dWB1)

plot(ineqx.dWB1, type="dMuP")
plot(ineqx.dWB1, type="dSigmaP")

plot(ineqx.dWB1, type="dPA")
ineqx.dWB1$dT[[1]]

plot(ineqx.dWB2, type="dPA")
ineqx.dWB2$dT[[1]]

plot(ineqx.dWB3, type="dPA")
ineqx.dWB3$dT[[1]]

# Comments:
# - ineqx() perfectly recovers changes due to dW and dB
# - We can see that variance, log-variance, and CV2 all predict a different trajectory of inequality
# - Both decomp=post and effect work

# ------------------------------------------------------------------------------------------------ #
# 3 Check dP ----
# ------------------------------------------------------------------------------------------------ #

# Next to the change of dB and dW that we keep for group 3 (same as in dat.dWB1),
# we now have group 2 increasing its pre-treatment mean after year 2.

dat.dD1 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c("0.2", "0.2", "0.2"),
    LP_mu = "10+100*(group==2)*(year>2)+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma = "1+10*(group==1)*x*t+10*(group==2)*x*t+10*(group==3)*x*t+10*(group==3)*x*t*year",
    x_const = T
  )

ineqx.dD11 <- ineqx(treat="i.x", post="i.t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.dD1)
plot(ineqx.dD11, type="dMuP")
plot(ineqx.dD11, type="dPA")
ineqx.dD11$dT[[1]]

# Variance decomposition:
# -> dP increases after year 2: 100*(group==2)*(year>2) - treatment effect of group 2 increases inequality even more than without pre-treatment increase of the mean of group 2
# -> dP decreases again after year 3 because the pre-treatment effect of group 2 is being offset by the increasing treatment effect of group 3
# -> dC unaffected

ineqx.dD12 <- ineqx(treat="i.x", post="i.t", y="c.inc", ystat="CV2", group = "i.group", time = "i.year", decom="post", ref=1, dat=dat.dD1)
plot(ineqx.dD12, type="dPA")
ineqx.dD12$dT[[1]]

# CV2 decomposition:
# -> Difference to variance is that the increasing treatment effect of group 3 increases the grand mean, which decreases inequality ~~ scale invariance

# ------------------------------------------------------------------------------------------------ #
# 4 Check dC ----
# ------------------------------------------------------------------------------------------------ #

# Next to the change of dB and dW that we keep for group 3 (same as in dat.dWB1),
# group 3 now also changes from 1000 to 200 individuals after year 1

dat.dC1 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000-(year>1)*800"),
    LP_x = c("0.2", "0.2", "0.2"),
    LP_mu = "10+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma = "1+10*(group==1)*x*t+10*(group==2)*x*t+10*(group==3)*x*t+10*(group==3)*x*t*year",
    x_const = F
  )

ineqx.dC11 <- ineqx(treat="i.x", post="i.t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.dC1)
plot(ineqx.dC11, type="dPA")
ineqx.dC11$dT[[1]]

# Variance decomposition
# -> dC indicates the reduction of observations in group 3 decreases inequality (as this group produced the most inequality)
#    However, this decrease is attenuated as the treatment effect of group 3 increases
# -> dP unaffected

ineqx.dC12 <- ineqx(treat="i.x", post="i.t", y="c.inc", ystat="CV2", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.dC1)
plot(ineqx.dC12, type="dPA")
ineqx.dC12$dT[[1]]

# CV2 decomposition:
# Same as variance

# ------------------------------------------------------------------------------------------------ #
# 5 Check dC + dP ----
# ------------------------------------------------------------------------------------------------ #

dat.dC2 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_x = c("0.1", "0.1", "0.1"),
    LP_n     = c("1000", "1000", "1000-(year>3)*800"),
    LP_mu = "10+100*(group==2)*(year>2)+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma = "1", # +10*(group==1)*x*t+10*(group==2)*x*t+10*(group==3)*x*t+10*(group==3)*x*t*year
    x_const = F
  )

ineqx.dC21 <- ineqx(treat="i.x", post="i.t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.dC2)
plot(ineqx.dC21, type="dPA")
ineqx.dC21$dT[[1]]

# dP: the change in pre-treatment earnings of group 2 after year 2 changes the effect of x
# dC after year 3, group 3 loses observations

# This changes the pre-treatment means (dP), the effect of x (dD)

# ------------------------------------------------------------------------------------------------ #
# 6 Check X distributions changes ----
# ------------------------------------------------------------------------------------------------ #

## x -----------------------------------------------------------------------------------------------

dat.WX2 <-
  crDat(
    N.T=2,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c("0.2", "0.9", "0.5"), # "0.2", "0.9", "0.1*i"
    LP_mu="10+100*(group==1)*x+150*(group==2)*x+20*(group==3)*x*year",
    LP_sigma="50", #+10*(group==1)*x*year
    x_const = F,
    sav = F,
    seed = 1
  ) %>% dplyr::filter(t==0)

# AME_mu <- list()
# AME_mu[[1]] <- data.frame(year=c(rep(1,3),rep(2,3)), group=rep(1:3,2), effect=c(100,150,20,100,150,40))
# AME_mu[[2]] <- data.frame(year=c(1,2), effect=c(90,96+2/3))
#
# AME_sigma <- list()
# AME_sigma[[1]] <- data.frame(year=c(rep(1,3),rep(2,3)), group=rep(1:3,2), effect=rep(0,6))
# AME_sigma[[2]] <- data.frame(year=c(1,2), effect=c(0,0))

# Variance decomposition
ineqx.WX21 <- ineqx(treat="i.x", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.WX2)
plot(ineqx.WX21, type="dPA")
ineqx.WX21$dT[[1]]

# Changing X distribution affects dC as changing x is reflected in a larger number of observations (=n) in group 3 who get the treatment

# CV2 decomposition
ineqx.WX22 <- ineqx(treat="i.x", y="c.inc", ystat="CV2", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.WX2)
plot(ineqx.WX22, type="dPA")
ineqx.WX22$dT[[1]]

# -> While dC is affected with CV2 as well, the majority of the change goes through dP
# -> Changing x distribution leads to higher post-treatment mean, which affects dP (effect of treatment via pre-treat)

## x and t -----------------------------------------------------------------------------------------

dat.WX3 <-
  crDat(
    N.T=9,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0.2', "0.9", "0.1*i"),
    LP_mu="10+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma="1+10*(group==1)*x*t*year",
    x_const = F,
    sav = F
  )

# Variance decomposition
ineqx.WX31 <- ineqx(treat="i.x", post="t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.WX3)
plot(ineqx.WX31, type="dPA")
ineqx.WX31$dT[[1]]

# CV2 decomposition
ineqx.WX32 <- ineqx(treat="i.x", post="t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="effect", ref=1, dat=dat.WX3)
plot(ineqx.WX32, type="dPA")
ineqx.WX32$dT[[1]]

# Works, especially for longer time series...

# ------------------------------------------------------------------------------------------------ #
# 7 ref are manual values ----
# ------------------------------------------------------------------------------------------------ #

dat.m <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("500", "500", "500"),
    LP_x = c("0.1", "0.1", "0.1"),
    LP_mu = "10+100*(group==1)*x*t+150*(group==2)*x*t+200*(group==3)*x*t*year",
    LP_sigma = "1",
    x_const = T
  )

ineqx.m1 <- ineqx(treat="x", post="t", y="inc", ystat="CV2", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.m)
plot(ineqx.m1, type = "dPA")
ineqx.m1$dT[[1]]



ineqx.m2 <- ineqx(treat="x", post="t", y="inc", ystat="CV2", group = "i.group", time = "i.year", decomp="post", ref=list(n=c(500,500,500), mu=c(10,10,10), sigma=c(1,1,1), beta=c(100,150,200), lambda=c(0,0,0)), dat=dat.m)
plot(ineqx.m2, type = "dPA")
ineqx.m2$dT[[1]]

## CONTINUE HERE ##
## dCV2T is NA WHY??? ##

# The two ways give the same results. Thus, the manual values work!

# ------------------------------------------------------------------------------------------------ #
# 8 Check averaging of treatment effects ----
# ------------------------------------------------------------------------------------------------ #

## x and x as quantitative variable ----------------------------------------------------------------

dat.AT1 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0', "0", "0"),
    LP_mu="10+100*(group==1)*t*year+150*(group==1)*t*year",
    LP_sigma="0.1",
    x_const = F,
    t_vals = c(0,1,2), # c(0,1,2,1000)
    sav = F) %>%
  dplyr::mutate(x=t) %>%
  dplyr::select(-t)

# Randomly remove one observation, so that each individual has just 2 observations, but some have
# 0->1 and others have 0->2
# dat.AT1 <-
#   dat.AT1 %>%
#   dplyr::mutate(z=runif(5*3*3*1000), z=case_when(t==0 ~ 0, TRUE ~ z)) %>%
#   group_by(year, group, id) %>%
#   dplyr::arrange(z, .by_group = T) %>%
#   dplyr::filter(row_number()<=2) %>%
#   ungroup() %>%
#   dplyr::select(-z)

# Variance decomposition
ineqx.AT1 <- ineqx(treat="i.x", post=NULL, y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.AT1)
plot(ineqx.AT1, type="dPA")
ineqx.AT1$dT[[1]]

# Averaging underestimates effect with x as quantitative variable, especially for high values of x
# No big changes with 2 instead of 3 observation per individual

## x and x as factorial variable -------------------------------------------------------------------

dat.AT2 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0', "0", "0"),
    LP_mu="10+100*(group==1)*(t==1)*year+150*(group==1)*(t==5)*year",
    LP_sigma="0.1",
    x_const = F,
    t_vals = c(0,1,5),
    sav = F) %>%
  dplyr::mutate(x=t) %>%
  dplyr::select(-t)


# Variance decomposition
ineqx.AT2 <- ineqx(treat="i.x", post=NULL, y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.AT2)
plot(ineqx.AT2, type="dPA")
ineqx.AT2$dT[[1]]

# dW within effect pops up
# The reason is that over years, the variance within group 1 grows

## x and t and t as quantitative variable ----------------------------------------------------------

dat.AT3 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0.1', "0.1", "0.1"),
    LP_mu="10+100*(group==1)*x*t*year+50*(group==1)*x*t*year-10*(group==3)*x*t*year-20*(group==3)*x*t*year",
    LP_sigma="1",
    x_const = F,
    t_vals = c(0,1,5),
    sav = F)

# Variance decomposition
ineqx.AT3 <- ineqx(treat="i.x", post="t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.AT3)
plot(ineqx.AT3, type="dPA")
ineqx.AT3$dT[[1]]

# Again, just underestimates the effect

## x and t and t as factorial variable -------------------------------------------------------------

dat.AT4 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0.1', "0.1", "0.1"),
    LP_mu="10+100*(group==1)*x*(t==1)*year+50*(group==1)*x*(t==5)*year-10*(group==3)*x*(t==1)*year-20*(group==3)*x*(t==5)*year",
    LP_sigma="1",
    x_const = F,
    t_vals = c(0,1,5),
    sav = F)

# dat.AT4 <-
#   dat.AT4 %>%
#   dplyr::mutate(z=runif(5*3*3*1000), z=case_when(t==0 ~ 0, TRUE ~ z)) %>%
#   group_by(year, group, id) %>%
#   dplyr::arrange(z, .by_group = T) %>%
#   dplyr::filter(row_number()<=2) %>%
#   ungroup() %>%
#   dplyr::select(-z)

# Variance decomposition
ineqx.AT4 <- ineqx(treat="i.x", post="t", y="c.inc", ystat="Var", group = "i.group", time = "i.year", decomp="post", ref=1, dat=dat.AT4)
plot(ineqx.AT4, type="dPA")
ineqx.AT4$dT[[1]]

# Again, when using t as factorial variable, a within effect pops up


# eof
