# Check for the ineqx package
# Benjamin Rosche
# July 2021

rm(list=ls())

library(ineqx)

source("./create/createExampleData.R")

# ------------------------------------------------------------------------------------------------ #
# Violation of the PD principle of transfers
# ------------------------------------------------------------------------------------------------ #

# Situation 1 ------------------------------------------------------------------------------------ #

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

wibe(dat.PD1, y="inc", groupvar="group", timevar="year", long=F)[[1]] %>% dplyr::select(year, group, mu) %>% pivot_wider(names_from=year, names_prefix = "y", values_from=mu)
wibe(dat.PD1, y="lninc", groupvar="group", timevar="year", long=F)[[1]] %>% dplyr::select(year, group, mu) %>% pivot_wider(names_from=year, names_prefix = "y", values_from=mu)
# -> We can see that both variance and log variance record a change that makes the poor poorer and the rich richer

wibe(dat.PD1, y="inc", groupvar="group", timevar="year", long=F)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
wibe(dat.PD1, y="lninc", groupvar="group", timevar="year", long=F)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
# -> The log variance, however, records a decrease of inequality
# -> violation of the PD principle of transfers

# Situation 2 ------------------------------------------------------------------------------------ #

dat.PD2 <-
  tibble(year=1, group=1, inc=c(0.05,0.05,0.05)) %>%
  add_row(data.frame(year=1, group=2, inc=c(2,2.8284,11.3137))) %>%
  add_row(data.frame(year=2, group=1, inc=c(0.0499,0.05,0.0501))) %>%
  add_row(data.frame(year=2, group=2, inc=c(2.1,2.6284,11.4137))) %>%
  dplyr::mutate(inc=inc*1000, lninc=log(inc))

wibe(dat.PD2, y="inc", groupvar="group", timevar="year", long=F)[[1]] %>% dplyr::select(year, group, mu, sigma2)
wibe(dat.PD2, y="lninc", groupvar="group", timevar="year", long=F)[[1]] %>% dplyr::select(year, group, mu, sigma2)
# -> Here we see that the change does not affect the between-group distances; mu's are unchanged both with variance and log-variance
# -> Moreover, the within-group variance and the within-group log-variance is larger at time point 2

wibe(dat.PD2, y="inc", groupvar="group", timevar="year", long=F)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
wibe(dat.PD2, y="lninc", groupvar="group", timevar="year", long=F)[[2]] %>% dplyr::select(year, VarW, VarB, VarT)
# -> Accordingly, the total variance goes up
# -> The log-variance, however, goes down!
# -> This is because the between-group variance is affected by the within-group changes in the case of the log-variance. As we can see, log-VarB goes down!
# -> Therefore, while the within/between-group decomposition seemingly works, it does not properly differentiate the two (for the arithmetic mean)

# More info:
# For how Gini violates decomposability: Cowell 2011, p.66
# For how log-variance violates decomposability: Cowell 1988

# ------------------------------------------------------------------------------------------------ #
# Check dW + dB
# ------------------------------------------------------------------------------------------------ #

dat.dWB1 <-
  crDat(
    N.T=10,
    N.G=3,
    LP_n = c("500", "500", "500"),
    LP_x = c("0.5", "0.5", "0.5"),
    LP_mu = "10+100*(group==1)*x+150*(group==2)*x+20*(group==3)*x*year",
    LP_sigma = "1+10*(group==1)*x+10*(group==2)*x+10*(group==3)*x+10*(group==3)*x*year",
    sav = F
  )

ggplot(data=
         wibe(dat.dWB1 %>% dplyr::filter(x==1), y="inc", groupvar="group", timevar="year", long=T)[[2]] %>% filter(variable %in% c("CV2B", "CV2W", "CV2T")),
       aes(x=year, y=value, group=variable, color=factor(variable))) +
  geom_line() +
  labs(x="", color="")

ineqx.dWB1 <- ineqx(x="x", y="c.inc",   ystat="CV2", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dWB1)
ineqx.dWB2 <- ineqx(x="x", y="c.inc",   ystat="Var", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dWB1)
ineqx.dWB3 <- ineqx(x="x", y="c.lninc", ystat="Var", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dWB1)

plot(ineqx.dWB1, type="dMuP")
plot(ineqx.dWB1, type="dSigmaP")

plot(ineqx.dWB1, type="dPA")
plot(ineqx.dWB2, type="dPA")
plot(ineqx.dWB3, type="dPA")

ineqx.dWB3$dT[[1]]

# Comments:
# - ineqx() perfectly recovers changes due to dW and dB
# - We can see that variance, log-variance, and CV2 all predict a different trajectory of inequality

# ------------------------------------------------------------------------------------------------ #
# Check dD + dO
# ------------------------------------------------------------------------------------------------ #

# Situation 1 ------------------------------------------------------------------------------------ #

dat.dD1 <-
  crDat(
    N.T=3,
    N.G=3,
    LP_x = c("0.5", "0.5", "0.5"),
    LP_n     = c("1000", "1000", "1000-(year>2)*500"),
    LP_mu    = " 100*(group==1) +        100*(group==2) +       100*(group==3) +       1000*(group==2)*(year>1) +
                -100*(group==1)*(x==1) + 300*(group==2)*(x==1) +  100*(group==3)*(x==1)",
    LP_sigma = "10",
    seed=T,
    sav = F
  )

# -> After year 1, group 2 has $1000 more at x=0
# -> After year 2, the size of group 3 (n) shrinks from 1000 to 500

# CV2
ggplot(data=
         wibe(dat.dD1 %>% dplyr::filter(x==1), y="inc", groupvar="group", timevar="year", long=T)[[2]] %>% filter(variable %in% c("CV2B", "CV2W", "CV2T")),
       aes(x=year, y=value, group=variable, color=factor(variable))) +
  geom_line() + scale_x_continuous(breaks=seq(1,10,1)) + labs(x="", y="CV2", color="")
# -> inequality goes up after year 1 and down after year 2

# Variance
ggplot(data=
         wibe(dat.dD1 %>% dplyr::filter(x==1), y="inc", groupvar="group", timevar="year", long=T)[[2]] %>% filter(variable %in% c("VarB", "VarW", "VarT")),
       aes(x=year, y=value, group=variable, color=factor(variable))) +
  geom_line() + scale_x_continuous(breaks=seq(1,10,1)) + labs(x="", y="Var", color="")
# -> inequality goes up after year and up after year 2

# Explanation for YEAR 1 to YEAR 2 --------------------------------------------------------------- #

# CV2

ineqx.dD11 <- ineqx(x="i.x", y="c.inc", ystat="CV2", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dD1)

ineqx.dD11$dW[[1]]
# -> unaffected

ineqx.dD11$dB[[1]]
# -> dO is positive because group 2 has $1000 more at x=0
# -> dD is negative because the impact of X weighs less heavily given the overall amount those households have ~~ scale invariance feature
# -> dB is positive because the positive effect of X contributes to inequality

# Var

ineqx.dD12 <- ineqx(x="i.x", y="c.inc", ystat="Var", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dD1)

ineqx.dD12$dW[[1]]
# -> unaffected

ineqx.dD12$dB[[1]]
# -> dO is positive because group 2 has $1000 more at x=0
# -> dD is positive because the variance is not scale invariant. That is, because group 2 has more at x=0, the effect of x on the variance will be stronger
# -> dB is positive because the positive effect of X contributes to inequality

# Explanation for YEAR 2 to YEAR 3 --------------------------------------------------------------- #

wibe(dat.dD1 %>% dplyr::filter(x==1), y="inc", groupvar="group", timevar="year", long=F)[[1]] %>% dplyr::filter(year==3)
wibe(dat.dD1 %>% dplyr::filter(x==1), y="inc", groupvar="group", timevar="year", long=F)[[2]]

# -> CV2 is going down while Var is going up
# -> Var goes up because the group-size weighted distance, VarB, is going up
# -> CV2, however, is normalized by the grand mean. The grand mean increases since the mean of group 2 (y_bar=1400) weighs more heavily now. Accordingly, CV2 goes down.
# -> VarW decreases marginally but not enough for the Var to go down.

# Intepretation of the different changes --------------------------------------------------------- #

plot(ineqx.dD11, type="dPA")
plot(ineqx.dD12, type="dPA")

# dB: The effect of X on between-group component
# dW: The effect of X on within-group component
# dD: The effect of n and mu that goes through X on within + between-group component
# dO: The effect of n and mu that does not go through X on within + between-group component
# dT: dB+dW+dD == total effect of X in a broad sense
# dB+dW == total effect of X in a narrow sense ("coefficient effect")
# dT+dO == total effect in the sense, the change of inequality at x=0 plus the changes from x=0 to x=1
# Comment:
# -> I don't think we can call dD "characteristic effect" since it is not the quantitative change of X but of n and mu

# ------------------------------------------------------------------------------------------------ #
# Different distributions of X
# ------------------------------------------------------------------------------------------------ #

dat.WX1 <-
  crDat(
    N.T=4,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0.2', "0.9", "0.1*i"),
    LP_mu="10+100*(group==1)*x*t+150*(group==2)*x*t+20*(group==3)*x*t*year",
    LP_sigma="10+10*(group==1)*x*t*year",
    X_const = F,
    sav = F
  )

ineqx.WX1 <- ineqx(x="i.x", t="t", y="c.inc", ystat="CV2", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.WX1)
plot(ineqx.WX1, type="dPA")

# > Difference in difference estimator works with different x distributions

dat.WX2 <-
  crDat(
    N.T=4,
    N.G=3,
    LP_n = c("1000", "1000", "1000"),
    LP_x = c('0.5', "0.5", "0.2*i"),
    LP_mu="10+100*(group==1)*x+150*(group==2)*x+20*(group==3)*x*year",
    LP_sigma="10+10*(group==1)*x*year",
    X_const = F,
    sav = F
  )

ineqx.WX2 <- ineqx(x="i.x", y="c.inc", ystat="CV2", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.WX2)
plot(ineqx.WX2, type="dPA")

# Between-group estimator also works!

# ------------------------------------------------------------------------------------------------ #
# Check x_type == "between"
# ------------------------------------------------------------------------------------------------ #

dat.dWB2 <-
  crDat(
    N.T=5,
    N.G=3,
    LP_n     = c("1000", "1000", "1000"),
    LP_mu="10+100*(group==1)*x+150*(group==2)*x+20*(group==3)*x*year",
    LP_sigma="10+10*(group==1)*x*year",
    X_type = "between",
    sav = F
  )

ineqx.dWB2 <- ineqx(x="i.x", y="c.inc", ystat="Var", groupvar = "i.group", timevar = "i.year", cf=1, dat=dat.dWB2)

plot(ineqx.dWB2, type="dMuP")
plot(ineqx.dWB2, type="dSigmaP")

plot(ineqx.dWB2, type="dPA")

ineqx.dWB2$dT[[1]]

# Comments:
# - The problem is that inequality in the values of X creates another grouping.
#   For the between effect, we then need to differentiate between the effect change and the endowment change.
#   For the within level, however, it is more complicated because coefficient changes at the between
#   level will create within change in each group because of difference in x within groups
