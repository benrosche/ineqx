# # ================================================================================================ #
# # Test of ineqx
# # ================================================================================================ #
#
# library(haven)
# library(dplyr)
#
# DIR <- "C:/Users/benja/OneDrive - Cornell University/GitHub/earningsinequality/"
#
# dat <- read_dta(paste0(DIR,"data/CPS/cps_final3.dta")) %>% dplyr::filter(N==2, bystart<=1, single==0, year2>2005)
# dat.ineqx <- dat %>% dplyr::mutate(by=byear*mother) %>% dplyr::select(year, f_SES, by, fearnings_wk, age, famsize, married, race, w_edu) %>% tidyr::drop_na()
#
# library(ineqx)
# library(tictoc)
#
# tic()
# i1 <- ineqx(dat.ineqx, x=by, y=fearnings_wk, groupvar=f_SES, timevar=year, controls = c("age", "famsize", "married", "race", "w_edu"))
# toc()
#
#
#
#
# library(ggthemes)
# library(ggpubr)
# library(ineqx)
#
# ggpubr2 <- theme_pubr() + theme(panel.grid.major.x = element_blank() , panel.grid.major.y = element_line(size=.1, color="grey", linetype = "dashed"),
#                                 legend.position="bottom", legend.margin=margin(-20,0,0,0), legend.box.margin=margin(0,0,0,0))
# theme_set(ggpubr2)
#
# # Descriptives ----------------------------------------------------------------------------------- #
#
# {
#
#   mu <-
#     ggplot(data=wibe(incdat %>% filter(by==1), y=inc, groupvar=SES, timevar=year)[[1]],
#            aes(x=year, y=mu, group=SES, color=factor(SES))) +
#     geom_line() +
#     ggpubr2 + scale_color_stata(labels=c("low", "medium", "high"))
#
#   sigma2 <-
#     ggplot(data=wibe(incdat %>% filter(by==1), y=inc, groupvar=SES, timevar=year)[[1]],
#            aes(x=year, y=sigma2, group=SES, color=factor(SES))) +
#     geom_line() +
#     ggpubr2 + scale_color_stata(labels=c("low", "medium", "high"))
#
#   n <-
#     ggplot(data=wibe(incdat, y=inc, groupvar=SES, timevar=year)[[1]] %>% group_by(year) %>% mutate(N=sum(n, na.rm = T), n=n/N),
#            aes(x=year, y=n*100, group=SES, color=factor(SES))) +
#     geom_line() +
#     ggpubr2 + scale_color_stata() + theme(legend.position="none") + labs(subtitle=expression(n), x="", y="%", color="")
#
#   CV2 <-
#     ggplot(data=wibe(incdat, y=inc, groupvar=SES, timevar=year)[[1]], aes(x=year, y=CV2, group=SES, color=factor(SES))) +
#     geom_line() +
#     ggpubr2 + scale_color_stata(labels=c("low", "medium", "high"))
#
#   CVWB <-
#     ggplot(data=
#              wibe(incdat %>% dplyr::filter(by==1), y=inc, groupvar=SES, timevar=year, long=T)[[2]] %>% filter(variable %in% c("CV2B", "CV2W", "CV2T")),
#            aes(x=year, y=value, group=variable, color=factor(variable))) +
#     geom_line() +
#     ggpubr2 + scale_color_stata(labels=c("Between", "Total", "Within"))
#
#   VarWB <-
#     ggplot(data=
#              wibe(incdat %>% filter(by==1), y=inc, groupvar=SES, timevar=year,  long=T)[[2]] %>% filter(variable %in% c("VarB", "VarW", "VarT")),
#            aes(x=year, y=value, group=variable, color=factor(variable))) +
#     geom_line() + geom_hline(yintercept=0.5) +
#     ggpubr2 + scale_color_stata(labels=c("Between", "Total", "Within"))
#
# }
#
# # ineqx ------------------------------------------------------------------------------------------ #
#
# {
#
#   ineqx.out1 <- ineqx(incdat, x=by, y=inc, ystat="CV2", groupvar = SES, timevar = year, cf=1)
#   ineqx.out2 <- ineqx(incdat, x=by, y=inc, ystat="Var", groupvar = SES, timevar = year, cf=1)
#
#   wibe(incdat %>% dplyr::filter(by==1), y=inc, groupvar=SES, timevar=year)[[2]]
#
#   ggplot(
#     ineqx.out1$dT[[1]] %>% mutate(dCV=CV2T-first(CV2T)),
#     aes(x=year)) +
#     geom_line(aes(y=dCV), color="black") +
#     geom_line(aes(y=dW), color="green") +
#     geom_line(aes(y=dB), color="blue") +
#     geom_line(aes(y=dD), color="orange") +
#     geom_line(aes(y=dO), color="grey") +
#     geom_line(aes(y=dD+dO), color="skyblue") +
#     geom_line(aes(y=dB+dW+dD+dO), color="red", linetype="dashed")
#
#   ggplot(
#     ineqx.out2$dT[[1]] %>% mutate(dVar=VarT-first(VarT)),
#     aes(x=year)) +
#     geom_line(aes(y=dVar), color="black") +
#     geom_line(aes(y=dW), color="green") +
#     geom_line(aes(y=dB), color="blue") +
#     geom_line(aes(y=dD), color="orange") +
#     geom_line(aes(y=dO), color="grey") +
#     geom_line(aes(y=dD+dO), color="skyblue") +
#     geom_line(aes(y=dB+dW+dD+dO), color="red", linetype="dashed")
#
#   plot(ineqx.out1, type = "dMuSigma")
#   plot(ineqx.out1, type = "dW")
#   plot(ineqx.out1, type = "dB")
#   plot(ineqx.out1, type = "dD")
#   plot(ineqx.out1, type = "dT")
#
# }
#
#
#
