# rm(list = ls())
#
# library(ineqx)
# library(haven)
# library(tidyverse)
# library(ggpubr)
# library(ggthemes)
# library(gridExtra)
# library(cowplot)
#
# setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")
#
# # ================================================================================================ #
# # Load data
# # ================================================================================================ #
#
# data <- read_dta("../parenthood/data/CPS/cps_final3.dta")
#
# dat.f <-
#   data %>%
#   dplyr::filter(N==2, bystart==0, single==0, year2>=1989) %>%
#   dplyr::select(id, year, year2, mother, bystart, byear, f_SES, fearnings_wk, lnfearnings_wk, earnwt, age, famsize, married, race, w_edu, single) %>%
#   drop_na() %>% # first drop na
#   group_by(id) %>% dplyr::mutate(n=byear+1, N=n()) %>% ungroup() %>% dplyr::filter(N>1) %>% ungroup() # recalculate n,N and drop incomplete cases
