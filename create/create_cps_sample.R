# To create a CPS sample with all young mothers but only 1000 nonmothers each year

library(haven)
library(dplyr)
library(ineqx)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/earningsinequality")

data <-
  read_dta("./data/cps_sheela.dta") %>%
  dplyr::filter(N==2, bystart<=0) %>%
  dplyr::mutate(year=year-1+minyear) %>%
  dplyr::mutate(across(everything(), as.vector))

nonmothers <-
  data %>%
  dplyr::select(id, year, n, mother) %>%
  dplyr::filter(mother==0, n==1) %>%
  group_by(year) %>%
  slice_sample(n=500) %>%
  ungroup() %>%
  dplyr::mutate(sample=1) %>%
  dplyr::select(id, sample)

cps_sample <-
  data %>%
  left_join(nonmothers, by="id") %>%
  dplyr::mutate(
    sample=case_when(
      is.na(sample) & mother == 1 ~ 1,
      is.na(sample) & mother == 0 ~ 0,
      TRUE ~ sample)) %>%
  dplyr::filter(sample==1) %>%
  dplyr::mutate(
    bym_age = byear*mother*age
  )

rm(data, nonmothers)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")

save(cps_sample, file="./data/cps_sample.RData")


