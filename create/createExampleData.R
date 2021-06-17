
library(plyr)
library(dplyr)
library(tidyr)
library(haven)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")

crDat <- function(N.T, N.G, LP_n, LP_mu, LP_sigma, seed=T, sav=F) {

  if(seed==T) set.seed(1)

  # Dataframe with number of observations per wave

  N <-
    tibble() %>%
    tidyr::expand(year=1:N.T)

  for(i in 1:N.G) {
    N <-
      N %>%
      dplyr::mutate("G{i}" := eval(parse(text=LP_n[i])))

  }

  # Create data structure ------------------------------------------------------------------------ #

  d <- replicate(N.T, vector("list", N.G), simplify = FALSE)

  for(i in 1:N.T) { # per wave

    for(j in 1:N.G) { # per group

      d[[i]][[j]] <- tibble() %>% tidyr::expand(id=1:(N %>% dplyr::filter(year==i) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()), year=i, group=j, x=0:1)

    }

  }

  # Create dataframe from list of lists
  dat <- ldply(d, .fun = ldply)

  n <- dim(dat)[1]

  # Create income -------------------------------------------------------------------------------- #

  incdat <-
    dat %>%
    tibble() %>%
    dplyr::mutate(
      z = rpois(n, 2),
      inc = rnorm(n, eval(parse(text=LP_mu)), eval(parse(text=LP_sigma))),
      lninc = log(inc-(min(inc)-1)*(min(inc)<0)),
      w = 1/n
    )

  if(sav=="R") save(incdat, file = "data/incdat.RData") else if (sav=="Stata") write_dta(incdat, path = "data/incdat.dta")

  return(incdat)

}
