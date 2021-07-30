
library(plyr)
library(dplyr)
library(tidyr)
library(haven)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")

crDat <- function(N.T, N.G, LP_n, LP_mu, LP_sigma, X_type, X_const=F, seed=T, sav=F) {

  # Arguments ------------------------------------------------------------------------------------ #
  # N.T      = Number of time points
  # N.G      = Number of groups
  # LP_n     = "Linear predictor" for n. This can be just numbers as characters
  #            LP_n = c("300", "300", "300") or functions LP_n = c("300-10*year", "300", "300")
  # LP_mu    = Linear predictor for mean, e.g.
  # LP_sigma = Linear predictor for standard deviation
  # X_type   = {"between", "within"}. If "between", causal effect is across ids. If "within",
  #            causal effect is within ids.
  # X_const  = {T|F} If True, distribution of X does not change across waves
  # seed     = {T|F}
  # sav      = {T|F}
  # ---------------------------------------------------------------------------------------------- #

  # N.T=5; N.G=3; LP_n = c("300", "300", "300"); LP_mu = "10*(group==1)*x + 20*(group==2)*x + 30*(group==3)*x + z"; LP_sigma = "1"; seed=F; sav=F; X_type="between"

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

  for(j in 1:N.G) { # per group

    # distribution of x constant across waves
    if(isTRUE(X_const)) {
      l <- N %>% dplyr::filter(year==1) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()
      x <- rpois(l, 2)
    }

    for(i in 1:N.T) { # per wave

      if(X_type=="between") {

        # distribution of x changes across waves
        if(isFALSE(X_const)) {
          l <- N %>% dplyr::filter(year==i) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()
          x <- rpois(l, 2)
        }

        d[[i]][[j]] <- tibble() %>% tidyr::expand(id=1:l, year=i, group=j) %>% dplyr::mutate(x=!!x, z=rpois(l, 0.1*x)) # treatment x, confounder z

      } else if(X_type=="within") {

        d[[i]][[j]] <- tibble() %>% tidyr::expand(id=1:(N %>% dplyr::filter(year==i) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()), year=i, group=j, x=0:1)

      }

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
      inc   = rnorm(n, eval(parse(text=LP_mu)), eval(parse(text=LP_sigma))),
      lninc = log(inc-(min(inc)-1)*(min(inc)<0))
    )

  if(sav=="R") save(incdat, file = "data/incdat.RData") else if (sav=="Stata") write_dta(incdat, path = "data/incdat.dta")

  return(incdat)

}
