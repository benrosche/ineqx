
library(plyr)
library(dplyr)
library(tidyr)
library(haven)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")

crDat <- function(N.T, N.G, LP_n, LP_x, LP_mu, LP_sigma, X_const=F, yp1=F, seed=T, sav=F) {

  # Arguments ------------------------------------------------------------------------------------ #
  # N.T      = Number of time points
  # N.G      = Number of groups
  # LP_n     = "Linear predictor" for n per group. This can be just numbers as characters
  #            LP_n = c("300", "300", "300") or functions LP_n = c("300-10*year", "300", "300")
  # LP_x     = Proportion of X=0 as "linear predictor". E.g., LP_x = c("0.1+0.1*i", "0.5", "0.9-0.1*i")
  # LP_mu    = Linear predictor for mean, e.g.
  # LP_sigma = Linear predictor for standard deviation
  # X_const  = {T|F} If True, distribution of X does not change across waves
  # yp1      = {T|F} IF True, {t=0, year}, {t=1, year+1}
  # seed     = {T|F}
  # sav      = {T|F}
  # ---------------------------------------------------------------------------------------------- #

  # N.T=3; N.G=3; LP_n = c("2", "2", "2"); LP_x <- c("0.5-0.1*i", "0.5", "0.5"); LP_mu = "10*(group==1)*x*t + 20*(group==2)*x*t + 30*(group==3)*x*t"; LP_sigma = "1"; X_const = F; seed=F; sav=F;

  if(seed==T) set.seed(1)
  if(isTRUE(sav)) stop("sav must be 'R', 'Stata', or FALSE.")

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

    for(i in 1:N.T) { # per time

      N.ij <- N %>% dplyr::filter(year==i) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()

      if(X_const & i > 1) d[[i]][[j]] <- d[[1]][[j]] %>% dplyr::mutate(year=i)
      else {
        d[[i]][[j]] <-
          tibble() %>%
          tidyr::expand(id=1:N.ij, year=i, group=j, x=1, t=0:1) %>% # t=0:1
          dplyr::mutate(x=case_when(id %in% sample(1:N.ij, N.ij*eval(parse(text=LP_x[j]))) & x==1 ~ 0, TRUE ~ as.double(x))) #x==1
      }

    }

  }

  # Create dataframe from list of lists
  dat <- ldply(d, .fun = ldply)

  # {t=0, year}, {t=1, year+1}?
  if(yp1) dat <- dat %>% dplyr::mutate(year=case_when(t>0 ~ as.double(year+1), TRUE ~ as.double(year)))

  n <- dim(dat)[1]

  # Create income -------------------------------------------------------------------------------- #

  incdat <-
    dat %>%
    tibble() %>%
    dplyr::mutate(
      inc   = rnorm(n, eval(parse(text=LP_mu)), eval(parse(text=LP_sigma))),
      lninc = log(inc-(min(inc)-1)*(min(inc)<0))
    )

  incdat <- incdat %>% dplyr::arrange(year, group, id, x, t) # sort

  if(sav=="R") save(incdat, file = "data/incdat.RData") else if (sav=="Stata") write_dta(incdat, path = "data/incdat.dta")

  return(incdat)

}


