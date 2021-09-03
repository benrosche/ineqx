
library(plyr)
library(dplyr)
library(tidyr)
library(haven)

setwd("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx")

crDat <- function(N.T, N.G, LP_n, LP_x, LP_mu, LP_sigma, X_const=F, seed=T, sav=F) {

  # Arguments ------------------------------------------------------------------------------------ #
  # N.T      = Number of time points
  # N.G      = Number of groups
  # LP_n     = "Linear predictor" for n per group. This can be just numbers as characters
  #            LP_n = c("300", "300", "300") or functions LP_n = c("300-10*year", "300", "300")
  # LP_x     = Proportion of X=0 as "linear predictor". E.g., LP_x = c("0.1+0.1*i", "0.5", "0.9-0.1*i")
  # LP_mu    = Linear predictor for mean, e.g.
  # LP_sigma = Linear predictor for standard deviation
  # X_const  = {T|F} If True, distribution of X does not change across waves
  # X_ctrl   = c([0-1],[0-1],[0-1]); proportion of data that is control per group
  # seed     = {T|F}
  # sav      = {T|F}
  # ---------------------------------------------------------------------------------------------- #

  # N.T=5; N.G=3; LP_n = c("10", "10", "10"); LP_x <- c("0.5-0.1*i", "0.5", "0.5"); LP_mu = "10*(group==1)*x + 20*(group==2)*x + 30*(group==3)*x + z"; LP_sigma = "1"; X_const = T; seed=F; sav=F;

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

    for(i in 1:N.T) { # per time

      N.ij <- N %>% dplyr::filter(year==i) %>% dplyr::select(paste0("G", j)) %>% unlist() %>% as.vector()

      if(X_const & i > 1) d[[i]][[j]] <- d[[1]][[j]] %>% dplyr::mutate(year=i)
      else {
        d[[i]][[j]] <-
          tibble() %>%
          tidyr::expand(id=1:N.ij, year=i, group=j, x=1, t=c(0,5)) %>% # t=0:1
          dplyr::mutate(x=case_when(id %in% sample(1:N.ij, N.ij*eval(parse(text=LP_x[j]))) & x==1 ~ 0, TRUE ~ as.double(x))) #x==1
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

# incdat <-
#   crDat(5, 3, c("100", "100", "100"), c("0.5", "0.5", "0.5"), "10+10*(group==1)*x*t + 20*(group==2)*x*t + 10*(group==3)*x*t*year", "1", X_const=T, seed=T, sav=F)
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
