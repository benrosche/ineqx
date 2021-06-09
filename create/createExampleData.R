library(dplyr)
library(tidyr)
library(plyr)

crDat <- function(N.obs, T.obs, LP_mu, LP_sigma, sav=F) {

  N <- tibble() %>% expand(year=1:T.obs) %>% mutate(SES1=N.obs/3, SES2=N.obs/3, SES3=N.obs/3)

  # Create data structure ------------------------------------------------------------------------ #

  d <- list()
  for(i in 1:T.obs) {
    d[[i]] <-
      tibble() %>% expand(id=1:(N %>% dplyr::filter(year==i) %>% .$SES1), child=0:1, SES=1, year=i) %>%
      add_row(tibble() %>% expand(id=1:(N %>% dplyr::filter(year==i) %>% .$SES2), child=0:1, SES=2, year=i)) %>%
      add_row(tibble() %>% expand(id=1:(N %>% dplyr::filter(year==i) %>% .$SES3), child=0:1, SES=3, year=i))
  }

  dat <- ldply(d)

  n <- dim(dat)[1]

  # Create income -------------------------------------------------------------------------------- #

  incdat <- dat %>% mutate(inc=rnorm(n, eval(parse(text=LP_mu)), eval(parse(text=LP_sigma)))) %>% tibble()

  if(sav==T) save(incdat, file = "data/incdat.RData")

  return(incdat)

}

# Good example dat for dB:
# incdat <-
#   crDat(
#     N.obs=1000,
#     T.obs=5,
#     LP_mu="10*(SES==1)+20*(SES==2)+30*(SES==3)+100*(SES==1)*child+150*(SES==2)*child+200*(SES==3)*child+10*(SES==3)*child*year",
#     LP_sigma="10+20*(SES==1)*child+20*(SES==2)*child+30*(SES==3)*child+2*(SES==3)*child*year",
#     sav = T
#   )


# Good example dat for dW:
# incdat <-
#   crDat(
#     N.obs=1000,
#     T.obs=5,
#     LP_mu="10*(SES==1)+20*(SES==2)+30*(SES==3)+100*(SES==1)*child+150*(SES==2)*child+200*(SES==3)*child",
#     LP_sigma="10+20*(SES==1)*child+20*(SES==2)*child+30*(SES==3)*child+20*(SES==3)*child*year",
#     sav = T
#   )
