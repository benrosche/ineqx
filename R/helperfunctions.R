
# ================================================================================================ #
# Function dW
# ================================================================================================ #

dW <- function(x, y, xpv, ystat, groupvar, timevar, cf, smoothDat=F, AME_mu, AME_sigma, dat) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1987, "w")
  # ---------------------------------------------------------------------------------------------- #

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})
  AME_mu <- AME_mu[[1]] %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }}) %>% dplyr::mutate(time=round(time)) # AME by group and time
  AME_sigma <- AME_sigma[[1]] %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }}) %>% dplyr::mutate(time=round(time))

  # Levels of group and time var
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as effectDat
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # wibe() --------------------------------------------------------------------------------------- #

  # @xpv[1]
  wibe.a <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[1]), smoothDat = smoothDat, rel = F)[[1]]

  # @xpv[2]
  wibe.b <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[2]), smoothDat = smoothDat, rel = F)[[2]]

  # Gather counterfactual data based on @xpv[1] -------------------------------------------------- #

  n <- wibe.a %>% filter(time==!!cf) %>% .$n

  mu <- wibe.a %>% filter(time==!!cf) %>% .$mu

  sigma <- wibe.a %>% filter(time==!!cf) %>% .$sigma

  beta <- (
    if(length(cf)==2) AME_mu %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else AME_mu %>% filter(time==!!cf) %>% .$beta %>% as.vector()
  )

  lambda <- (
    if(length(cf)==2) AME_sigma %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else AME_sigma %>% filter(time==!!cf) %>% .$lambda %>% as.vector()
  )

  # Create counterfactual dataset -------------------------------------------------------------- #

  dat.cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
    inner_join(tibble(group=group_levels, sigma.cf=sigma), by="group") %>%
    inner_join(tibble(group=group_levels, beta.cf=beta), by="group") %>%
    inner_join(tibble(group=group_levels, lambda.cf=lambda), by="group")

  # beta.cf: by gender
  if(length(cf)==2) {
    dat.cf <-
      dat.cf %>%
      dplyr::select(-beta.cf) %>%
      inner_join(
        AME_mu %>%
          dplyr::select(-!!paste0("beta_", cf[2]), -beta) %>% # remove beta and beta_x
          inner_join(tibble(group=group_levels, !!paste0("beta_", cf[2]):=beta), by="group") %>% # add cf beta_x
          dplyr::mutate(beta.cf=beta_w+beta_m) %>% # calculate cf beta
          dplyr::select(time, group, beta.cf),
        by=c("time", "group"))
  }

  # lambda.cf: by gender
  if(length(cf)==2) {
    dat.cf <-
      dat.cf %>%
      dplyr::select(-lambda.cf) %>%
      inner_join(
        AME_sigma %>%
          dplyr::select(-!!paste0("lambda_", cf[2]), -lambda) %>% # remove lambda and lambda_x
          inner_join(tibble(group=group_levels, !!paste0("lambda_", cf[2]):=lambda), by="group") %>% # add cf lambda_x
          dplyr::mutate(lambda.cf=lambda_w+lambda_m) %>% # calculate cf lambda
          dplyr::select(time, group, lambda.cf),
        by=c("time", "group"))
  }

  # Merge factual and counterfactual data
  dat.f_cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    left_join(AME_mu, by=c("time", "group")) %>%
    left_join(AME_sigma, by=c("time", "group")) %>%
    inner_join(wibe.a %>% dplyr::select(time, group, n, mu, sigma), by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f =  mu, sigma.f = sigma, beta.f = beta, lambda.f = lambda)

  # Calculate impact ----------------------------------------------------------------------------- #

  if(ystat=="CV2") {

    # Calculates the partial change in CV2W by group and time #
    impact.within <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~deltaCV2W is executed by year
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat="CV2W", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="CV2W", partial = T))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="CV2W", partial = T))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      ungroup() %>%
      dplyr::select(time, group, dW, dD, dO)

    # Calculates the total change in CV2W by time #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat="CV2W", partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="CV2W", partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="CV2W", partial = F))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dW, dD, dO)

  } else if (ystat=="Var") {

    # Calculates the partial change in VarW by SES #
    impact.within <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat="VarW", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="VarW", partial = T))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="VarW", partial = T))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      ungroup() %>%
      dplyr::select(time, group, dW, dD, dO)

    # Calculates the total change in VarW #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat="VarW", partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="VarW", partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="VarW", partial = F))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dW, dD, dO)

  }

  # Rename variables back to their original names ------------------------------------------------ #

  impact.within <-
    impact.within %>%
    inner_join(wibe.b %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
      !!enquo(groupvar) := group,
    )

  impact.within.total <-
    impact.within.total %>%
    inner_join(wibe.b %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
    )

  return(list(impact.within, impact.within.total))

}


# ================================================================================================ #
# Function dB
# ================================================================================================ #

dB <- function(x, y, xpv, ystat="CV2", groupvar, timevar, cf, smoothDat=F, AME_mu, dat) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1, "beta_w")
  # ---------------------------------------------------------------------------------------------- #

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})
  AME_mu <- AME_mu[[1]] %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }}) %>% dplyr::mutate(time=round(time)) # AME by group and time

  # Levels of group and time
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as AME_mu
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # wibe() --------------------------------------------------------------------------------------- #

  # @xpv[1]
  wibe.a <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[1]), smoothDat = smoothDat, rel = F)[[1]]

  # @xpv[2]
  wibe.b <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[2]), smoothDat = smoothDat, rel = F)[[2]]

  # Gather counterfactual data based on @xpv[1] -------------------------------------------------- #

  n  <- wibe.a %>% filter(time==!!cf) %>% .$n

  mu <- wibe.a %>% filter(time==!!cf) %>% .$mu

  beta <- (
    if(length(cf)==2) AME_mu %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else AME_mu %>% filter(time==!!cf) %>% .$beta %>% as.vector()
  )

  # Create counterfactual dataset ---------------------------------------------------------------- #

  dat.cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
    inner_join(tibble(group=group_levels, beta.cf=beta), by="group")

  # beta.cf: by gender
  if(length(cf)==2) {
    dat.cf <-
      dat.cf %>%
      dplyr::select(-beta.cf) %>%
      inner_join(
        AME_mu %>%
          dplyr::select(-!!paste0("beta_", cf[2]), -beta) %>% # remove beta and beta_x
          inner_join(tibble(group=group_levels, !!paste0("beta_", cf[2]):=beta), by="group") %>% # add cf beta_x
          dplyr::mutate(beta.cf=beta_w+beta_m) %>% # calculate cf beta
          dplyr::select(time, group, beta.cf),
        by=c("time", "group"))
  }

  # Merge factual and counterfactual data
  dat.f_cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    left_join(AME_mu, by=c("time", "group")) %>%
    inner_join(wibe.a %>% dplyr::select(time, group, n, mu), by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f =  mu, beta.f = beta)

  # Calculate impact ----------------------------------------------------------------------------- #

  # deltaCV2B calculates the effect of x (i.e.)

  if(ystat=="CV2") {

    # Calculates the partial change (i.e. computes the effect of one group while keeping the other group effects at 0) in CV2B by group and time #
    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~deltaCV2B is executed by year
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat="CV2B", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="CV2B", partial = T))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="CV2B", partial = T))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      ungroup() %>%
      dplyr::select(time, group, dB, dD, dO)

    # Calculates the total change (i.e. total impact of change in each group) in CV2B by time #
    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat="CV2B", partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="CV2B", partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="CV2B", partial = F))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dB, dD, dO)


  } else if (ystat=="Var") {

    # Calculates the partial change in VarB by group and time #
    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat="VarB", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="VarB", partial = T))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="VarB", partial = T))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::select(time, group, dB, dD, dO)

    # Calculates the total change in VarB by time #
    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat="VarB", partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdD(.x, ystat="VarB", partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdO(.x, ystat="VarB", partial = F))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::select(time, dB, dD, dO)

  }

  # Rename variables back to their original names ------------------------------------------------ #

  impact.partial <-
    impact.partial %>%
    inner_join(wibe.b %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
      !!enquo(groupvar) := group,
    )

  impact.total <-
    impact.total %>%
    inner_join(wibe.b %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
    )

  return(list(impact.partial, impact.total))

}

# ================================================================================================ #
# Function dD
# ================================================================================================ #


dD <- function(dW.out, dB.out) {

  ret <-
    purrr::map2(
      .x=dB.out, .y=dW.out,
      ~ .x %>%
        dplyr::select(dD, dO, matches("CV2T|Var2T")) %>%
        dplyr::mutate(dD= .x$dD+.y$dD,
                      dO= .x$dO+.y$dO,

        ) %>%
        ungroup()
    )

  return(list(tibble(dB.out[[1]][,1:2],ret[[1]]), tibble(dB.out[[2]][,1], ret[[2]])))

}

# ================================================================================================ #
# Function dT
# ================================================================================================ #

dT <- function(dW.out, dB.out, cf) {

  # Extract vars
  timevar  <- names(dW.out[[1]])[1]
  ystat <- names(dW.out[[1]])[7]
  ystat <- substr(ystat,1,nchar(ystat)-1) # remove T

  # Combine effects ------------------------------------------------------------------------------ #

  total <-
    dW.out[[2]] %>%
    dplyr::rename(dD.W=dD, dO.W=dO) %>%
    inner_join(
      dB.out[[2]] %>%
        dplyr::rename(dD.B=dD, dO.B=dO),
      by=c(timevar, paste0(ystat, "T"))) %>%
    dplyr::mutate(
      dD=dD.B+dD.W,
      dO=dO.B+dO.W,
      dT=dW+dB+dD,
      "dT+dO"=dT+dO
    ) %>%
    dplyr::select(-dD.W, -dD.B, -dO.W, -dO.B) %>%
    dplyr::relocate(paste0(ystat, c("W", "B", "T")), .after = last_col())


  # Create long dat
  total <- if(ystat == "CV2") total %>% dplyr::mutate(dCV2T=CV2T-CV2T[eval(parse(text = timevar)) == cf]) else  total %>% dplyr::mutate(dVarT=VarT-VarT[eval(parse(text = timevar)) == cf])

  # Calculate shares ----------------------------------------------------------------------------- #

  shares <-
    total %>%
    dplyr::select(!!timevar, dW, dB, dD, dO) %>%
    dplyr::mutate(abssum=abs(dW)+abs(dB)+abs(dD)+abs(dO)) %>%
    dplyr::mutate(across(c(dW, dB, dD, dO), ~abs(.x)/abssum)) %>%
    dplyr::select(-abssum) %>%
    pivot_longer(cols=-!!timevar, names_to="d", values_to="share")

  return(list(total, shares))

}


# ================================================================================================ #
# Function dYdX
# ================================================================================================ #

dYdX <- function(dat, ystat, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the group-specific effect of X on the within- or between-group component

    if(isTRUE(partial)) {

      # Partial or total effect?

      beta.star    <- beta.cf
      beta.star[i] <- beta.f[i]

      if(wb=="W") lambda.star    <- lambda.cf
      if(wb=="W") lambda.star[i] <- lambda.f[i]

    } else {

      beta.star <- beta.f
      if(wb=="W") lambda.star <- lambda.f

    }

    if(ystat=="CV2B") {

      dydx_t1 <- CV2B(n.cf, mu.cf + beta.star) - CV2B(n.cf, mu.cf)
      dydx_t0 <- CV2B(n.cf, mu.cf + beta.cf)   - CV2B(n.cf, mu.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="VarB") {

      dydx_t1 <- VarB(n.cf, mu.cf + beta.star) - VarB(n.cf, mu.cf)
      dydx_t0 <- VarB(n.cf, mu.cf + beta.cf)   - VarB(n.cf, mu.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="CV2W") {

      dydx_t1 <- CV2W(n.cf, mu.cf + beta.star, sigma.cf + lambda.star)                 - CV2W(n.cf, mu.cf, sigma.cf)
      dydx_t0 <- CV2W(n.cf, mu.cf + beta.cf, sigma.cf + lambda.cf) - CV2W(n.cf, mu.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="VarW") {

      dydx_t1 <- VarW(n.cf, sigma.cf + lambda.star)         - VarW(n.cf, sigma.cf)
      dydx_t0 <- VarW(n.cf, sigma.cf + lambda.cf) - VarW(n.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdD
# ================================================================================================ #

dYdD <- function(dat, ystat, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the group-specific effect of n and mu that is realized through X

    if(isTRUE(partial)) {

      # Partial or total effect?

      n.star     <- n.cf
      n.star[i]  <- n.f[i]

      mu.star    <- mu.cf
      mu.star[i] <- mu.f[i]

      if(wb=="W") sigma.star    <- sigma.cf
      if(wb=="W") sigma.star[i] <- sigma.f[i]

    } else {

      n.star  <- n.f
      mu.star <- mu.f
      if(wb=="W") sigma.star <- sigma.f

    }

    if(ystat=="CV2B") {

      dydx_t1 <- CV2B(n.star, mu.star + beta.f) - CV2B(n.star, mu.star)
      dydx_t0 <- CV2B(n.cf, mu.cf + beta.f)     - CV2B(n.cf, mu.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="VarB") {

      dydx_t1 <- VarB(n.star, mu.star + beta.f) - VarB(n.star, mu.star)
      dydx_t0 <- VarB(n.cf, mu.cf + beta.f)     - VarB(n.cf, mu.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="CV2W") {

      dydx_t1 <- CV2W(n.star, mu.star + beta.f, sigma.star + lambda.f) - CV2W(n.star, mu.star, sigma.star)
      dydx_t0 <- CV2W(n.cf, mu.cf + beta.f, sigma.cf + lambda.f)       - CV2W(n.cf, mu.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="VarW") {

      dydx_t1 <- VarW(n.star, sigma.star + lambda.f) - VarW(n.cf, sigma.star)
      dydx_t0 <- VarW(n.cf, sigma.cf + lambda.f)     - VarW(n.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdO
# ================================================================================================ #

dYdO <- function(dat, ystat, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the impact of n and mu that does not go through X

    if(partial==T) {

      # Partial or total effect?

      n.star     <- n.cf
      n.star[i]  <- n.f[i]

      mu.star    <- mu.cf
      mu.star[i] <- mu.f[i]

      if(wb=="W") sigma.star    <- sigma.cf
      if(wb=="W") sigma.star[i] <- sigma.f[i]

    } else {

      n.star  <- n.f
      mu.star <- mu.f
      if(wb=="W") sigma.star <- sigma.f

    }

    if(ystat=="VarW") {

      delta[i] <- VarW(n.star, sigma.star) - VarW(n.cf, sigma.cf)

    } else if(ystat=="VarB") {

      delta[i] <- VarB(n.star, mu.star) - VarB(n.cf, mu.cf)

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n.star, mu.star, sigma.star) - CV2W(n.cf, mu.cf, sigma.cf)

    } else if(ystat=="CV2B") {

      delta[i] <- CV2B(n.star, mu.star) - CV2B(n.cf, mu.cf)

    }

  }

  return(as.numeric(delta))

}

# ================================================================================================ #
# Function VarW
# ================================================================================================ #

VarW <- function(n, sigma) {

  w <- n/sum(n)

  return (sum(w * sigma^2))
}

# ================================================================================================ #
# Function VarB
# ================================================================================================ #

VarB <- function(n, mu) {

  w <- n/sum(n)
  gmu <- sum(w * mu)

  return( sum(w*(mu-gmu)^2) )
}

# ================================================================================================ #
# Function CV2B
# ================================================================================================ #

CV2B <- function(n, mu) {

  w <- n/sum(n)
  gmu <- sum(w * mu)

  return( sum(w * (mu-gmu)^2)/gmu^2 )
}

# ================================================================================================ #
# Function CV2W
# ================================================================================================ #

CV2W <- function(n, mu, sigma) {

  w <- n/sum(n)

  return( (sum(w * sigma^2)) / (sum(w * mu))^2 )
}
