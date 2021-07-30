# ================================================================================================ #
# Function dissectVar
# ================================================================================================ #

dissectVar <- function(x, cicheck="ci") {

  if(cicheck=="ci" & !grepl("(c|i).?\\.", x)) stop(paste0(x," must be specified as either c.", x, " (continuous) or i.", x, " (categorical)."))
  if(cicheck=="i"  & grepl("(c).?\\.", x)) stop(paste0(x," must be i.", as.symbol(sub("(c|i).?\\.", "", x)), " (categorical)"))
  if(cicheck=="c"  & grepl("(i).?\\.", x)) stop(paste0(x," must be c.", as.symbol(sub("(c|i).?\\.", "", x)), " (continuous)"))

  xpoly <- suppressWarnings(as.numeric(substr(x, 2,2)))
  if(is.na(xpoly)) xpoly <- 1

  xcont <- grepl("c.?\\.", x)

  x <- as.symbol(sub("(c|i).?\\.", "", x))

  return(list("var"=x, "xcont"=xcont, "xpoly"=xpoly))

}

# ================================================================================================ #
# Function calcAME
# ================================================================================================ #

calcAME <- function(x=NULL, xpv, groupvar, timevar, what, vfr1, vfr2, dat) {

  if(!is.null(x)) {

    nT <- length(unique(dat$time))

    # Two counterfactuals
    dat.CX <- dat %>% dplyr::mutate(x=rep(dat %>% dplyr::filter(time==!!cf) %>% .$x, nT)) # x constant across time

    # Predictions from variance function regression
    if(what=="mu") {

      pred <-
        dat %>%
        dplyr::select(time, group, x) %>%
        dplyr::mutate(mu   =predict(vfr1, what = "mu", type = "response", data = dat)) %>%
        dplyr::mutate(mu.CX=predict(vfr1, what = "mu", type = "response", newdata = dat.CX, data = dat)) %>%
        dplyr::mutate(mu.CB=predict(vfr2, what = "mu", type = "response", newdata = dat, data = dat %>% dplyr::filter(time==!!cf)))

    } else if(what=="sigma") {

      pred <-
        dat %>%
        dplyr::select(time, group, x) %>%
        dplyr::mutate(sigma   =predict(vfr1, what = "sigma", type = "response", data = dat)) %>%
        dplyr::mutate(sigma.CX=predict(vfr1, what = "sigma", type = "response", newdata = dat.CX, data = dat)) %>%
        dplyr::mutate(sigma.CL=predict(vfr2, what = "sigma", type = "response", newdata = dat, data = dat %>% dplyr::filter(time==!!cf)))

    }

    # Calculate Average Marginal Effects

    AME <- list()

    # Summarize by group and time
    AME[[1]] <-
      pred %>%
      group_by(time, group) %>%
      dplyr::summarise(across(everything(), ~ mean(.))) %>%
      ungroup() %>%
      dplyr::select(-x) %>%
      dplyr::rename(
        !!enquo(groupvar) := group,
        !!enquo(timevar) := time
      )

    # Summarize by time
    AME[[2]] <-
      pred %>%
      group_by(time) %>%
      dplyr::summarise(across(everything(), ~ mean(.))) %>%
      ungroup() %>%
      dplyr::select(-x) %>%
      dplyr::rename(
        !!enquo(timevar) := time
      )

  } else {

    AME <- list()

    AME[[1]] <-
      dat %>%
      dplyr::select(time, group) %>%
      group_by(time, group) %>%
      dplyr::summarise(beta=0) %>%
      dplyr::rename(
        !!enquo(groupvar) := group,
        !!enquo(timevar) := time
      )

    AME[[2]] <-
      dat %>%
      dplyr::select(time) %>%
      group_by(time) %>%
      dplyr::summarise(beta=0) %>%
      dplyr::rename(
        !!enquo(timevar) := time
      )

  }

  return(AME)

}


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

  if(!is.null(x)) {
    # wibe.a <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[1]), smoothDat = smoothDat, rel = F)[[1]] # @xpv[1]
    # wibe.b <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[2]), smoothDat = smoothDat, rel = F)[[2]] # @xpv[2]
    wibe.a <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[1]]
    wibe.b <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[2]]
  } else {
    wibe.a <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[1]]
    wibe.b <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[2]]
  }

  # Gather counterfactual data based on @xpv[1] -------------------------------------------------- #

  n     <- wibe.a %>% dplyr::filter(time==!!cf) %>% .$n

  mu    <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu

  mu.CX <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu.CX

  mu.CB <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu.CB

  sigma    <- AME_sigma %>% dplyr::filter(time==!!cf) %>% .$sigma

  sigma.CX <- AME_sigma %>% dplyr::filter(time==!!cf) %>% .$sigma.CX

  sigma.CL <- AME_sigma %>% dplyr::filter(time==!!cf) %>% .$sigma.CS

  # Create counterfactual dataset -------------------------------------------------------------- #

  dat.cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
    inner_join(tibble(group=group_levels, mu.CX.cf=mu.CX), by="group") %>%
    inner_join(tibble(group=group_levels, mu.CB.cf=mu.CB), by="group") %>%
    inner_join(tibble(group=group_levels, sigma.cf=sigma), by="group") %>%
    inner_join(tibble(group=group_levels, sigma.CX.cf=sigma.CX), by="group") %>%
    inner_join(tibble(group=group_levels, sigma.CL.cf=sigma.CL), by="group")

  # Merge factual and counterfactual data
  dat.f_cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(wibe.a %>% dplyr::select(time, group, n), by=c("time", "group")) %>%
    left_join(AME_mu, by=c("time", "group")) %>%
    left_join(AME_sigma, by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f = mu, mu.CX.f = mu.CX, mu.CB.f = mu.CB, sigma.f = sigma, sigma.CX.f = sigma.CX, sigma.CL.f = sigma.CL)


  # Calculate impact ----------------------------------------------------------------------------- #

  if(!is.null(x)) {

    # dWdX
    # (partial change = effect of one group while keeping the other group effects at 0)

    impact.within <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~ fcts are executed by year
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="CX", partial = T))) %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="CL", partial = T))) %>%
      dplyr::mutate(dI = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="I", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x, ystat=paste0(ystat, "W"), partial = T))) %>%
      unnest(cols = c(data, dW, dX, dI, dD)) %>%
      ungroup() %>%
      dplyr::select(time, group, dW, dX, dI, dD)

    # Calculates the total change in CV2W by time #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="CX", partial = F))) %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="CL", partial = F))) %>%
      dplyr::mutate(dI = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "W"), what="I", partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x, ystat=paste0(ystat, "W"), partial = F))) %>%
      unnest(cols = c(data, dW, dX, dI, dD)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dW, dX, dI, dD)

  } else {

    # Calculates the partial change in CV2W by group and time #
    impact.within <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=paste0(ystat, "W"), partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=paste0(ystat, "W"), partial = T))) %>%
      unnest(cols = c(data, dW, dD)) %>%
      ungroup() %>%
      dplyr::select(time, group, dW, dD)

    # Calculates the total change in CV2W by time #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=paste0(ystat, "W"), partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=paste0(ystat, "W"), partial = F))) %>%
      unnest(cols = c(data, dW, dD)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dW, dD)

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

  # smoothDat=F

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1, "beta_w")
  # ---------------------------------------------------------------------------------------------- #

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})

  # Levels of group and time
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as AME_mu
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))
  AME_mu <- AME_mu[[1]] %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }}) %>% dplyr::mutate(time=round(time)) # AME by group and time

  # wibe() --------------------------------------------------------------------------------------- #

  if(!is.null(x)) {
    # wibe.a <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[1]), smoothDat = smoothDat, rel = F)[[1]] # @xpv[1]
    # wibe.b <- wibe(y, group, time, dat %>% dplyr::filter(x == xpv[2]), smoothDat = smoothDat, rel = F)[[2]] # @xpv[2]
    wibe.a <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[1]]
    wibe.b <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[2]]
  } else {
    wibe.a <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[1]]
    wibe.b <- wibe(y, group, time, dat, smoothDat = smoothDat, rel = F)[[2]]
  }

  # Gather counterfactual data based on @xpv[1] -------------------------------------------------- #

  n     <- wibe.a %>% dplyr::filter(time==!!cf) %>% .$n

  mu    <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu

  mu.CX <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu.CX

  mu.CB <- AME_mu %>% dplyr::filter(time==!!cf) %>% .$mu.CB

  # Create counterfactual dataset ---------------------------------------------------------------- #

  dat.cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
    inner_join(tibble(group=group_levels, mu.CX.cf=mu.CX), by="group") %>%
    inner_join(tibble(group=group_levels, mu.CB.cf=mu.CB), by="group")

  # Merge factual and counterfactual data
  dat.f_cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(wibe.a %>% dplyr::select(time, group, n), by=c("time", "group")) %>%
    left_join(AME_mu, by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f = mu, mu.CX.f =  mu.CX, mu.CB.f = mu.CB)

  # Calculate impact ----------------------------------------------------------------------------- #

  if(!is.null(x)) {

    # Effect = dYdX
    # (partial change = effect of one group while keeping the other group effects at 0)

    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~ fcts are executed by year
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "CX", partial = T))) %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "CB", partial = T))) %>%
      dplyr::mutate(dI = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "I", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x, ystat=paste0(ystat, "B"), partial = T))) %>%
      unnest(cols = c(data, dB, dX, dI, dD)) %>%
      ungroup() %>%
      dplyr::select(time, group, dB, dX, dI, dD)

    # Calculates the total change (i.e. total impact of change in each group) in CV2B by time #
    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "CX", partial = T))) %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "CB", partial = T))) %>%
      dplyr::mutate(dI = purrr::map(.x = data, ~ dYdX(.x, ystat=paste0(ystat, "B"), what = "I", partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x, ystat=paste0(ystat, "B"), partial = T))) %>%
      unnest(cols = c(data, dB, dX, dI, dD)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dB, dX, dI, dD)

  } else {

    # Effect = dYdMuSigmaN

    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=paste0(ystat, "B"), partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=paste0(ystat, "B"), partial = T))) %>%
      unnest(cols = c(data, dB, dD)) %>%
      ungroup() %>%
      dplyr::select(time, group, dB, dD)

    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=paste0(ystat, "B"), partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=paste0(ystat, "B"), partial = F))) %>%
      unnest(cols = c(data, dB, dD)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dB, dD)

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

dD <- function(dW.out, dB.out, dO=T) {

  if(isTRUE(dO)) {

    ret <-
      purrr::map2(
        .x=dB.out, .y=dW.out,
        ~ .x %>%
          dplyr::select(dD, dO, matches("CV2T|Var2T")) %>%
          dplyr::mutate(
            dD=.x$dD+.y$dD,
            dO=.x$dO+.y$dO
          ) %>%
          ungroup()
      )

  } else {

    ret <-
      purrr::map2(
        .x=dB.out, .y=dW.out,
        ~ .x %>%
          dplyr::select(dD, matches("CV2T|Var2T")) %>%
          dplyr::mutate(dD=.x$dD+.y$dD) %>%
          ungroup()
      )

  }

  return(list(tibble(dB.out[[1]][,1:2],ret[[1]]), tibble(dB.out[[2]][,1], ret[[2]])))

}


# ================================================================================================ #
# Function dT
# ================================================================================================ #

dT <- function(dW.out, dB.out, x, cf, ystat) {

  # Extract vars
  timevar  <- names(dW.out[[1]])[1]

  # Combine effects ------------------------------------------------------------------------------ #

  if(!is.null(x)) {

    total <-
      dW.out[[2]] %>%
      dplyr::rename(dD.W=dD, dO.W=dO) %>%
      inner_join(
        dB.out[[2]] %>%
          dplyr::rename(dD.B=dD, dO.B=dO) %>%
          dplyr::select(-paste0(ystat, "T")),
        by=c(timevar)) %>%
      dplyr::mutate(
        dD=dD.W+dD.B,
        dT=dW+dB+dD,
        dO=dO.W+dO.B
      ) %>%
      dplyr::select(!!timevar, dW, dB, dD, dT, dO, paste0(ystat, c("W", "B", "T")))

  } else {

    total <-
      dW.out[[2]] %>%
      dplyr::rename(dD.W=dD) %>%
      inner_join(
        dB.out[[2]] %>%
          dplyr::rename(dD.B=dD) %>%
          dplyr::select(-paste0(ystat, "T")),
        by=c(timevar)) %>%
      dplyr::mutate(
        dD=dD.W+dD.B,
        dT=dW+dB+dD
      ) %>%
      dplyr::select(!!timevar, dW, dB, dD, dT, paste0(ystat, c("W", "B", "T")))

  }

  # Create long dat
  total <- if(ystat == "CV2") total %>% dplyr::mutate(dCV2T=CV2T-CV2T[eval(parse(text = timevar)) == cf]) else  total %>% dplyr::mutate(dVarT=VarT-VarT[eval(parse(text = timevar)) == cf])

  # Calculate shares ----------------------------------------------------------------------------- #

  shares <-
    total %>%
    dplyr::select(-matches("Var"), -matches("CV2"), -"dT") %>%
    dplyr::mutate(across(everything(), ~abs(.x))) %>%
    dplyr::mutate(abssum=rowSums(across(c(everything(), -!!timevar)))) %>%
    dplyr::mutate(across(c(everything(), -!!timevar), ~.x/abssum)) %>%
    dplyr::select(-abssum) %>%
    pivot_longer(cols=-!!timevar, names_to="d", values_to="share")

  return(list(total, shares))

}


# ================================================================================================ #
# Function dYdX
# This function manipulates beta and lambda
# ================================================================================================ #

dYdX <- function(dat, ystat, what, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the group-specific effect of X on the within- or between-group component

    if(isTRUE(partial)) {

      # Partial or total effect?

      beta.star    <- 0 # beta.cf
      beta.star[i] <- 0 # beta.f[i]

      if(wb=="W") lambda.star    <- lambda.cf
      if(wb=="W") lambda.star[i] <- lambda.f[i]

    } else {

      beta.star <- 0 #beta.f
      if(wb=="W") lambda.star <- lambda.f

    }

    if(ystat=="CV2B") {


      if(what=="CX") {

        delta[i] <- CV2B(n.f, mu.CX.f)  - CV2B(n.cf, mu.CX.cf)

      } else if(what=="CB") {

        delta[i] <- CV2B(n.f, mu.CB.f)  - CV2B(n.cf, mu.CB.cf)

      } else if(what=="I") {

        delta[i] <- CV2B(n.f, mu.f)  - CV2B(n.cf, mu.cf) - (CV2B(n.f, mu.CX.f)  - CV2B(n.cf, mu.CX.cf) + CV2B(n.f, mu.CB.f)  - CV2B(n.cf, mu.CB.cf))

      }


    } else if(ystat=="VarB") {

      dydx_t1 <- VarB(n.cf, mu.cf + beta.star) - VarB(n.cf, mu.cf)
      dydx_t0 <- VarB(n.cf, mu.cf + beta.cf)   - VarB(n.cf, mu.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="CV2W") {

      dydx_t1 <- CV2W(n.cf, mu.cf + beta.star, sigma.cf + lambda.star) - CV2W(n.cf, mu.cf, sigma.cf)
      dydx_t0 <- CV2W(n.cf, mu.cf + beta.cf, sigma.cf + lambda.cf)     - CV2W(n.cf, mu.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    } else if(ystat=="VarW") {

      dydx_t1 <- VarW(n.cf, sigma.cf + lambda.star) - VarW(n.cf, sigma.cf)
      dydx_t0 <- VarW(n.cf, sigma.cf + lambda.cf)   - VarW(n.cf, sigma.cf)

      delta[i] <- dydx_t1 - dydx_t0

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdXdD
# This function manipulates Mu, Sigma, and N in their effect through X
# ================================================================================================ #

dYdXdD <- function(dat, ystat, partial=F) {

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
# Function dYdD
# This function manipulates Mu, Sigma, and N
# ================================================================================================ #

dYdD <- function(dat, ystat, partial=F) {

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
# Function dYdMuSigma
# This function only manipulates Mu and Sigma, not N
# ================================================================================================ #

dYdMuSigma <- function(dat, ystat, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the impact of mu and sigma when no x is specified

    if(partial==T) {

      # Partial or total effect?

      mu.star    <- mu.cf
      mu.star[i] <- mu.f[i]

      if(wb=="W") sigma.star    <- sigma.cf
      if(wb=="W") sigma.star[i] <- sigma.f[i]

    } else {

      mu.star <- mu.f
      if(wb=="W") sigma.star <- sigma.f

    }

    if(ystat=="VarW") {

      delta[i] <- VarW(n.f, sigma.star) - VarW(n.f, sigma.cf)

    } else if(ystat=="VarB") {

      delta[i] <- VarB(n.f, mu.star) - VarB(n.f, mu.cf)

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n.f, mu.star, sigma.star) - CV2W(n.f, mu.cf, sigma.cf)

    } else if(ystat=="CV2B") {

      delta[i] <- CV2B(n.f, mu.star) - CV2B(n.f, mu.cf)

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdN
# This function manipulates only N, not Mu and Sigma
# ================================================================================================ #

dYdN <- function(dat, ystat, partial=F) {

  delta <- c() # container
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the impact of mu and sigma when no x is specified

    if(partial==T) {

      # Partial or total effect?

      n.star     <- n.cf
      n.star[i]  <- n.f[i]

    } else {

      n.star  <- n.f

    }

    if(ystat=="VarW") {

      delta[i] <- VarW(n.star, sigma.f) - VarW(n.cf, sigma.f)

    } else if(ystat=="VarB") {

      delta[i] <- VarB(n.star, mu.f) - VarB(n.cf, mu.f)

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n.star, mu.f, sigma.f) - CV2W(n.cf, mu.f, sigma.f)

    } else if(ystat=="CV2B") {

      delta[i] <- CV2B(n.star, mu.f) - CV2B(n.cf, mu.f)

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
