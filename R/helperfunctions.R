# ================================================================================================ #
# Function dissectVar
# ================================================================================================ #

dissectVar <- function(x, cicheck="ci") {

  if(!is.null(x)) {

    if(cicheck=="ci" & !grepl("(c|i).?\\.", x)) stop(paste0(x," must be specified as either c.", x, " (continuous) or i.", x, " (categorical)."))
    if(cicheck=="i"  & grepl("(c).?\\.", x)) stop(paste0(x," must be i.", as.symbol(sub("(c|i).?\\.", "", x)), " (categorical)"))
    if(cicheck=="c"  & grepl("(i).?\\.", x)) stop(paste0(x," must be c.", as.symbol(sub("(c|i).?\\.", "", x)), " (continuous)"))

    xpoly <- suppressWarnings(as.numeric(substr(x, 2,2)))
    if(is.na(xpoly)) xpoly <- 1

    xcont <- grepl("c.?\\.", x)

    x <- as.symbol(sub("(c|i).?\\.", "", x))

    return(list("var"=x, "xcont"=xcont, "xpoly"=xpoly, "nox"=F))

  } else {

    return(list("var"=NULL, "xcont"=NULL, "xpoly"=NULL, "nox"=T))

  }

}


# ================================================================================================ #
# Function createCF
# ================================================================================================ #

createCF <- function(groupvar, timevar, ref, wibe.xpv0, AME_mu, AME_sigma, dat) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # ref = 1987 | ref = c(1987, "effect1") | ref=list(n=c(), mu=c(), sigma=c(), beta=c(), lambda=c())
  # ---------------------------------------------------------------------------------------------- #

  # Levels of group and time var
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as effectDat
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # If notime, we take 0 as manual reference
  if(wibe.xpv0 %>% .$time %>% unique() %>% length() == 1) {
    ref <- list(n=wibe.xpv0 %>% .$n, mu=wibe.xpv0 %>% .$mu, sigma=wibe.xpv0 %>% .$sigma, beta=c(0,0,0), lambda=c(0,0,0))
  }

  # Counterfactual reference --------------------------------------------------------------------- #

  if(length(ref)==5) {

    refm <- TRUE # manual values

    n <- ref$n
    mu <- ref$mu
    sigma <- ref$sigma
    beta <- ref$beta
    lambda <- ref$lambda

    ref <- c(0, "effect")

  } else if(between(length(ref),1,2)) {

    refm <- FALSE # values from reference time

    if(length(ref)==1) ref <- c(ref, "effect")

    n <- wibe.xpv0 %>% dplyr::filter(time==!!ref[1]) %>% .$n
    mu <- wibe.xpv0 %>% dplyr::filter(time==!!ref[1]) %>% .$mu
    sigma <- wibe.xpv0 %>% dplyr::filter(time==!!ref[1]) %>% .$sigma
    beta <- AME_mu[[1]] %>% dplyr::filter(time==!!ref[1]) %>% dplyr::select(!!ref[2]) %>% unlist()
    lambda <- AME_sigma[[1]] %>% dplyr::filter(time==!!ref[1]) %>% dplyr::select(!!ref[2]) %>% unlist()

  } else stop("ref must either be one value (reference time) (recommended), two values (reference time of a specific column), or a list of 5 vectors (n, mu, sigma, beta, lambda)")

  # Create counterfactual dataset ---------------------------------------------------------------- #

  dat.cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
    inner_join(tibble(group=group_levels, sigma.cf=sigma), by="group") %>%
    inner_join(
      AME_mu[[1]] %>%
        dplyr::select(-!!ref[2]) %>% # replace ref[2] ...
        inner_join(tibble(group=group_levels, !!ref[2]:=beta), by="group") %>% # ... with counterfactual values
        dplyr::mutate(beta.cf = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% # and then calculate beta.cf as sum over all columns
        dplyr::select(time, group, beta.cf)
      , by=c("time", "group")) %>%
    inner_join(
      AME_sigma[[1]] %>%
        dplyr::select(-!!ref[2]) %>%
        inner_join(tibble(group=group_levels, !!ref[2]:=lambda), by="group") %>%
        dplyr::mutate(lambda.cf = rowSums(dplyr::across(c(everything(), -time, -group)))) %>%
        dplyr::select(time, group, lambda.cf),
      by=c("time", "group"))

  # Re: Beta and Lambda
  # I allow AME_mu/AME_sigma to contain several betas that together sum to the total beta
  # so that I can calculate the effect of changing only one column
  # If there is just one column, the procedure is equivalent to
  # inner_join(tibble(group=group_levels, beta.cf=beta), by="group")

  # Merge factual and counterfactual data

  dat.f_cf <-
    tibble() %>%
    expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(wibe.xpv0 %>% dplyr::select(time, group, n, mu, sigma) %>% dplyr::rename(n.f = n, mu.f =  mu, sigma.f = sigma), by=c("time", "group")) %>% # factual n, mu, sigma
    left_join(AME_mu[[1]] %>% dplyr::mutate(effect = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% dplyr::rename(beta.f=effect), by=c("time", "group")) %>% # factual beta
    left_join(AME_sigma[[1]] %>% dplyr::mutate(effect = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% dplyr::rename(lambda.f=effect), by=c("time", "group")) %>% # factual lambda
    inner_join(dat.cf, by=c("time", "group")) %>% # counterfactual data
    dplyr::select(time, group, ends_with(".f"), ends_with(".cf"))

  if(refm) {

    # If ref are manual values, add zero row with those values

    row_zero <-
      tibble() %>%
      expand(time=0, group=group_levels) %>%
      inner_join(tibble(group=group_levels, n.f=n), by="group") %>%
      inner_join(tibble(group=group_levels, mu.f=mu), by="group") %>%
      inner_join(tibble(group=group_levels, sigma.f=sigma), by="group") %>%
      inner_join(tibble(group=group_levels, beta.f=beta), by="group") %>%
      inner_join(tibble(group=group_levels, lambda.f=lambda), by="group") %>%
      inner_join(tibble(group=group_levels, n.cf=n), by="group") %>%
      inner_join(tibble(group=group_levels, mu.cf=mu), by="group") %>%
      inner_join(tibble(group=group_levels, sigma.cf=sigma), by="group") %>%
      inner_join(tibble(group=group_levels, beta.cf=beta), by="group") %>%
      inner_join(tibble(group=group_levels, lambda.cf=lambda), by="group")

    dat.f_cf <-
      dat.f_cf %>%
      add_row(row_zero) %>%
      arrange(time)

  }

  return(dat.f_cf)

}


# ================================================================================================ #
# Function dX
# ================================================================================================ #

dX <- function(nox, ystat, dat.f_cf) {

  # The changing effect of X --------------------------------------------------------------------- #

  if(!nox) {

    # Partial change = Effect of X for one group keeping the effect of X of other groups at 0
    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~ fcts are executed by year
      nest() %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x,   ystat=ystat, partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdXdD(.x, ystat=ystat, partial = T))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdD(.x,   ystat=ystat, partial = T))) %>%
      unnest(cols = c(data, dX, dD, dO)) %>%
      ungroup() %>%
      dplyr::select(time, group, dX, dD, dO)

    # Total change = Total effect of X
    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdX(.x,   ystat=ystat, partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdXdD(.x, ystat=ystat, partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ dYdD(.x,   ystat=ystat, partial = F))) %>%
      unnest(cols = c(data, dX, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dX, dD, dO)

  } else {

    # The changing effect of Mu and Sigma ---------------------------------------------------------- #

    impact.partial <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=ystat, partial = T))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=ystat, partial = T))) %>%
      unnest(cols = c(data, dX, dD)) %>%
      ungroup() %>%
      dplyr::select(time, group, dX, dD)

    impact.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dX = purrr::map(.x = data, ~ dYdMuSigma(.x, ystat=ystat, partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ dYdN(.x,       ystat=ystat, partial = F))) %>%
      unnest(cols = c(data, dX, dD)) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, dX, dD)
  }

  return(list(impact.partial, impact.total))

}


# ================================================================================================ #
# Function dD
# ================================================================================================ #

dD <- function(nox, dW.out, dB.out) {

  if(!nox) {

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

dT <- function(nox, dW.out, dB.out, ystat) {

  # Combine effects ------------------------------------------------------------------------------ #

  if(!nox) {

    total <-
      dW.out[[2]] %>%
      dplyr::rename(dD.W=dD, dO.W=dO) %>%
      inner_join(
        dB.out[[2]] %>%
          dplyr::rename(dD.B=dD, dO.B=dO) %>%
          dplyr::select(-paste0(ystat, "T")),
        by=c("time")) %>%
      dplyr::mutate(
        dD=dD.W+dD.B,
        dT=dW+dB+dD,
        dO=dO.W+dO.B
      ) %>%
      dplyr::select(time, dW, dB, dD, dT, dO, paste0(ystat, c("W", "B", "T")))

  } else {

    total <-
      dW.out[[2]] %>%
      dplyr::rename(dD.W=dD) %>%
      inner_join(
        dB.out[[2]] %>%
          dplyr::rename(dD.B=dD) %>%
          dplyr::select(-paste0(ystat, "T")),
        by=c("time")) %>%
      dplyr::mutate(
        dD=dD.W+dD.B,
        dT=dW+dB+dD
      ) %>%
      dplyr::select(time, dW, dB, dD, dT, paste0(ystat, c("W", "B", "T")))

  }

  # Calculate shares ----------------------------------------------------------------------------- #

  shares <-
    total %>%
    dplyr::select(-matches("Var"), -matches("CV2"), -"dT") %>%
    dplyr::mutate(across(everything(), ~abs(.x))) %>%
    dplyr::mutate(abssum=rowSums(across(c(everything(), -time)))) %>%
    dplyr::mutate(across(c(everything(), -time), ~.x/abssum)) %>%
    dplyr::select(-abssum) %>%
    pivot_longer(cols=-time, names_to="d", values_to="share")

  return(list(total, shares))

}


# ================================================================================================ #
# Function dYdX
# ...
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
