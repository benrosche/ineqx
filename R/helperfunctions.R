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

calcAME <- function(nox, xpv, what, vfr, dat) {

  if(!nox) {

    # Data with x = xpv[1] and x = xpv[2]
    newdat.a <- dat %>% dplyr::mutate(x=xpv[1])
    newdat.b <- dat %>% dplyr::mutate(x=xpv[2])

    # Predictions from variance function regression
    if(what=="mu") {
      pred <-
        dat %>%
        dplyr::select(time, group, x) %>%
        dplyr::mutate(pa=predict(vfr, what = "mu", type = "response", newdata = newdat.a, data = dat)) %>%
        dplyr::mutate(pb=predict(vfr, what = "mu", type = "response", newdata = newdat.b, data = dat))
    } else if(what=="sigma") {
      pred <-
        dat %>%
        dplyr::select(time, group, x) %>%
        dplyr::mutate(pa=predict(vfr, what = "sigma", type = "response", newdata = newdat.a, data = dat)) %>%
        dplyr::mutate(pb=predict(vfr, what = "sigma", type = "response", newdata = newdat.b, data = dat))
    }

    # Calculate Average Marginal Effects

    AME <- list()

    # Summarize by group and time
    AME[[1]] <-
      pred %>%
      group_by(time, group, x) %>%
      dplyr::summarise(pa=mean(pa), pb=mean(pb), beta=pb-pa) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(-x, -pa, -pb)

    # Summarize by time
    AME[[2]] <-
      pred %>%
      group_by(time, x) %>%
      dplyr::summarise(pa=mean(pa), pb=mean(pb), beta=pb-pa) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(-x, -pa, -pb)

  } else {

    AME <- list()

    AME[[1]] <-
      dat %>%
      dplyr::select(time, group) %>%
      group_by(time, group) %>%
      dplyr::summarise(beta=0)

    AME[[2]] <-
      dat %>%
      dplyr::select(time) %>%
      group_by(time) %>%
      dplyr::summarise(beta=0)

  }

  # Rename beta to lambda if what is sigma
  if(what=="sigma") AME <- lapply(AME, FUN=function(x) x %>% dplyr::rename(lambda=beta))

  return(AME)

}


# ================================================================================================ #
# Function createCF
# ================================================================================================ #

createCF <- function(x, y, ystat, groupvar, timevar, cf, wibe.xpv1, AME_mu, AME_sigma, dat) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1987, "w") | cf=list(n=c(), mu=c(), sigma=c(), beta=c(), lambda=c())
  # ---------------------------------------------------------------------------------------------- #

  # Levels of group and time var
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as effectDat
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # Counterfactual reference --------------------------------------------------------------------- #

  if(length(cf)==5) {

    # manual values

    n <- cf$n
    mu <- cf$mu
    sigma <- cf$sigma
    beta <- cf$beta
    lambda <- cf$lambda

  } else if(length(cf)==2) {

    # values from time point, gender specific

    n <- wibe.xpv1 %>% dplyr::filter(time==!!cf[1]) %>% .$n
    mu <- wibe.xpv1 %>% dplyr::filter(time==!!cf[1]) %>% .$mu
    sigma <- wibe.xpv1 %>% dplyr::filter(time==!!cf[1]) %>% .$sigma
    beta <- AME_mu[[1]] %>% dplyr::filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    lambda <- AME_sigma[[1]] %>% dplyr::filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()

  }  else if(length(cf)==1) {

    # values from time point

    n <- wibe.xpv1 %>% dplyr::filter(time==!!cf) %>% .$n
    mu <- wibe.xpv1 %>% dplyr::filter(time==!!cf) %>% .$mu
    sigma <- wibe.xpv1 %>% dplyr::filter(time==!!cf) %>% .$sigma
    beta <- AME_mu[[1]] %>% dplyr::filter(time==!!cf) %>% .$beta %>% as.vector()
    lambda <- AME_sigma[[1]] %>% dplyr::filter(time==!!cf) %>% .$lambda %>% as.vector()

  } else stop("cf must either be one value (reference time), two values (reference time + ...), or a list of 5 vectors (n, mu, sigma, beta, lambda)")

  # Create counterfactual dataset ---------------------------------------------------------------- #

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
        AME_mu[[1]] %>%
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
        AME_sigma[[1]] %>%
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
    inner_join(wibe.xpv1 %>% dplyr::select(time, group, n, mu, sigma), by=c("time", "group")) %>%
    left_join(AME_mu[[1]], by=c("time", "group")) %>%
    left_join(AME_sigma[[1]], by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f =  mu, sigma.f = sigma, beta.f = beta, lambda.f = lambda)

  return(dat.f_cf)

}


# ================================================================================================ #
# Function dWB
# ================================================================================================ #

dWB <- function(nox, ystat, dat.f_cf) {

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
      unnest(cols = c(data, dB, dD)) %>%
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

dT <- function(nox, dW.out, dB.out, ystat, cf) {

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

  # Add actual change in inequality
  if(length(cf)==5) {
    if(ystat == "CV2") {
      total <- total %>% dplyr::mutate(dCV2T=CV2T-(CV2W(cf$n, cf$mu, cf$sigma)+CV2B(cf$n, cf$mu)))
    } else {
      total <- total %>% dplyr::mutate(dVarT=VarT-(VarW(cf$n, cf$sigma)+VarB(cf$n, cf$mu)))
    }
  } else {
    total <- if(ystat == "CV2") total %>% dplyr::mutate(dCV2T=CV2T-CV2T[time == cf[1]]) else  total %>% dplyr::mutate(dVarT=VarT-VarT[time == cf[1]])
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
