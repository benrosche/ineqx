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

    return(list("var"=x, "cont"=xcont, "poly"=xpoly, "novar"=F))

  } else {

    return(list("var"=NULL, "cont"=NULL, "poly"=NULL, "novar"=T))

  }

}


# ================================================================================================ #
# Function createCF
# ================================================================================================ #

calc01 <- function(group, time, ref, wibe.pre, AME_mu, AME_sigma, notime, dat) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # ref = 1987 | ref = c(1987, "effect1") | ref=list(n=c(), mu=c(), sigma=c(), beta=c(), lambda=c())
  # ---------------------------------------------------------------------------------------------- #

  # Levels of group and time var
  group_levels <- dat %>% pull(group) %>% unique() %>% sort()
  time_levels  <- dat %>% pull(time) %>% unique() %>% sort()

  # Take subset of dat that has the same time window as effectDat
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # If notime, we take pre as reference
  if(notime) {
    ref <- list(n=wibe.pre %>% pull(n), mu=wibe.pre %>% pull(mu), sigma=wibe.pre %>% pull(sigma), beta=rep(0, length(group_levels)), lambda=rep(0, length(group_levels)))
  }

  # Get reference values ------------------------------------------------------------------------- #

  if(is.null(names(ref)) & between(length(ref), 1, 2)) {

    refm <- FALSE # time reference

    if(length(ref)==1) ref <- c(ref, "effect")

    n <- wibe.pre %>% dplyr::filter(time==!!ref[1]) %>% pull(n)
    mu <- wibe.pre %>% dplyr::filter(time==!!ref[1]) %>% pull(mu)
    sigma <- wibe.pre %>% dplyr::filter(time==!!ref[1]) %>% pull(sigma)
    beta <- AME_mu[[1]] %>% dplyr::filter(time==!!ref[1]) %>% pull(!!ref[2])
    lambda <- AME_sigma[[1]] %>% dplyr::filter(time==!!ref[1]) %>% pull(!!ref[2])

  } else if(!is.null(names(ref)) & all(names(ref) %in% c("n", "mu", "sigma", "beta", "lambda"))) {

    refm <- TRUE # manual values as reference

    n      <- if(any(names(ref) %in% "n")) ref$n else NA
    mu     <- if(any(names(ref) %in% "mu")) ref$mu else NA
    sigma  <- if(any(names(ref) %in% "sigma")) ref$sigma else NA
    beta   <- if(any(names(ref) %in% "beta")) ref$beta else NA
    lambda <- if(any(names(ref) %in% "lambda")) ref$lambda else NA

    ref <- c(0, "effect")

  } else stop("ref must either be a reference time (numeric value) (recommended), a reference time of a specific column, e.g. c(1987, 'effect1') (experimental), or a list of vectors that can be either or all of: list(n=c(), mu=c(), sigma=c(), beta=c(), lambda=c()).")

  # Create counterfactual dataset ---------------------------------------------------------------- #

  dat0 <-
    tibble() %>%
    tidyr::expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(tibble(group=group_levels, n0=n), by="group") %>%
    inner_join(tibble(group=group_levels, mu0=mu), by="group") %>%
    inner_join(tibble(group=group_levels, sigma0=sigma), by="group") %>%
    inner_join(
      AME_mu[[1]] %>%
        dplyr::select(-!!ref[2]) %>% # replace ref[2] ...
        inner_join(tibble(group=group_levels, !!ref[2]:=beta), by="group") %>% # ... values at t=0
        dplyr::mutate(beta0 = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% # and then calculate beta0 as sum over all columns
        dplyr::select(time, group, beta0),
      by=c("time", "group")) %>%
    inner_join(
      AME_sigma[[1]] %>%
        dplyr::select(-!!ref[2]) %>%
        inner_join(tibble(group=group_levels, !!ref[2]:=lambda), by="group") %>%
        dplyr::mutate(lambda0 = rowSums(dplyr::across(c(everything(), -time, -group)))) %>%
        dplyr::select(time, group, lambda0),
      by=c("time", "group"))

  # Re: Beta and Lambda
  # I allow AME_mu/AME_sigma to contain several betas that together sum to the total beta
  # so that I can calculate the effect of changing only one column
  # If there is just one column, the procedure is equivalent to
  # inner_join(tibble(group=group_levels, beta0=beta), by="group")

  # Merge factual and counterfactual data -------------------------------------------------------- #

  dat01 <-
    tibble() %>%
    tidyr::expand(time=min(time_levels):max(time_levels), group=group_levels) %>%
    inner_join(wibe.pre %>% dplyr::select(time, group, n, mu, sigma) %>% dplyr::rename(n1 = n, mu1 =  mu, sigma1 = sigma), by=c("time", "group")) %>% # factual n, mu, sigma
    left_join(AME_mu[[1]] %>% dplyr::mutate(effect = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% dplyr::rename(beta1=effect), by=c("time", "group")) %>% # factual beta
    left_join(AME_sigma[[1]] %>% dplyr::mutate(effect = rowSums(dplyr::across(c(everything(), -time, -group)))) %>% dplyr::rename(lambda1=effect), by=c("time", "group")) %>% # factual lambda
    inner_join(dat0, by=c("time", "group")) %>% # counterfactual data
    dplyr::select(time, group, ends_with("1"), ends_with("0")) %>%
    dplyr::mutate(
      n0 = case_when(is.na(n0) ~ as.numeric(n1), TRUE ~ as.numeric(n0)),
      mu0 = case_when(is.na(mu0) ~ as.numeric(mu1), TRUE ~ as.numeric(mu0)),
      sigma0 = case_when(is.na(sigma0) ~ as.numeric(sigma1), TRUE ~ as.numeric(sigma0)),
      beta0 = case_when(is.na(beta0) ~ as.numeric(beta1), TRUE ~ as.numeric(beta0)),
      lambda0 = case_when(is.na(lambda0) ~ as.numeric(lambda1), TRUE ~ as.numeric(lambda0))
    )

  # Re last mutate: If manual values are only provided to a subset of the parameters,
  #                 then the other parameters take on the values at t=1

  if(refm) {

    # If ref are manual values, add zero row with those values

    row_zero <-
      dat01 %>%
      dplyr::filter(time==min(time)) %>%
      dplyr::mutate(time=0, n1=n0, mu1=mu0, sigma1=sigma0, beta1=beta0, lambda1=lambda0)

    dat01 <-
      dat01 %>%
      add_row(row_zero) %>%
      arrange(time)

  }

  return(dat01)

}


# ================================================================================================ #
# Function deltaWB
# ...
# ================================================================================================ #

deltaWB <- function(ystat, wb, decomp, notreat, dat01) {

  # The changing effect of treat on within or between-group inequality --------------------------- #

  wb <- toupper(wb)

  if(!notreat) {

    # Partial change = Effect of treat for one group keeping the effect of treat of other groups at 0
    impact.partial <-
      dat01 %>%
      group_by(time) %>% # so that ~ fcts are executed by year
      tidyr::nest() %>%
      dplyr::mutate("d{wb}" := purrr::map(.x = data, ~ dYdX(ystat=paste0(ystat,wb), partial = T, .x))) %>%
      dplyr::mutate(dC = purrr::map(.x = data, ~ dYdXdC(ystat=paste0(ystat,wb), decomp=decomp, partial = T, .x))) %>%
      dplyr::mutate(dP = purrr::map(.x = data, ~ dYdXdP(ystat=paste0(ystat,wb), decomp=decomp, partial = T, .x))) %>%
      tidyr::unnest(cols = c(data, starts_with("d"))) %>%
      ungroup() %>%
      dplyr::select(time, group, starts_with("d"))

    # Total change = Total effect of treat
    impact.total <-
      dat01 %>%
      group_by(time) %>%
      tidyr::nest() %>%
      dplyr::mutate("d{wb}" := purrr::map(.x = data, ~ dYdX(ystat=paste0(ystat,wb), partial = F, .x))) %>%
      dplyr::mutate(dC = purrr::map(.x = data, ~ dYdXdC(ystat=paste0(ystat,wb), decomp=decomp, partial = F, .x))) %>%
      dplyr::mutate(dP = purrr::map(.x = data, ~ dYdXdP(ystat=paste0(ystat,wb), decomp=decomp, partial = F, .x))) %>%
      tidyr::unnest(cols = c(data, starts_with("d"))) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, starts_with("d"))

  } else {

    # The changing effect of Mu and Sigma ---------------------------------------------------------- #

    impact.partial <-
      dat01 %>%
      group_by(time) %>%
      tidyr::nest() %>%
      dplyr::mutate("d{wb}" := purrr::map(.x = data, ~ dYdP(ystat=paste0(ystat,wb), partial = T, .x))) %>%
      dplyr::mutate(dC = purrr::map(.x = data, ~ dYdC(ystat=paste0(ystat,wb), partial = T, .x))) %>%
      tidyr::unnest(cols = c(data, starts_with("d"))) %>%
      ungroup() %>%
      dplyr::select(time, group, starts_with("d"))

    impact.total <-
      dat01 %>%
      group_by(time) %>%
      tidyr::nest() %>%
      dplyr::mutate("d{wb}" := purrr::map(.x = data, ~ dYdP(ystat=paste0(ystat,wb), partial = F, .x))) %>%
      dplyr::mutate(dC = purrr::map(.x = data, ~ dYdC(ystat=paste0(ystat,wb), partial = F, .x))) %>%
      tidyr::unnest(cols = c(data, starts_with("d"))) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(time, starts_with("d"))
  }

  return(list(impact.partial, impact.total))

}


# ================================================================================================ #
# Function deltaC
# ================================================================================================ #

deltaC <- function(notreat, deltaW.out, deltaB.out) {

  ret <-
    purrr::map2(
      .x=deltaB.out, .y=deltaW.out,
      ~ .x %>%
        dplyr::select(dC, matches("VarT|CV2T")) %>%
        dplyr::mutate(dC=.x$dC+.y$dC) %>%
        ungroup()
    )

  return(list(tibble(deltaB.out[[1]][,1:2],ret[[1]]), tibble(deltaB.out[[2]][,1], ret[[2]])))

}


# ================================================================================================ #
# Function deltaP
# ================================================================================================ #

deltaP <- function(notreat, deltaW.out, deltaB.out) {

  ret <-
    purrr::map2(
      .x=deltaB.out, .y=deltaW.out,
      ~ .x %>%
        dplyr::select(dP, matches("VarT|CV2T")) %>%
        dplyr::mutate(dP=.x$dP+.y$dP) %>%
        ungroup()
      )

  return(list(tibble(deltaB.out[[1]][,1:2],ret[[1]]), tibble(deltaB.out[[2]][,1], ret[[2]])))

}


# ================================================================================================ #
# Function deltaT
# ================================================================================================ #

deltaT <- function(notreat, deltaW.out, deltaB.out) {

  # Combine effects ------------------------------------------------------------------------------ #

  if(!notreat) {

    total <-
      deltaW.out[[2]] %>%
      dplyr::rename(dC.W=dC, dP.W=dP) %>%
      inner_join(
        deltaB.out[[2]] %>%
          dplyr::rename(dC.B=dC, dP.B=dP) %>%
          dplyr::select(-matches("VarT|CV2T")),
        by=c("time")) %>%
      dplyr::mutate(
        dC=dC.W+dC.B,
        dP=dP.W+dP.B,
        dT=dW+dB+dC+dP
      ) %>%
      dplyr::select(time, dW, dB, dC, dP, dT, matches("VarW|CV2W"), matches("VarB|CV2B"), matches("VarT|CV2T"))

  } else {

    total <-
      deltaW.out[[2]] %>%
      dplyr::rename(dC.W=dC) %>%
      inner_join(
        deltaB.out[[2]] %>%
          dplyr::rename(dC.B=dC) %>%
          dplyr::select(-matches("VarT|CV2T")),
        by=c("time")) %>%
      dplyr::mutate(
        dC=dC.W+dC.B,
        dT=dW+dB+dC
      ) %>%
      dplyr::select(time, dW, dB, dC, dT, matches("VarW|CV2W"), matches("VarB|CV2B"), matches("VarT|CV2T"))

  }

  # Calculate shares ----------------------------------------------------------------------------- #

  shares <-
    total %>%
    dplyr::select(-matches("Var"), -matches("CV2"), -"dT") %>%
    dplyr::mutate(across(everything(), ~abs(.x))) %>%
    dplyr::mutate(abssum=rowSums(across(c(everything(), -time)))) %>%
    dplyr::mutate(across(c(everything(), -time), ~.x/abssum)) %>%
    dplyr::select(-abssum) %>%
    tidyr::pivot_longer(cols=-time, names_to="delta", values_to="share")

  return(list(total, shares))

}


# ================================================================================================ #
# Function dYdX
# ...
# ================================================================================================ #

dYdX <- function(ystat, partial=F, dat) {

  delta <- c() # container
  ystat <- toupper(ystat)
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the group-specific effect of X on the within- or between-group component

    if(isTRUE(partial)) {

      # Partial or total effect?

      beta.star    <- beta0
      beta.star[i] <- beta1[i]

      if(wb=="W") lambda.star    <- lambda0
      if(wb=="W") lambda.star[i] <- lambda1[i]

    } else {

      beta.star <- beta1
      if(wb=="W") lambda.star <- lambda1

    }

    if(ystat=="CV2B") {

      delta[i] <- CV2B(n0, mu0 + beta.star) - CV2B(n0, mu0) - ( CV2B(n0, mu0 + beta0) - CV2B(n0, mu0) )

    } else if(ystat=="VARB") {

      delta[i] <- VarB(n0, mu0 + beta.star) - VarB(n0, mu0) - ( VarB(n0, mu0 + beta0) - VarB(n0, mu0) )

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n0, mu0 + beta.star, sigma0 + lambda.star) - CV2W(n0, mu0, sigma0) - ( CV2W(n0, mu0 + beta0, sigma0 + lambda0) - CV2W(n0, mu0, sigma0) )

    } else if(ystat=="VARW") {

      delta[i] <- VarW(n0, sigma0 + lambda.star) - VarW(n0, sigma0) - ( VarW(n0, sigma0 + lambda0) - VarW(n0, sigma0) )

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdXdC
# This function manipulates n in its effect through X
# ================================================================================================ #

dYdXdC <- function(ystat, decomp, partial=F, dat) {

  delta <- c() # container
  ystat <- toupper(ystat)
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Partial or total effect?

    if(isTRUE(partial)) {

      n.star     <- n0
      n.star[i]  <- n1[i]

      mu.star    <- mu0
      mu.star[i] <- mu1[i]

      if(wb=="W") sigma.star    <- sigma0
      if(wb=="W") sigma.star[i] <- sigma1[i]

    } else {

      n.star  <- n1
      mu.star <- mu1
      if(wb=="W") sigma.star <- sigma1

    }

    # Calculates the group-specific effect of n realized through X

    if(ystat=="CV2B") {

      if(decomp=="effect") delta[i] <- CV2B(n.star, mu0 + beta1) - CV2B(n0, mu0 + beta1) - ( CV2B(n1, mu1) - CV2B(n0, mu1) )
      if(decomp=="post")   delta[i] <- CV2B(n.star, mu0 + beta1) - CV2B(n0, mu0 + beta1)

    } else if(ystat=="VARB") {

      if(decomp=="effect") delta[i] <- VarB(n.star, mu0 + beta1) - VarB(n0, mu0 + beta1) - ( VarB(n1, mu1) - VarB(n0, mu1) )
      if(decomp=="post")   delta[i] <- VarB(n.star, mu0 + beta1) - VarB(n0, mu0 + beta1)

    } else if(ystat=="CV2W") {

      if(decomp=="effect") delta[i] <- CV2W(n.star, mu0 + beta1, sigma0 + lambda1)  - CV2W(n0, mu0 + beta1, sigma0 + lambda1) - ( CV2W(n1, mu1, sigma1) - CV2W(n0, mu1, sigma1) )
      if(decomp=="post")   delta[i] <- CV2W(n.star, mu0 + beta1, sigma0 + lambda1)  - CV2W(n0, mu0 + beta1, sigma0 + lambda1)

    } else if(ystat=="VARW") {

      if(decomp=="effect") delta[i] <- VarW(n.star, sigma0 + lambda1) - VarW(n0, sigma0 + lambda1) - ( VarW(n1, sigma1) - VarW(n0, sigma1) )
      if(decomp=="post")   delta[i] <- VarW(n.star, sigma0 + lambda1) - VarW(n0, sigma0 + lambda1)

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdXdP
# This function manipulates Mu, Sigma in their effect through X
# ================================================================================================ #

dYdXdP <- function(ystat, decomp, partial=F, dat) {

  delta <- c() # container
  ystat <- toupper(ystat)
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Partial or total effect?

    if(isTRUE(partial)) {

      n.star     <- n0
      n.star[i]  <- n1[i]

      mu.star    <- mu0
      mu.star[i] <- mu1[i]

      if(wb=="W") sigma.star    <- sigma0
      if(wb=="W") sigma.star[i] <- sigma1[i]

    } else {

      n.star  <- n1
      mu.star <- mu1
      if(wb=="W") sigma.star <- sigma1

    }

    # Calculates the group-specific effect of Mu and Sigma that is realized through X

    if(ystat=="CV2B") {

      if(decomp=="effect") delta[i] <- CV2B(n1, mu.star + beta1) - CV2B(n1, mu0 + beta1) - ( CV2B(n0, mu.star) - CV2B(n0, mu0) )
      if(decomp=="post")   delta[i] <- CV2B(n1, mu.star + beta1) - CV2B(n1, mu0 + beta1)

    } else if(ystat=="VARB") {

      if(decomp=="effect") delta[i] <- VarB(n1, mu.star + beta1) - VarB(n1, mu0 + beta1) - ( VarB(n0, mu.star) - VarB(n0, mu0) )
      if(decomp=="post")   delta[i] <- VarB(n1, mu.star + beta1) - VarB(n1, mu0 + beta1)

    } else if(ystat=="CV2W") {

      if(decomp=="effect") delta[i] <- CV2W(n1, mu.star + beta1, sigma.star + lambda1) - CV2W(n1, mu0 + beta1, sigma0 + lambda1) - ( CV2W(n0, mu.star, sigma.star) - CV2W(n0, mu0, sigma0) )
      if(decomp=="post")   delta[i] <- CV2W(n1, mu.star + beta1, sigma.star + lambda1) - CV2W(n1, mu0 + beta1, sigma0 + lambda1)

    } else if(ystat=="VARW") {

      if(decomp=="effect") delta[i] <-  VarW(n1, sigma.star + lambda1) - VarW(n1, sigma0 + lambda1) - ( VarW(n0, sigma.star) - VarW(n0, sigma0) )
      if(decomp=="post")   delta[i] <-  VarW(n1, sigma.star + lambda1) - VarW(n1, sigma0 + lambda1)

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdC
# This function manipulates only n
# ================================================================================================ #

dYdC <- function(ystat, partial=F, dat) {

  delta <- c() # container
  ystat <- toupper(ystat)
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the impact of mu and sigma when no x is specified

    if(partial==T) {

      # Partial or total effect?

      n.star     <- n0
      n.star[i]  <- n1[i]

    } else {

      n.star  <- n1

    }

    if(ystat=="VARW") {

      delta[i] <- VarW(n.star, sigma1) - VarW(n0, sigma1)

    } else if(ystat=="VARB") {

      delta[i] <- VarB(n.star, mu1) - VarB(n0, mu1)

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n.star, mu1, sigma1) - CV2W(n0, mu1, sigma1)

    } else if(ystat=="CV2B") {

      delta[i] <- CV2B(n.star, mu1) - CV2B(n0, mu1)

    }

  }

  return(as.numeric(delta))

}


# ================================================================================================ #
# Function dYdP
# This function manipulates only Mu and Sigma
# ================================================================================================ #

dYdP <- function(ystat, partial=F, dat) {

  delta <- c() # container
  ystat <- toupper(ystat)
  wb <- substr(ystat, 4,4) # within/between
  l <- dim(dat)[1] # number of groups
  list2env(dat, envir=environment()) # extract variables from dataframe

  for(i in 1:l) {

    # Calculates the impact of mu and sigma when no x is specified

    if(partial==T) {

      # Partial or total effect?

      mu.star    <- mu0
      mu.star[i] <- mu1[i]

      if(wb=="W") sigma.star    <- sigma0
      if(wb=="W") sigma.star[i] <- sigma1[i]

    } else {

      mu.star <- mu1
      if(wb=="W") sigma.star <- sigma1

    }

    if(ystat=="VARW") {

      delta[i] <- VarW(n0, sigma.star) - VarW(n0, sigma0)

    } else if(ystat=="VARB") {

      delta[i] <- VarB(n0, mu.star) - VarB(n0, mu0)

    } else if(ystat=="CV2W") {

      delta[i] <- CV2W(n0, mu.star, sigma.star) - CV2W(n0, mu0, sigma0)

    } else if(ystat=="CV2B") {

      delta[i] <- CV2B(n0, mu.star) - CV2B(n0, mu0)

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

# eof
