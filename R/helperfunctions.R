
# ================================================================================================ #
# Function dB
# ================================================================================================ #

dB <- function(dat, pred_mu, x, y, ystat="CV2", groupvar, timevar, cf, smoothDat=F) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1, "beta_w")
  # ---------------------------------------------------------------------------------------------- #

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})
  pred_mu <- pred_mu %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }})

  # Levels of group and time
  group_levels <- pred_mu %>% .$group %>% unique() %>% sort()
  time_levels  <- pred_mu %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as pred_mu
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # wibe() --------------------------------------------------------------------------------------- #

  # dat.between must source pre-birth n and mu
  dat.between <-
    tibble() %>% expand(time=min(time_levels):max(time_levels)) %>%
    left_join(wibe(dat %>% dplyr::filter(x == 0), y, group, time, smoothDat = smoothDat, rel = F)[[1]], by="time")

  # dat.total must source post-birth CV2B
  dat.total <-
    tibble() %>% expand(time=min(time_levels):max(time_levels)) %>%
    left_join(wibe(dat %>% filter(x == 1), y, group, time, smoothDat = smoothDat, rel = F)[[2]], by="time")

  # Gather counterfactual data ------------------------------------------------------------------- #

  n  <- dat.between %>% filter(time==!!cf) %>% .$n

  mu <- dat.between %>% filter(time==!!cf) %>% .$mu

  beta <- (
    if(length(cf)==2) pred_mu %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else pred_mu %>% filter(time==!!cf) %>% .$beta %>% as.vector()
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
        pred_mu %>%
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
    left_join(pred_mu, by=c("time", "group")) %>%
    inner_join(dat.between %>% dplyr::select(time, group, n, mu), by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f =  mu, beta.f = beta)

  # Calculate impact ----------------------------------------------------------------------------- #

  if(ystat=="CV2") {

    # Calculates the partial change in CV2B by group and time #
    impact.between <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~deltaCV2B is executed by year
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.f) - deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.cf))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaCV2B(.x$n.f, .x$mu.f, .x$beta.f) - deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.f))) %>%
      dplyr::mutate(dO=purrr::map(.x = data, ~ CV2B(.x$n.f, .x$mu.f) - CV2B(.x$n.cf, .x$mu.cf))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::select(time, group, dB, dD, dO)

    # Calculates the total change in CV2B by time #
    impact.between.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.f, partial=F) - deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.cf, partial=F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaCV2B(.x$n.f, .x$mu.f, .x$beta.f, partial=F) - deltaCV2B(.x$n.cf, .x$mu.cf, .x$beta.f, partial=F))) %>%
      dplyr::mutate(dO=purrr::map(.x = data, ~ CV2B(.x$n.f, .x$mu.f) - CV2B(.x$n.cf, .x$mu.cf))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::select(time, dB, dD, dO)

  } else if (ystat=="Var") {

    # Calculates the partial change in VarB by SES #
    impact.between <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.f) - deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.cf))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaSigma2B(.x$n.f, .x$mu.f, .x$beta.f) - deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.f))) %>%
      dplyr::mutate(dO=purrr::map(.x = data, ~ VarB(.x$n.f, .x$mu.f) - VarB(.x$n.cf, .x$mu.cf))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::select(time, group, dB, dD, dO)

    # Calculates the total change in VarB #
    impact.between.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dB = purrr::map(.x = data, ~ deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.f, partial = F) - deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.cf, partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaSigma2B(.x$n.f, .x$mu.f, .x$beta.f, partial = F) - deltaSigma2B(.x$n.cf, .x$mu.cf, .x$beta.f, partial = F))) %>%
      dplyr::mutate(dO=purrr::map(.x = data, ~ VarB(.x$n.f, .x$mu.f) - VarB(.x$n.cf, .x$mu.cf))) %>%
      unnest(cols = c(data, dB, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::select(time, dB, dD, dO)

  }

  # Rename variables back to their original names ------------------------------------------------ #

  impact.between <-
    impact.between %>%
    inner_join(dat.total %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
      !!enquo(groupvar) := group,
    )

  impact.between.total <-
    impact.between.total %>%
    inner_join(dat.total %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
    )

  return(list(impact.between, impact.between.total))

}

# ================================================================================================ #
# Function dW
# ================================================================================================ #

dW <- function(dat, pred_mu, pred_sigma, x, y, ystat, groupvar, timevar, cf, smoothDat=F) {

  # ---------------------------------------------------------------------------------------------- #
  # Function arguments
  # cf = 1987 | cf = c(1987, "w")
  # ---------------------------------------------------------------------------------------------- #

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})
  pred_mu <- pred_mu %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }})
  pred_sigma <- pred_sigma %>% dplyr::rename(group={{ groupvar }}, time={{ timevar }})

  # Levels of group and time var
  group_levels <- pred_sigma %>% .$group %>% unique() %>% sort()
  time_levels  <- pred_sigma %>% .$time %>% unique() %>% sort()

  # Take subset of dat that has the same time window as effectDat
  dat <- dat %>% filter(between(time, min(time_levels), max(time_levels)))

  # wibe() --------------------------------------------------------------------------------------- #

  # dat.between must source pre-birth n, mu, and sigma2
  dat.between <-
    tibble() %>% expand(time=min(time_levels):max(time_levels)) %>%
    left_join(wibe(dat %>% dplyr::filter(x == 0), y, group, time, smoothDat = smoothDat, rel = F)[[1]], by="time")

  # dat.total must source post-birth CV2
  dat.total <-
    tibble() %>% expand(time=min(time_levels):max(time_levels)) %>%
    left_join(wibe(dat %>% filter(x == 1), y, group, time, smoothDat = smoothDat, rel = F)[[2]], by="time")

  # Gather counterfactual data ------------------------------------------------------------------- #

  n <- dat.between %>% filter(time==!!cf) %>% .$n

  mu <- dat.between %>% filter(time==!!cf) %>% .$mu

  sigma <- dat.between %>% filter(time==!!cf) %>% .$sigma

  beta <- (
    if(length(cf)==2) pred_mu %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else pred_mu %>% filter(time==!!cf) %>% .$beta %>% as.vector()
  )

  lambda <- (
    if(length(cf)==2) pred_sigma %>% filter(time==!!cf[1]) %>% dplyr::select(!!cf[2]) %>% unlist()
    else pred_sigma %>% filter(time==!!cf) %>% .$lambda %>% as.vector()
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
        pred_mu %>%
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
        pred_sigma %>%
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
    left_join(pred_mu, by=c("time", "group")) %>%
    left_join(pred_sigma, by=c("time", "group")) %>%
    inner_join(dat.between %>% dplyr::select(time, group, n, mu, sigma), by=c("time", "group")) %>%
    inner_join(dat.cf, by=c("time", "group")) %>%
    dplyr::rename(n.f = n, mu.f =  mu, sigma.f = sigma, beta.f = beta, lambda.f = lambda)

  # Calculate impact ----------------------------------------------------------------------------- #

  if(ystat=="CV2") {

    # Calculates the partial change in CV2W by group and time #
    impact.within <-
      dat.f_cf %>%
      group_by(time) %>% # so that ~deltaCV2W is executed by year
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.f, .x$lambda.f) - deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.cf, .x$lambda.cf))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaCV2W(.x$n.f,.x$mu.f, .x$sigma.f, .x$beta.f, .x$lambda.f) - deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.f, .x$lambda.f))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ CV2W(.x$n.f, .x$mu.f, .x$sigma.f) - CV2W(.x$n.cf, .x$mu.cf, .x$sigma.cf))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::select(time, group, dW, dD, dO)

    # Calculates the total change in CV2W by time #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.f, .x$lambda.f, partial=F) - deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.cf, .x$lambda.cf, partial=F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaCV2W(.x$n.f,.x$mu.f, .x$sigma.f, .x$beta.f, .x$lambda.f, partial=F) - deltaCV2W(.x$n.cf,.x$mu.cf, .x$sigma.cf, .x$beta.f, .x$lambda.f, partial=F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ CV2W(.x$n.f, .x$mu.f, .x$sigma.f) - CV2W(.x$n.cf, .x$mu.cf, .x$sigma.cf))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::select(time, dW, dD, dO)

  } else if (ystat=="Var") {

    # Calculates the partial change in VarW by SES #
    impact.within <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.f) - deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.cf))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaSigma2W(.x$n.f, .x$sigma.f, .x$lambda.f) - deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.f))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ VarW(.x$n.f, .x$sigma.f) - VarW(.x$n.cf, .x$sigma.cf))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::select(time, group, dW, dD, dO)

    # Calculates the total change in VarW #
    impact.within.total <-
      dat.f_cf %>%
      group_by(time) %>%
      nest() %>%
      dplyr::mutate(dW = purrr::map(.x = data, ~ deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.f, partial = F) - deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.cf, partial = F))) %>%
      dplyr::mutate(dD = purrr::map(.x = data, ~ deltaSigma2W(.x$n.f, .x$sigma.f, .x$lambda.f, partial = F) - deltaSigma2W(.x$n.cf, .x$sigma.cf, .x$lambda.f, partial = F))) %>%
      dplyr::mutate(dO = purrr::map(.x = data, ~ VarW(.x$n.f, .x$sigma.f) - VarW(.x$n.cf, .x$sigma.cf))) %>%
      unnest(cols = c(data, dW, dD, dO)) %>%
      dplyr::filter(row_number()==1) %>%
      dplyr::select(time, dW, dD, dO)

  }

  # Rename variables back to their original names ------------------------------------------------ #

  impact.within <-
    impact.within %>%
    inner_join(dat.total %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
      !!enquo(groupvar) := group,
    )

  impact.within.total <-
    impact.within.total %>%
    inner_join(dat.total %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
    ungroup() %>%
    dplyr::rename(
      !!enquo(timevar) := time,
    )

  return(list(impact.within, impact.within.total))

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

  return(list(tibble(dB.out[[1]][,1:2],ret[[1]]), tibble(dB.out[[2]][,1:2], ret[[2]])))

}

# ================================================================================================ #
# Function dT
# ================================================================================================ #

dT <- function(dW.out, dB.out, ystat="CV2") {

  ystatvar <- as.symbol(paste0(ystat, "T"))

  timevar <- names(dW.out[[1]])[1]

  # Combine effects ------------------------------------------------------------------------------ #

  total <-
    dW.out[[2]] %>%
    dplyr::rename(dD.W=dD, dO.W=dO) %>%
    inner_join(
      dB.out[[2]] %>%
        dplyr::rename(dD.B=dD, dO.B=dO),
      by=c("year", paste0(ystat, "T"))) %>%
    dplyr::mutate(
      dD=dD.B+dD.W,
      dT=dB+dW+dD,
      dO=dO.B+dO.W
    ) %>%
    dplyr::select(-dD.W, -dD.B, -dO.W, -dO.B) %>%
    dplyr::relocate(paste0(ystat, c("W", "B", "T")), .after = last_col())

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
# Function deltaSigma2W
# ================================================================================================ #

# Change of within-group variance
deltaSigma2W <- function(n, sigma, lambda, partial=T) {

  delta <- c()

  l <- length(n)

  for(i in 1:l) {

    # Partial or total effect
    if(partial==T) {
      lambda.star <- rep(0, l)
      lambda.star[i] <- lambda[i]
    } else {
      lambda.star <- lambda
    }

    delta[i] <-   VarW(n, sigma+lambda.star) - VarW(n, sigma)
  }

  return(as.numeric(delta))

}

# Calculates the within-group variance
VarW <- function(n, sigma) {

  w <- n/sum(n)

  return (sum(w * sigma^2))
}

# ================================================================================================ #
# Function deltaSigma2B
# ================================================================================================ #

# Change of between-group variance
deltaSigma2B <- function(n, mu, beta, partial=T) {

  delta <- c()

  l <- length(n)

  for(i in 1:l) {

    # Partial or total effect
    if(partial==T) {
      beta.star <- rep(0, l)
      beta.star[i] <- beta[i]
    } else {
      beta.star <- beta
    }

    delta[i] <-  VarB(n, mu+beta.star) - VarB(n, mu)

  }

  return(as.numeric(delta))

}

# Calculates the between-group variance
VarB <- function(n, mu) {

  w <- n/sum(n)
  gmu <- sum(w * mu)

  return( sum(w*(mu-gmu)^2) )
}

# ================================================================================================ #
# Function deltaCV2B
# ================================================================================================ #

# Change of between-group squared coefficient of variation
deltaCV2B<- function(n, mu, beta, partial=T) {

  # The function does not differentiate by year!

  # Check each lambda for its effect on the CV2B and fix the others at 0
  delta <- c()
  l <- length(n)

  for(i in 1:l) {

    # Partial or total effect
    if(partial==T) {
      beta.star <- rep(0, l)
      beta.star[i] <- beta[i]
    } else {
      beta.star <- beta
    }

    delta[i] <-  CV2B(n, mu+beta.star) - CV2B(n, mu)

  }

  return(as.numeric(delta))
}

# Calculates the between-group squared coefficient of variation
CV2B <- function(n, mu) {

  w <- n/sum(n)
  gmu <- sum(w * mu)

  return( sum(w * (mu-gmu)^2)/gmu^2 )
}

# ================================================================================================ #
# Function deltaCV2W
# ================================================================================================ #

# Change of within-group squared coefficient of variation
deltaCV2W <- function(n, mu, sigma, beta, lambda, partial=T) {

  # The function does not differentiate by year!

  delta <- c()
  l <- length(n)

  for(i in 1:l) {

    # Partial or total effect
    if(partial==T) {

      lambda.star <- rep(0, l)
      lambda.star[i] <- lambda[i]

      beta.star <- rep(0, l)
      beta.star[i] <- beta[i]

    } else {
      lambda.star <- lambda
      beta.star <- beta
    }

    delta[i] <- CV2W(n, mu+beta.star, sigma+lambda.star) - CV2W(n, mu, sigma)

  }

  return(delta)

}

# Calculates the within-group squared coefficient of variation
CV2W <- function(n, mu, sigma) {

  w <- n/sum(n)

  return( (sum(w * sigma^2)) / (sum(w * mu))^2 )
}
