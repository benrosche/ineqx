#' @title Calculate AME
#'
#' @description [...]
#'
#' @param treat Character string. Treatment variable.
#' @param group Grouping variable
#' @param time Time variable
#' @param what "Mu" or "Sigma"
#' @param vfr gamlss output
#' @param dat Dataframe
#'
#' @return List of length 2. Element 1 returns AME by group and time. Elements 2 returns the AME by time.
#'
#' @examples data(incdat)
#' vfr <- gamlss()
#' AME_mu <- calcAME("treat", group="SES", time="year", what="mu", vfr, incdat)
#'
#' @export calcAME
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

calcAME <- function(treat, group, time, what, vfr, dat) {

  # treat="tp"; group="group"; time="time"; what="mu"

  # ---------------------------------------------------------------------------------------------- #
  # Checks
  # ---------------------------------------------------------------------------------------------- #

  if(!as.character(treat) %in% names(dat)) stop(paste0(treat, " not in dataset"))
  if(!as.character(group) %in% names(dat)) stop(paste0(group, " not in dataset"))
  if(!as.character(time) %in% names(dat)) stop(paste0(time, " not in dataset"))

  treat_levels <- dat %>% pull(treat) %>% unique() %>% sort()
  if(!0 %in% treat_levels) stop(paste0("Values of ", treat, " must contain 0."))
  treat_levels <- treat_levels[treat_levels!=0]
  if(length(treat_levels)>1) warning("Effect of treat was calculated as weighted average of ", paste0("0->", treat_levels, sep=" "), ".")

  # ---------------------------------------------------------------------------------------------- #
  # Predictions from variance function regression
  # ---------------------------------------------------------------------------------------------- #

  if(what=="mu") {

    # y_hat for each observation if treat=0
    pred <-
      dat %>%
      dplyr::select({{ time }}, {{ group }} , {{ treat }}) %>%
      dplyr::mutate(pa=predict(vfr, what = "mu", type = "response", newdata = dat %>% dplyr::mutate({{ treat }} :=0), data = dat))

    # y_hat for each observation and each of the treat_levels
    for(i in 1:length(treat_levels)) {
      pred <-
        pred %>%
        dplyr::mutate("pb.{i}":=predict(vfr, what = "mu", type = "response", newdata = dat %>% dplyr::mutate({{ treat }} :=!!treat_levels[i]), data = dat))
    }

  } else if(what=="sigma") {

    # y_hat for each observation if treat=0
    pred <-
      dat %>%
      dplyr::select({{ time }}, {{ group }} , {{ treat }}) %>%
      dplyr::mutate(pa=predict(vfr, what = "sigma", type = "response", newdata = dat %>% dplyr::mutate({{ treat }} :=0), data = dat))

    # y_hat for each observation and each of the treat_levels
    for(i in 1:length(treat_levels)) {
      pred <-
        pred %>%
        dplyr::mutate("pb.{i}":=predict(vfr, what = "sigma", type = "response", newdata = dat %>% dplyr::mutate({{ treat }} :=!!treat_levels[i]), data = dat))
    }

  }

  # ---------------------------------------------------------------------------------------------- #
  # Calculate Average Marginal Effects
  # ---------------------------------------------------------------------------------------------- #

  AME <- list()

  # By time and group ---------------------------------------------------------------------------- #

  if(length(treat_levels)>1) {

    # Relative frequency of each treat_levels by time and group
    treat_props <-
      dat %>%
      dplyr::filter(!!rlang::sym(treat)!=0) %>% # freq. of 0s don't need to be counted
      group_by(across(all_of(c(time, group, treat)))) %>%
      dplyr::summarise(n=n()) %>%
      group_by(across(all_of(c(time, group)))) %>%
      dplyr::summarise(n=n/sum(n)) %>%
      dplyr::mutate(rn=row_number()) %>%
      tidyr::pivot_wider(names_prefix="n.", names_from = rn, values_from = n)

    # Calculate effect by time and group
    AME[[1]] <-
      pred %>%
      group_by(across(all_of(c(time, group)))) %>%
      dplyr::summarise(across(starts_with("pb"), ~ mean(.x-pa))) %>% # AME for each effect (0->treat_levels[1], 0->treat_levels[2], ...)
      ungroup() %>%
      inner_join(treat_props, by=c(time, group)) %>%
      rowwise() %>% dplyr::mutate(effect=sum(c_across(starts_with("n"))*c_across(starts_with("pb")))) %>% # total effect = n.1*pb.1 + n.2*pb.2 + ...
      ungroup() %>%
      dplyr::select(-starts_with("pb"), -starts_with("n"))

  } else {

    AME[[1]] <-
      pred %>%
      group_by(across(all_of(c(time, group)))) %>%
      dplyr::summarise(effect = mean(pb.1-pa)) %>% # AME
      ungroup()

  }

  # By time -------------------------------------------------------------------------------------- #

  if(length(treat_levels)>1) {

    # Relative frequency of each treat_levels by time
    treat_levels <-
      dat %>%
      dplyr::filter(!!rlang::sym(treat)!=0) %>%
      group_by(across(all_of(c(time, treat)))) %>%
      dplyr::summarise(n=n()) %>%
      group_by(across(all_of(c(time)))) %>%
      dplyr::summarise(n=n/sum(n)) %>%
      dplyr::mutate(rn=row_number()) %>%
      tidyr::pivot_wider(names_prefix="n.", names_from = rn, values_from = n)

    # Calculate effect by time
    AME[[2]] <-
      pred %>%
      group_by(across(all_of(c(time)))) %>%
      dplyr::summarise(across(starts_with("pb"), ~ mean(.x-pa))) %>% # AME
      ungroup() %>%
      inner_join(treat_levels, by=c(time)) %>%
      rowwise() %>% dplyr::mutate(effect=sum(c_across(starts_with("n"))*c_across(starts_with("pb")))) %>%
      ungroup() %>%
      dplyr::select(-starts_with("pb"), -starts_with("n"))

  } else {

    AME[[2]] <-
      pred %>%
      group_by(across(all_of(c(time)))) %>%
      dplyr::summarise(effect = mean(pb.1-pa)) %>% # AME
      ungroup()

  }

  return(AME)

}
