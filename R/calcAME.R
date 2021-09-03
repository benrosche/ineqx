#' @title Calculate AME
#'
#' @description [...]
#'
#' @param x Character string. Treatment variable.
#' @param groupvar Grouping variable
#' @param timevar Time variable
#' @param what "Mu" or "Sigma"
#' @param vfr gamlss output
#' @param dat Dataframe
#'
#' @return List of length 2. Element 1 returns AME by group and time. Elements 2 returns the AME by time.
#'
#' @examples data(incdat)
#' vfr <- gamlss()
#' AME_mu <- calcAME("x", groupvar=SES, timevar=year, what="mu", vfr, incdat)
#'
#' @export calcAME
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

calcAME <- function(x, groupvar, timevar, what, vfr, dat) {

  # x="xt"; groupvar="group"; timevar="time"; what="mu"

  # ---------------------------------------------------------------------------------------------- #
  # Checks
  # ---------------------------------------------------------------------------------------------- #

  if(!as.character(x) %in% names(dat)) stop(paste0(x, " not in dataset"))
  if(!as.character(groupvar) %in% names(dat)) stop(paste0(groupvar, " not in dataset"))
  if(!as.character(timevar) %in% names(dat)) stop(paste0(timevar, " not in dataset"))

  x.vals <- dat[x] %>% unique() %>% unlist() %>% as.vector() %>% sort()
  if(!0 %in% x.vals) stop(paste0("Values of ", x, " must contain 0."))
  x.vals <- x.vals[x.vals!=0]
  if(length(x.vals)>1) warning("Effect of x was calculated as weighted average of ", paste0("0->", x.vals, sep=" "), ".")

  # ---------------------------------------------------------------------------------------------- #
  # Predictions from variance function regression
  # ---------------------------------------------------------------------------------------------- #

  if(what=="mu") {

    # y_hat for each observation if x=0
    pred <-
      dat %>%
      dplyr::select({{ timevar }}, {{ groupvar }} , {{ x }}) %>%
      dplyr::mutate(pa=predict(vfr, what = "mu", type = "response", newdata = dat %>% dplyr::mutate({{ x }} :=0), data = dat))

    # y_hat for each observation and each of the x.vals
    for(i in 1:length(x.vals)) {
      pred <-
        pred %>%
        dplyr::mutate("pb.{i}":=predict(vfr, what = "mu", type = "response", newdata = dat %>% dplyr::mutate({{ x }} :=!!x.vals[i]), data = dat))
    }

  } else if(what=="sigma") {

    # y_hat for each observation if x=0
    pred <-
      dat %>%
      dplyr::select({{ timevar }}, {{ groupvar }} , {{ x }}) %>%
      dplyr::mutate(pa=predict(vfr, what = "sigma", type = "response", newdata = dat %>% dplyr::mutate({{ x }} :=0), data = dat))

    # y_hat for each observation and each of the x.vals
    for(i in 1:length(x.vals)) {
      pred <-
        pred %>%
        dplyr::mutate("pb.{i}":=predict(vfr, what = "sigma", type = "response", newdata = dat %>% dplyr::mutate({{ x }} :=!!x.vals[i]), data = dat))
    }

  }

  # ---------------------------------------------------------------------------------------------- #
  # Calculate Average Marginal Effects
  # ---------------------------------------------------------------------------------------------- #

  AME <- list()

  # By timevar and groupvar ---------------------------------------------------------------------- #

  if(length(x.vals)>1) {

    # Relative frequency of each x.vals by timevar and groupvar
    x.props <-
      dat %>%
      dplyr::filter(!!rlang::sym(x)!=0) %>% # freq. of 0s don't need to be counted
      group_by(across(all_of(c(timevar, groupvar, x)))) %>%
      dplyr::summarise(n=n()) %>%
      group_by(across(all_of(c(timevar, groupvar)))) %>%
      dplyr::summarise(n=n/sum(n)) %>%
      dplyr::mutate(rn=row_number()) %>%
      pivot_wider(names_prefix="n.", names_from = rn, values_from = n)

    # Calculate effect by timevar and groupvar
    AME[[1]] <-
      pred %>%
      group_by(across(all_of(c(timevar, groupvar)))) %>%
      dplyr::summarise(across(starts_with("pb"), ~ mean(.x-pa))) %>% # AME for each effect (0->x.vals[1], 0->x.vals[2], ...)
      ungroup() %>%
      inner_join(x.props, by=c(timevar, groupvar)) %>%
      rowwise() %>% dplyr::mutate(effect=sum(c_across(starts_with("n"))*c_across(starts_with("pb")))) %>% # total effect = n.1*pb.1 + n.2*pb.2 + ...
      ungroup() %>%
      dplyr::select(-starts_with("pb"), -starts_with("n"))

  } else {

    AME[[1]] <-
      pred %>%
      group_by(across(all_of(c(timevar, groupvar)))) %>%
      dplyr::summarise(effect = mean(pb.1-pa)) %>% # AME
      ungroup()

  }

  # By timevar ----------------------------------------------------------------------------------- #

  if(length(x.vals)>1) {

    # Relative frequency of each x.vals by timevar
    x.props <-
      dat %>%
      dplyr::filter(!!rlang::sym(x)!=0) %>%
      group_by(across(all_of(c(timevar, x)))) %>%
      dplyr::summarise(n=n()) %>%
      group_by(across(all_of(c(timevar)))) %>%
      dplyr::summarise(n=n/sum(n)) %>%
      dplyr::mutate(rn=row_number()) %>%
      pivot_wider(names_prefix="n.", names_from = rn, values_from = n)

    # Calculate effect by timevar
    AME[[2]] <-
      pred %>%
      group_by(across(all_of(c(timevar)))) %>%
      dplyr::summarise(across(starts_with("pb"), ~ mean(.x-pa))) %>% # AME
      ungroup() %>%
      inner_join(x.props, by=c(timevar)) %>%
      rowwise() %>% dplyr::mutate(effect=sum(c_across(starts_with("n"))*c_across(starts_with("pb")))) %>%
      ungroup() %>%
      dplyr::select(-starts_with("pb"), -starts_with("n"))

  } else {

    AME[[2]] <-
      pred %>%
      group_by(across(all_of(c(timevar)))) %>%
      dplyr::summarise(effect = mean(pb.1-pa)) %>% # AME
      ungroup()

  }

  return(AME)

}
