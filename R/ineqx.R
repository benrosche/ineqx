#' @title Variance function regression + decomposition of the variance and sq. coefficient of variation
#'
#' @description The ineqx package implements Rosche (2021). [...]
#'
#' @param dat Dataframe
#' @param x Explanatory variable
#' @param y Dependent variable
#' @param ystat "Var" or "CV2"
#' @param controls Character vector with additional control variables
#' @param groupvar To create within/between division
#' @param timevar To analyze change over time
#' @param cf Either FALSE or numeric time
#'
#' @return Returns a list: dW, dB, dD, dT
#'
#' @examples data(incdat)
#' decomp1 <- ineqx(child, inc, ystat="CV2", groupvar=SES, timevar=year, cf=1, incdat)
#'
#' @export ineqx
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

ineqx <- function(x, y, ystat="CV2", controls=NULL, groupvar=NULL, timevar=NULL, timesp=T, cf=1, dat) {

  # dat = incdat; x=substitute(child); y=substitute(inc); ystat="CV2"; groupvar=substitute(SES); timevar=substitute(year); cf=1; controls=NULL
  # dat = dat.ineqx; x=substitute(child); y=substitute(fearnings_wk); ystat="CV2"; groupvar=substitute(f_SES); timevar=substitute(year); timesp=T; cf=1; controls=NULL

  # controls = c("age", "famsize", "married", "race", "w_edu")

  # Rename variables
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, group={{ groupvar }}, time={{ timevar }})

  # Control variables
  if(is.null(controls)) controls <- 0

  # Splines?

  if(timesp==T) {

    form_mu <- "y ~ as.factor(x)*as.factor(group)*pb(time) +"
    form_sigma <- "  ~ as.factor(x)*as.factor(group)*pb(time) +"

  } else {

    form_mu <- "y ~ as.factor(x)*as.factor(group)*as.factor(time) +"
    form_sigma <- "  ~ as.factor(x)*as.factor(group)*as.factor(time) +"

  }

  # Variance function regression
  vfr <- gamlss(formula       = as.formula(paste0(form_mu, paste0(controls, collapse= " + "))),
                sigma.formula = as.formula(paste0(form_sigma, paste0(controls, collapse= " + "))),
                data = dat)

  # Predicted values
  pred_mu <-
    ggpredict(vfr, c("x", "group", "time"), what = "mu", ci.lvl = NA, data=dat) %>%
    as.data.frame() %>%
    dplyr::rename(time=facet) %>%
    dplyr::mutate(across(everything(), ~as.numeric(as.character(.)))) %>% # as.character to prevent 0/1 factor levels to turn to 1/2 values
    pivot_wider(id_cols=c(group, time), names_from = x, names_prefix = "p", values_from = predicted) %>% dplyr::mutate(beta=p1-p0) %>% dplyr::select(-p0, -p1) %>%
    dplyr::rename(
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

  pred_sigma <-
    ggpredict(vfr, c("x", "group", "time"), what = "sigma", ci.lvl = NA, data=dat) %>% # predicts sigma, not (!) sigma2
    as.data.frame() %>%
    dplyr::rename(time=facet) %>%
    dplyr::mutate(across(everything(), ~as.numeric(as.character(.)))) %>%
    pivot_wider(id_cols=c(group, time), names_from = x, names_prefix = "p", values_from = predicted) %>% dplyr::mutate(lambda=p1-p0) %>% dplyr::select(-p0, -p1) %>%
    dplyr::rename(
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

  dat <-
    dat %>%
    dplyr::rename(
      !!enquo(x) := x,
      !!enquo(y) := y,
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

  # Output
  dW.out <- dW(dat, pred_mu, pred_sigma, {{ x }}, {{ y }}, ystat, {{ groupvar }}, {{ timevar }}, cf=cf)
  dB.out <- dB(dat, pred_mu, {{ x }}, {{ y }}, ystat, {{ groupvar }}, {{ timevar }}, cf=cf)
  dD.out <- dD(dW.out, dB.out)
  dT.out <- dT(dW.out, dB.out, ystat=ystat)

  vars <-  rlang::enexprs(x, y, groupvar, timevar, ystat) %>% as.character()
  names(vars) <- c("x", "y", "groupvar", "timevar", "ystat")

  out <- list("vars"=vars, "dMu"=pred_mu, "dSigma"=pred_sigma, "dW"=dW.out, "dB"=dB.out, "dD"=dD.out, "dT"=dT.out)
  class(out) <- "ineqx"

  return(out)

}
