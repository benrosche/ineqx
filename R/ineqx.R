#' @title Variance function regression + decomposition of the variance and sq. coefficient of variation
#'
#' @description The ineqx package implements Rosche (2021). [...]
#'
#' @param x Character string. Explanatory variable. c.x or i.x. c2.x will specify quadratic polynomials for x.
#' @param y Character string. Dependent variable. Can only be c.x
#' @param xpv Numeric vector of length 2. Choose the values of x that are used for prediction, e.g. c(0,1)
#' @param ystat Character string. Either "Var" or "CV2". Choose to analyze the effect of x on the variance or squared coefficient of variation
#' @param controls Character vector with additional control variables. E.g. c("c.age", "i.sex", ...)
#' @param groupvar Character string. Must be i.x. Grouping variable to decompose variance into within- and between-group components.
#' @param timevar Character string. c.x with specify penalized splines, i.x will specify dummies. Time variable to analyze change over time.
#' @param cf Number. Counterfactual reference point.
#' @param dat Dataframe
#'
#' @return List with six elements: pred_mu, pred_sigma, dW, dB, dD, dT. See details.
#'
#' @examples data(incdat)
#' decomp1 <- ineqx(child, inc, ystat="CV2", groupvar=SES, timevar=year, cf=1, incdat)
#'
#' @export ineqx
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @details
#' The main function of \code{ineqx}. [...]

ineqx <- function(x=NULL, y, xpv=c(0,1), ystat="CV2", weights=NULL, controls=NULL, groupvar=NULL, timevar=NULL, cf=NULL, dat) {

  # dat = incdat %>% na.omit(); x="i.x"; xpv=c(0,1); y="c.inc"; ystat="Var"; controls=NULL; weights=NULL; groupvar="i.group"; timevar="c.year"; cf=0
  # dat = dat.ineqx; x="c.x"; xpv=c(0,1); y="c.fearnings_wk"; ystat="CV2"; groupvar="f_SES"; timevar="c.year"; cf=27; controls=c("i.mother", "i.byear", "c2.age", "c2.famsize", "i.married", "i.race", "i.w_edu")
  # dat = dat1; x="c3.x"; xpv=c(0,1); y="inc"; ystat="Var"; groupvar="i.group"; timevar="i.year"; cf=1; controls=NULL; weights=NULL

  # ---------------------------------------------------------------------------------------------- #
  # Dissect input
  # ---------------------------------------------------------------------------------------------- #

  # X
  if(!is.null(x)) {
    c(x, xcont, xpoly) %<-% dissectVar(x, cicheck="ci")
    form_x <- if(xcont) paste0("I(x^", 1:xpoly, ")", collapse = " + ") else paste0("as.factor(x)")
    form_x <- paste0("(", form_x, ") *")
    nox <- F
  } else {
    form_x <- ""
    nox <- T
  }

  # Y
  y <- dissectVar(y, cicheck = "c")[[1]]

  # Weights
  if(!is.null(weights)) w <- as.symbol(weights) else w <- NULL

  # Controls
  if(!is.null(controls)) {
    ctrl_list <- lapply(controls, FUN=function(x) dissectVar(x))
    form_ctrl <- paste0(lapply(ctrl_list, FUN=function(x) { if(x[[2]]) paste0("I(", x[[1]], "^", 1:x[[3]], ")", collapse = " + ") else paste0("as.factor(", x[[1]], ")") }), collapse = " + ")
    form_ctrl <- paste0(" + ", form_ctrl)
  } else form_ctrl <- ""

  # Groupvar
  groupvar <- dissectVar(groupvar, cicheck = "i")[[1]]

  # Timevar
  time_list <- dissectVar(timevar, cicheck = "ci")
  timevar <- time_list[[1]]
  form_time <- if(time_list[[2]] & time_list[[3]]==1) "time" else if(time_list[[2]] & time_list[[3]]>1) paste0("bs(time, df=", time_list[[3]], ")") else "as.factor(time)"

  # Check variables and values
  if(!nox) if(!as.character(x) %in% names(dat)) stop(paste0(x, " not in dataset"))
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset"))
  if(!as.character(groupvar) %in% names(dat)) stop(paste0(groupvar, " not in dataset"))
  if(!as.character(timevar) %in% names(dat)) stop(paste0(timevar, " not in dataset"))
  if(!is.null(weights)) if(!as.character(weights) %in% names(dat)) stop(paste0(weights, " not in dataset"))
  if(!is.null(controls)) lapply(ctrl_list, FUN=function(x) { if(!as.character(x[[1]]) %in% names(dat)) stop(paste0(x[[1]], " not in dataset")) } )
  if(length(cf) <= 2) if(!cf[1] %in% (dat[timevar] %>% unique() %>% unlist() %>% as.vector())) stop(paste0("cf not observed in ", timevar))

  # Rename variables in dat
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, w = {{ w }}, group={{ groupvar }}, time={{ timevar }})

  # ---------------------------------------------------------------------------------------------- #
  # Variance function regression
  # ---------------------------------------------------------------------------------------------- #

  message("Running variance function regression ...")
  vfr <- gamlss(formula       = as.formula( paste0("y ~ ", form_x, "(as.factor(group) * ", form_time, ")", form_ctrl) ), # x group time x#group x#time group#time x#group#time
                sigma.formula = as.formula( paste0( " ~ ", form_x, "(as.factor(group) * ", form_time, ")", form_ctrl) ),
                weights = w,
                data = dat)

  # ---------------------------------------------------------------------------------------------- #
  # Compute average marginal effects
  # ---------------------------------------------------------------------------------------------- #

  message("Computing average marginal effects ...")
  # 2do: add SE with predict(,se.fit = T)

  AME_mu    <- calcAME(nox, xpv, what="mu", vfr, dat)
  AME_sigma <- calcAME(nox, xpv, what="sigma", vfr, dat)

  # ---------------------------------------------------------------------------------------------- #
  # Decomposition
  # ---------------------------------------------------------------------------------------------- #

  message("Performing decomposition ...")

  # Decompose data at @xpv[1] (pre-treatment) and @xpv[2] (post-treatment)
  if(!is.null(x)) {
    wibe.xpv1 <- wibe("y", "group", "time", dat %>% dplyr::filter(x == xpv[1]), smoothDat = F, rel = F)[[1]] # @xpv[1]
    wibe.xpv2 <- wibe("y", "group", "time", dat %>% dplyr::filter(x == xpv[2]), smoothDat = F, rel = F)[[2]] # @xpv[2]
  } else {
    wibe.xpv1 <- wibe("y", "group", "time", dat, smoothDat = F, rel = F)[[1]]
    wibe.xpv2 <- wibe("y", "group", "time", dat, smoothDat = F, rel = F)[[2]]
  }

  # Gather factual and counterfactual data
  dat.f_cf <- createCF("x", "y", ystat, "group", "time", cf, wibe.xpv1, AME_mu, AME_sigma, dat)

  # Calculate impact ----------------------------------------------------------------------------- #

  # Within
  dW.out <-
    calcImpact(nox, paste0(ystat, "W"), dat.f_cf) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          inner_join(wibe.xpv2 %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
          dplyr::rename(dW=dX) %>%
          ungroup()
      })


  # Between
  dB.out <-
    calcImpact(nox, paste0(ystat, "B"), dat.f_cf) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          inner_join(wibe.xpv2 %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
          dplyr::rename(dB=dX) %>%
          ungroup()
      })

  # Demographic
  dD.out <- dD(nox, dW.out, dB.out)

  # Total
  dT.out <- dT(nox, dW.out, dB.out, ystat, cf)

  # ---------------------------------------------------------------------------------------------- #
  # Rename variables again and return output
  # ---------------------------------------------------------------------------------------------- #

  # Rename

  rnm <- function(input, timevar, groupvar) {
    return(list(input[[1]] %>% dplyr::rename(!!enquo(timevar) := time, !!enquo(groupvar) := group), input[[2]] %>% dplyr::rename(!!enquo(timevar) := time)))
  }

  AME_mu    <- rnm(AME_mu, {{ timevar }}, {{ groupvar }})
  AME_sigma <- rnm(AME_sigma, {{ timevar }}, {{ groupvar }})
  dW.out <- rnm(dW.out, {{ timevar }}, {{ groupvar }})
  dB.out <- rnm(dB.out, {{ timevar }}, {{ groupvar }})
  dD.out <- rnm(dD.out, {{ timevar }}, {{ groupvar }})
  dT.out <- dT.out %>% purrr::map(.f=function(x) x %>% dplyr::rename(!!enquo(timevar) := time))

  # Return

  vars <-  rlang::enexprs(x, y, groupvar, timevar, ystat, cf) %>% as.character()
  names(vars) <- c("x", "y", "groupvar", "timevar", "ystat", "cf")

  out <- list("vars"=vars, "dMu"=AME_mu, "dSigma"=AME_sigma, "dW"=dW.out, "dB"=dB.out, "dD"=dD.out, "dT"=dT.out)
  class(out) <- "ineqx"

  message("Done.")

  return(out)

}
