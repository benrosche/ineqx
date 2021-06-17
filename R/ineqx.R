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

ineqx <- function(x, y, xpv=c(0,1), ystat="CV2", weights=NULL, controls=NULL, groupvar=NULL, timevar=NULL, cf=NULL, dat) {

  # dat = incdat %>% na.omit(); x="c2.child"; xpv=c(0,1); y="c.inc"; ystat="Var"; weights=NULL; controls=NULL; groupvar="i.SES"; timevar="c.year"; cf=1
  # dat = dat.ineqx; x="c.by"; xpv=c(0,1); y="c.fearnings_wk"; ystat="CV2"; groupvar="f_SES"; timevar="c.year"; cf=27; controls=c("i.mother", "i.byear", "c2.age", "c2.famsize", "i.married", "i.race", "i.w_edu")

  # ---------------------------------------------------------------------------------------------- #
  # Prepare input
  # ---------------------------------------------------------------------------------------------- #

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

  # Dissect variables ---------------------------------------------------------------------------- #

  # X
  c(x, xcont, xpoly) %<-% dissectVar(x, cicheck="ci")
  form_x <- if(xcont) paste0("I(x^", 1:xpoly, ")", collapse = " + ") else paste0("as.factor(x)")

  # Y
  y <- dissectVar(y, cicheck = "c")[[1]]

  # Weights
  if(!is.null(weights)) w <- as.symbol(weights) else w <- NULL

  # Controls
  if(!is.null(controls)) {
    ctrl_list <- lapply(controls, FUN=function(x) dissectVar(x))
    form_ctrl <- paste0(lapply(ctrl_list, FUN=function(x) { if(x[[2]]) paste0("I(", x[[1]], "^", 1:x[[3]], ")", collapse = " + ") else paste0("as.factor(", x[[1]], ")") }), collapse = " + ")
  } else form_ctrl <- "0"

  # Groupvar
  groupvar <- dissectVar(groupvar, cicheck = "i")[[1]]

  # Timevar
  time_list <- dissectVar(timevar, cicheck = "ci")
  timevar <- time_list[[1]]
  form_time <- if(time_list[[2]] & time_list[[3]]==1) "time" else if(time_list[[2]] & time_list[[3]]>1) paste0("bs(time, df=", time_list[[3]], ")") else "as.factor(time)"

  # Create formula for mu, and sigma
  form_mu   <- paste0("y ~ (", form_x, ") * (as.factor(group) * ", form_time, ") + ", form_ctrl) # x group time x#group x#time group#time x#group#time
  form_sigma <- paste0(" ~ (", form_x, ") * (as.factor(group) * ", form_time, ") + ", form_ctrl) # "

  # Check variables and values
  if(!as.character(x) %in% names(dat)) stop(paste0(x, " not in dataset"))
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset"))
  if(!as.character(groupvar) %in% names(dat)) stop(paste0(groupvar, " not in dataset"))
  if(!as.character(timevar) %in% names(dat)) stop(paste0(timevar, " not in dataset"))
  if(!is.null(weights)) if(!as.character(weights) %in% names(dat)) stop(paste0(weights, " not in dataset"))
  if(!is.null(controls)) lapply(ctrl_list, FUN=function(x) { if(!as.character(x[[1]]) %in% names(dat)) stop(paste0(x[[1]], " not in dataset")) } )
  if(!cf %in% (dat[timevar] %>% unique() %>% unlist() %>% as.vector())) stop(paste0("cf not observed in ", timevar))

  # Rename variables in dat
  dat <- dat %>% dplyr::rename(x={{ x }}, y={{ y }}, w = {{ w }}, group={{ groupvar }}, time={{ timevar }})

  # ---------------------------------------------------------------------------------------------- #
  # Variance function regression
  # ---------------------------------------------------------------------------------------------- #

  message("Running variance function regression ...")
  vfr <- gamlss(formula       = as.formula(form_mu),
                sigma.formula = as.formula(form_sigma),
                weights = w,
                data = dat)

  # ---------------------------------------------------------------------------------------------- #
  # Compute average marginal effects
  # ---------------------------------------------------------------------------------------------- #

  # Predict at values of xpv[1] and xpv[2]
  newdat.a <- dat %>% dplyr::mutate(x=xpv[1])
  newdat.b <- dat %>% dplyr::mutate(x=xpv[2])

  message("Computing average marginal effects (Mu) ...")
  # 2do: add SE with predict(,se.fit = T)

  pred_mu <-
    dat %>%
    dplyr::select(time, group, x) %>%
    dplyr::mutate(pa=predict(vfr, what = "mu", type = "response", newdata = newdat.a, data = dat)) %>%
    dplyr::mutate(pb=predict(vfr, what = "mu", type = "response", newdata = newdat.b, data = dat))

  AME_mu <- list()

  # Summarize by group and time
  AME_mu[[1]] <-
    pred_mu %>%
    group_by(time, group, x) %>%
    dplyr::summarise(pa=mean(pa), pb=mean(pb), beta=pb-pa) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup() %>%
    dplyr::select(-x, -pa, -pb) %>%
    dplyr::rename(
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

    # Summarize by time
    AME_mu[[2]] <-
      pred_mu %>%
      group_by(time, x) %>%
      dplyr::summarise(pa=mean(pa), pb=mean(pb), beta=pb-pa) %>%
      dplyr::filter(row_number()==1) %>%
      ungroup() %>%
      dplyr::select(-x, -pa, -pb) %>%
      dplyr::rename(
        !!enquo(timevar) := time
      )

  message("Computing average marginal effects (Sigma) ...")
  # 2do: add SE with predict(,se.fit = T)

  pred_sigma <-
    dat %>%
    dplyr::select(time, group, x) %>%
    dplyr::mutate(pa=predict(vfr, what = "sigma", type = "response", newdata = newdat.a, data = dat)) %>%
    dplyr::mutate(pb=predict(vfr, what = "sigma", type = "response", newdata = newdat.b, data = dat))

  AME_sigma <- list()

  # Summarize by group and time
  AME_sigma[[1]] <-
    pred_sigma %>%
    group_by(time, group, x) %>%
    dplyr::summarise(pa=mean(pa), pb=mean(pb), lambda=pb-pa) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup() %>%
    dplyr::select(-x, -pa, -pb) %>%
    dplyr::rename(
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

  # Summarize by time
  AME_sigma[[2]] <-
    pred_sigma %>%
    group_by(time, x) %>%
    dplyr::summarise(pa=mean(pa), pb=mean(pb), lambda=pb-pa) %>%
    dplyr::filter(row_number()==1) %>%
    ungroup() %>%
    dplyr::select(-x, -pa, -pb) %>%
    dplyr::rename(
      !!enquo(timevar) := time
    )

  # Rename variables again

  dat <-
    dat %>%
    dplyr::rename(
      !!enquo(x) := x,
      !!enquo(y) := y,
      !!enquo(groupvar) := group,
      !!enquo(timevar) := time
    )

  if(!is.null(weights)) dat <- dat %>% dplyr::rename(!!enquo(w) := w)

  # ---------------------------------------------------------------------------------------------- #
  # Decomposition
  # ---------------------------------------------------------------------------------------------- #

  message("Performing decomposition ...")

  # Output
  dW.out <- dW({{ x }}, {{ y }}, xpv, ystat, {{ groupvar }}, {{ timevar }}, cf=cf, smoothDat=F, AME_mu, AME_sigma, dat)
  dB.out <- dB({{ x }}, {{ y }}, xpv, ystat, {{ groupvar }}, {{ timevar }}, cf=cf, smoothDat=F, AME_mu, dat)
  dD.out <- dD(dW.out, dB.out)
  dT.out <- dT(dW.out, dB.out, cf)

  # ---------------------------------------------------------------------------------------------- #
  # Return
  # ---------------------------------------------------------------------------------------------- #

  vars <-  rlang::enexprs(x, y, groupvar, timevar, ystat, cf) %>% as.character()
  names(vars) <- c("x", "y", "groupvar", "timevar", "ystat", "cf")

  out <- list("vars"=vars, "dMu"=AME_mu, "dSigma"=AME_sigma, "dW"=dW.out, "dB"=dB.out, "dD"=dD.out, "dT"=dT.out)
  class(out) <- "ineqx"

  message("Done.")

  return(out)

}
