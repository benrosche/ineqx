#' @title Variance function regression + decomposition of the variance and sq. coefficient of variation
#'
#' @description The ineqx package implements Rosche (2021). [...]
#'
#' @param x Character string. Treatment variable. Values must be binary (0/1)
#' @param t Character string. Before/after variable. Values must be binary (0/1)
#' @param y Character string. Dependent variable. Can only be c.x
#' @param ystat Character string. Either "Var" or "CV2". Choose to analyze the effect of x on the variance or squared coefficient of variation
#' @param groupvar Character string. Must be i.x. Grouping variable to decompose variance into within- and between-group components.
#' @param timevar Character string. c.x with specify penalized splines, i.x will specify dummies. Time variable to analyze change over time.
#' @param ref Number, vector, or list. Counterfactual reference point. See details.
#' @param controls Character vector with additional control variables. E.g. c("c.age", "i.sex", ...)
#' @param weights Character string. Weight variable.
#' @param AME_mu Dataframe with average marginal effects (Mu)
#' @param AME_sigma Dataframe with average marginal effects (Sigma)
#' @param dat Dataframe
#'
#' @return List with six elements: pred_mu, pred_sigma, dW, dB, dCD, dT. See details.
#'
#' @examples data(incdat)
#' decomp1 <- ineqx()
#'
#' @export ineqx
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @details
#' The main function of \code{ineqx}. [...]
#'
#' \bold{Counterfactual reference point (ref)}
#'
#' \code{ref} can either be a number, c(number, "variable name"), or list(n=c(...), mu=c(...), sigma=c(...), beta=c(...), lambda=c(...)). [...]
#'
#' \bold{Return}
#'
#' ...

ineqx <- function(x=NULL, t=NULL, y, ystat="CV2", groupvar=NULL, timevar=NULL, weights=NULL, controls=NULL, ref=NULL, AME_mu=NULL, AME_sigma=NULL, dat) {

  # dat = incdat; x="i.x"; t="i.t"; y="inc"; ystat="CV2"; groupvar="group"; timevar="i.year"; ref=list(beta=c(0,0,0)); controls=NULL; weights=NULL; AME_mu=NULL; AME_sigma=NULL
  # dat = dat.f1; x=NULL; t=NULL; y="earnwk_hh"; ystat="CV2"; groupvar="f_SES"; timevar="i.year2"; ref=1990; controls=NULL; weights="earnwt"; AME_mu=NULL; AME_sigma=NULL

  # ---------------------------------------------------------------------------------------------- #
  # Dissect input ----
  # ---------------------------------------------------------------------------------------------- #

  # x
  c(x, xcont, xpoly, nox) %<-% dissectVar(x, cicheck="i")

  # t
  c(t, tcont, tpoly, not) %<-% dissectVar(t, cicheck="i")

  # y
  y <- dissectVar(y, cicheck = "c")[[1]]

  # controls
  if(!is.null(controls)) {
    ctrl_list <- lapply(controls, FUN=function(x) dissectVar(x, cicheck="ci"))
    controls  <- lapply(ctrl_list, FUN=function(x) x[[1]]) %>% as.character()
  }

  # groupvar
  groupvar <- dissectVar(groupvar, cicheck = "i")[[1]]

  # timevar
  c(timevar, timecont, timepoly, notime) %<-% dissectVar(timevar, cicheck = "ci")

  # Are manual values provided as reference?
  refm <- (!is.null(names(ref)) & all(names(ref) %in% c("n", "mu", "sigma", "beta", "lambda")))

  # Check variables and values
  if(nox & !not) stop("t can only be specified if x is specified as well.")
  if(!nox) {
    if(!as.character(x) %in% names(dat)) stop(paste0(x, " not in dataset."))
    x_levels <- dat[x] %>% unique() %>% unlist() %>% as.vector() %>% sort()
    if(!0 %in% x_levels | length(x_levels)==1) stop(paste0(x, " must contain 0 (not treated) and at least one other value (treated)."))
    if(!not & any(!x_levels %in% c(0,1))) stop("If t is specified, x values must be 0 (not treated) or 1 (treated).")
    if(not & length(x_levels)>2) warning("Treatment effect was calculated as weighted average of ", paste0("0->", x_levels[x_levels!=0], sep=" "), " as x has multiple treatment values.")
  }
  if(!not) {
    if(!as.character(t) %in% names(dat)) stop(paste0(t, " not in dataset."))
    t_levels <- dat[t] %>% unique() %>% unlist() %>% as.vector() %>% sort()
    if(!0 %in% t_levels | length(t_levels)==1) stop(paste0(t, " must contain 0 (pre) and at least one other value (post)."))
    if(length(t_levels)>2) warning("Treatment effect was calculated as weighted average of ", paste0("0->", t_levels[t_levels!=0], sep=" "), "as t has multiple post-treatment values.")
  }
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset."))
  if(!as.character(groupvar) %in% names(dat)) stop(paste0(groupvar, " not in dataset."))
  if(!is.null(weights)) if(!as.character(weights) %in% names(dat)) stop(paste0(weights, " not in dataset."))
  if(!is.null(controls)) if(any(!controls %in% names(dat))) stop(paste0(paste0(controls[!controls %in% names(dat)], collapse = " + "), " not in dataset."))
  if(!notime) {
    if(!as.character(timevar) %in% names(dat)) stop(paste0(timevar, " not in dataset."))
    time_levels <- dat[timevar] %>% unique() %>% unlist() %>% as.vector() %>% sort() # time_levels
    if(is.null(ref)) stop("Counterfactual reference point must be given.")
    if(!refm) if(!ref[1] %in% time_levels) stop(paste0("ref not observed in ", timevar, "."))
    if(refm)  if(0 %in% time_levels) stop(paste0(timevar, "must not contain 0 if manual values are provided for ref."))
  }

  # Rename variables in dat
  dat <-
    dat %>%
    dplyr::select({{ x }}, {{ t }}, {{ y }}, {{ groupvar }}, {{ timevar }}, {{ weights }}, !!controls) %>%
    rename_with(~paste0("c",1:length(controls)), .cols = controls) %>%
    dplyr::rename(x={{ x }}, t={{ t }}, y={{ y }}, w = {{ weights }}, group={{ groupvar }}, time={{ timevar }})

  # Add variable: xt (= treatment * before/after)
  dat <- if(!is.null(t)) dat %>% dplyr::mutate(xt=x*t) else dat %>% dplyr::mutate(xt=x)

  # Add variable: time=1 if is.null(timevar)
  if(notime) dat <- dat %>% dplyr::mutate(time=1)

  # Levels of group and time
  group_levels <- dat %>% .$group %>% unique() %>% sort()
  time_levels  <- dat %>% .$time %>% unique() %>% sort()

  # Weights
  if(!is.null(weights)) w <- "w" else w <- NULL

  # ---------------------------------------------------------------------------------------------- #
  # Create formulas ----
  # ---------------------------------------------------------------------------------------------- #

  if(!notime) {
    form_time <- if(timecont & timepoly==1) " * time" else if(timecont & timepoly>1) paste0(" * bs(time, df=", timepoly, ")") else " * as.factor(time)"
  } else form_time <- ""

  if(!is.null(controls)) {
    form_ctrl <- paste0(lapply(1:length(controls), FUN=function(i) { if(ctrl_list[[i]]$xcont & ctrl_list[[i]]$xpoly==1) paste0("c", i) else if (ctrl_list[[i]]$xcont & ctrl_list[[i]]$xpoly>1) paste0("bs(c", i, ", df=", ctrl_list[[i]]$xpoly, ")") else paste0("as.factor(c", i, ")") }), collapse = " + ")
    form_ctrl <- paste0(" + ", form_ctrl)
  } else form_ctrl <- ""

  if(!is.null(x)) form_x <- paste0(" xt *") else form_x <- ""

  if(!is.null(t)) form_t <- paste0(" + x * (as.factor(group)", form_time, ") + t") else form_t <- ""

  # Is AME_mu and AME_sigma provided? ------------------------------------------------------------ #

  if(!nox & is.null(AME_mu) & is.null(AME_sigma)) {

    # -------------------------------------------------------------------------------------------- #
    # Variance function regression ----
    # -------------------------------------------------------------------------------------------- #

    message("Running variance function regression ...")
    vfr <- gamlss(formula       = as.formula( paste0("y ~", form_x, " (as.factor(group)", form_time, ")", form_t, form_ctrl) ), # x t xt group time x#group x#time x#group#time xt#group xt#time xt#group#time group#time
                  sigma.formula = as.formula( paste0( " ~", form_x, " (as.factor(group)", form_time, ")", form_t, form_ctrl) ),
                  weights = w,
                  data = dat)
    # -------------------------------------------------------------------------------------------- #
    # Compute average marginal effects ----
    # -------------------------------------------------------------------------------------------- #

    message("Computing average marginal effects ...")
    # 2do: add SE with predict(,se.fit = T)

    AME_mu    <- suppressWarnings( calcAME("xt", "group", "time", what="mu", vfr, dat) )
    AME_sigma <- suppressWarnings( calcAME("xt", "group", "time", what="sigma", vfr, dat) )

  } else if(nox & is.null(AME_mu) & is.null(AME_sigma)) {

    # If nox, AME is empty

    AME_mu <- list()
    AME_mu[[1]] <-
      tibble() %>%
      expand(time=time_levels, group=group_levels) %>%
      dplyr::mutate(effect=0)

    AME_mu[[2]] <-
      tibble() %>%
      expand(time=time_levels) %>%
      dplyr::mutate(effect=0)

    AME_sigma <- AME_mu

  } else {

    AME_mu[[1]]    <- AME_mu[[1]] %>% dplyr::rename(time = {{ timevar }}, group = {{ groupvar }})
    AME_mu[[2]]    <- AME_mu[[2]] %>% dplyr::rename(time = {{ timevar }})
    AME_sigma[[1]] <- AME_sigma[[1]] %>% dplyr::rename(time = {{ timevar }}, group = {{ groupvar }})
    AME_sigma[[2]] <- AME_sigma[[2]] %>% dplyr::rename(time = {{ timevar }})

  }

  # ---------------------------------------------------------------------------------------------- #
  # Decomposition ----
  # ---------------------------------------------------------------------------------------------- #

  message("Performing decomposition ...")

  # Decompose data pre- and post treatment
  if(!is.null(x)) {

    # Filter expression
    if(not) {
      f0 <- "x==0" # untreated
      f1 <- "x>0"  # treated
    } else {
      f0 <- "x>0&xt==0" # treated, pre-treatment
      f1 <- "x>0&xt>0"  # treated, post-treatment
    }

    n.pt <- dat %>% dplyr::filter(eval(parse(text=f1))) %>% group_by(time, group) %>% dplyr::summarise(n=n()) # post-treatment n

    wibe.xpv0 <-
      wibe(y="y", groupvar="group", timevar="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f0))))[[1]] %>%
      dplyr::select(-n) %>%
      inner_join(n.pt, by=c("time", "group")) %>%
      dplyr::relocate(n, .after = "group")

    # > Pre-treatment mu + sigma but post-treatment n.
    # > This is because AME_mu and AME_sigma calculate the treatment effect as pre-treatment mu + beta
    # > However, only those who get the treatment, will affect inequality, which is why we need post-treatment n
    # > mutate(n=...) is important for x=x, t=NULL. When t is also specified, the mutation of n does
    #   not make a difference as f0/f1 both lead to the same n.

    # ------ #

    wibe.xpv1 <- wibe(y="y", groupvar="group", timevar="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f1))))[[2]]

  } else {

    wibe.xpv0 <- wibe(y="y", groupvar="group", timevar="time", weights=w, dat=dat)[[1]]
    wibe.xpv1 <- wibe(y="y", groupvar="group", timevar="time", weights=w, dat=dat)[[2]]

  }

  ## Gather factual and counterfactual data --------------------------------------------------------

  dat.f_cf <- createCF("group", "time", ref, wibe.xpv0, AME_mu, AME_sigma, dat)

  # Add time = 0 if manual references are provided
  if(refm) {

    ref <- 0

    dat0 <- dat.f_cf %>% dplyr::filter(time==0) %>% dplyr::select(time, group, ends_with(".f"))

    # Calculate inequality at time 0
    wibe.xpv1 <-
      wibe.xpv1 %>%
      add_row(
        tibble(time=0,
               N=sum(dat0$n.f),
               gmu=sum(dat0$n.f/sum(dat0$n.f)*(dat0$mu.f+dat0$beta.f)),
               VarW=VarW(dat0$n.f, dat0$sigma.f+dat0$lambda.f),
               VarB=VarB(dat0$n.f, dat0$mu.f+dat0$beta.f),
               CV2W=CV2W(dat0$n.f, dat0$mu.f+dat0$beta.f, dat0$sigma.f+dat0$lambda.f),
               CV2B=CV2B(dat0$n.f, dat0$mu.f+dat0$beta.f),
               VarWBRatio=VarW/VarB,
               CV2WBRatio=CV2W/CV2B,
               VarT=VarW+VarB,
               CV2T=CV2W+CV2B)
      ) %>%
      arrange(time)

  }

  ## Calculate impact ------------------------------------------------------------------------------

  # Within
  dW.out <-
    dX(nox, paste0(ystat, "W"), dat.f_cf) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.xpv1 %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
          dplyr::rename(dW=dX) %>%
          ungroup()
      })


  # Between
  dB.out <-
    dX(nox, paste0(ystat, "B"), dat.f_cf) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.xpv1 %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
          dplyr::rename(dB=dX) %>%
          ungroup()
      })

  # Compositional + Demographic
  dCD.out <- dCD(nox, dW.out, dB.out)

  # Total
  dT.out <- dT(nox, dW.out, dB.out, ystat)
  dT.out[[1]] <- if(ystat == "CV2") dT.out[[1]] %>% dplyr::mutate(dCV2T=CV2T-CV2T[time == ref[1]]) else dT.out[[1]] %>% dplyr::mutate(dVarT=VarT-VarT[time == ref[1]]) # add actual change in inequality

  # ---------------------------------------------------------------------------------------------- #
  # Rename variables again and return output ----
  # ---------------------------------------------------------------------------------------------- #

  # Rename

  rnm <- function(input, timevar, groupvar) {

    if(!is.null(groupvar)) input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(groupvar) := group)
    if(!is.null(timevar))  input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(timevar) := time)

    if(!is.null(timevar))  input[[2]] <- input[[2]] %>% dplyr::rename(!!enquo(timevar) := time)

    return(input)

  }

  AME_mu    <- rnm(AME_mu, {{ timevar }}, {{ groupvar }})
  AME_sigma <- rnm(AME_sigma, {{ timevar }}, {{ groupvar }})
  dW.out    <- rnm(dW.out, {{ timevar }}, {{ groupvar }})
  dB.out    <- rnm(dB.out, {{ timevar }}, {{ groupvar }})
  dCD.out   <- rnm(dCD.out, {{ timevar }}, {{ groupvar }})
  dT.out    <- rnm(dT.out, {{ timevar }}, NULL)

  # Return

  vars <-  rlang::enexprs(x, t, y, groupvar, timevar, ystat, ref) %>% as.character()
  names(vars) <- c("x", "t", "y", "groupvar", "timevar", "ystat", "ref")

  out <- list("vars"=vars, "dMu"=AME_mu, "dSigma"=AME_sigma, "dW"=dW.out, "dB"=dB.out, "dCD"=dCD.out, "dT"=dT.out)
  class(out) <- "ineqx"

  message("Done.")

  return(out)

}
