#' @title Descriptive and causal variance decompositions
#'
#' @description The ineqx package implements Rosche (202X). [...]
#'
#' @param treat Character string. Treatment variable. Values must be 0/1
#' @param post Character string. Before/after variable. Values must be 0/1
#' @param y Character string. Dependent variable. Variable must be continuous.
#' @param ystat Character string. Either "Var" (default) or "CV2". Choose to analyze (the effect of x on) the variance or squared coefficient of variation.
#' @param group Character string. Variable must be categorical. Grouping variable to decompose variance into within- and between-group components.
#' @param time Character string. c.x with specify penalized splines, i.x will specify dummies. Time variable to analyze change over time.
#' @param weights Character string. Weight variable.
#' @param controls Character vector with additional control variables. E.g. c("c.age", "i.sex", ...)
#' @param decomp Character string. Either "post" (default) or "effect".
#' @param ref Number, vector, or list. Counterfactual reference point. See details.
#' @param AME_mu Dataframe with average marginal effects (Mu)
#' @param AME_sigma Dataframe with average marginal effects (Sigma)
#' @param dat Dataframe
#'
#' @return List with six elements: dMu, dSigma, dW, dB, dCP, dT. See details.
#'
#' @examples data(incdat)
#' decomp1 <- ineqx(...)
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

ineqx <- function(treat=NULL, post=NULL, y, ystat="Var", group=NULL, time=NULL, weights=NULL, controls=NULL, decomp="post", ref=NULL, AME_mu=NULL, AME_sigma=NULL, dat) {

  # dat = cps_sample; treat="mother"; post="byear"; y="earnweekf"; ystat="CV2"; group="SES"; time="i.year10"; ref=1980; decomp="effect"; controls=NULL; weights="earnwtf"; AME_mu=NULL; AME_sigma=NULL

  # ---------------------------------------------------------------------------------------------- #
  # Dissect input ----
  # ---------------------------------------------------------------------------------------------- #

  # treat
  c(treat, treatcont, treatpoly, notreat) %<-% dissectVar(treat, cicheck="i")

  # post
  c(post, postcont, postpoly, nopost) %<-% dissectVar(post, cicheck="i")

  # y
  y <- dissectVar(y, cicheck = "c")[[1]]

  # controls
  if(!is.null(controls)) {
    ctrl_list <- lapply(controls, FUN=function(x) dissectVar(x, cicheck="ci"))
    controls  <- lapply(ctrl_list, FUN=function(x) x[[1]]) %>% as.character()
  }

  # group
  group <- dissectVar(group, cicheck = "i")[[1]]

  # time
  c(time, timecont, timepoly, notime) %<-% dissectVar(time, cicheck = "ci")

  # Are manual values provided as reference?
  refm <- (!is.null(names(ref)) & all(names(ref) %in% c("n", "mu", "sigma", "beta", "lambda")))

  # Check variables and values
  if(notreat & !nopost) stop("post can only be specified if treat is specified as well.")
  if(!notreat) {
    if(!as.character(treat) %in% names(dat)) stop(paste0(treat, " not in dataset."))
    treat_levels <- dat[treat] %>% unique() %>% unlist() %>% as.vector() %>% sort()
    if(!0 %in% treat_levels | length(treat_levels)==1) stop(paste0(treat, " must contain 0 (not treated) and at least one other value (treated)."))
    if(!nopost & any(!treat_levels %in% c(0,1))) stop("If post is specified, treat values must be 0 (not treated) or 1 (treated).")
    if(nopost & length(treat_levels)>2) warning("Treatment effect was calculated as weighted average of ", paste0("0->", treat_levels[treat_levels!=0], sep=" "), " as treat has multiple treatment values.")
  }
  if(!nopost) {
    if(!as.character(post) %in% names(dat)) stop(paste0(post, " not in dataset."))
    post_levels <- dat[post] %>% unique() %>% unlist() %>% as.vector() %>% sort()
    if(!0 %in% post_levels | length(post_levels)==1) stop(paste0(post, " must contain 0 (pre) and at least one other value (post)."))
    if(length(post_levels)>2) warning("Treatment effect was calculated as weighted average of ", paste0("0->", post_levels[post_levels!=0], sep=" "), "as post has multiple post-treatment values.")
  }
  if(!as.character(decomp) %in% c("post", "effect")) stop("decomp must be 'post' or 'effect'.")
  if(!as.character(y) %in% names(dat)) stop(paste0(y, " not in dataset."))
  if(!as.character(group) %in% names(dat)) stop(paste0(group, " not in dataset."))
  if(!is.null(weights)) if(!as.character(weights) %in% names(dat)) stop(paste0(weights, " not in dataset."))
  if(!is.null(controls)) if(any(!controls %in% names(dat))) stop(paste0(paste0(controls[!controls %in% names(dat)], collapse = " + "), " not in dataset."))
  if(!notime) {
    if(!as.character(time) %in% names(dat)) stop(paste0(time, " not in dataset."))
    time_levels <- dat[time] %>% unique() %>% unlist() %>% as.vector() %>% sort() # time_levels
    if(is.null(ref)) stop("Counterfactual reference point must be given.")
    if(!refm) if(!ref[1] %in% time_levels) stop(paste0("ref not observed in ", time, "."))
    if(refm)  if(0 %in% time_levels) stop(paste0(time, "must not contain 0 if manual values are provided for ref."))
  }

  # Rename variables in dat
  dat <-
    dat %>%
    dplyr::mutate(across(everything(), as.vector)) %>% # remove labels and other attributes
    dplyr::select({{ treat }}, {{ post }}, {{ y }}, {{ group }}, {{ time }}, {{ weights }}, all_of(controls)) %>%
    dplyr::rename_with(~paste0("c",1:length(controls)), .cols = all_of(controls)) %>%
    dplyr::rename(treat={{ treat }}, post={{ post }}, y={{ y }}, w = {{ weights }}, group={{ group }}, time={{ time }}) %>%
    drop_na()

  # Add variable: treat*post
  dat <- if(!is.null(post)) dat %>% dplyr::mutate(tp=treat*post) else dat %>% dplyr::mutate(tp=treat)

  # Add variable: time=1 if is.null(time)
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
    form_ctrl <- paste0(lapply(1:length(controls), FUN=function(i) { if(ctrl_list[[i]]$cont & ctrl_list[[i]]$poly==1) paste0("c", i) else if (ctrl_list[[i]]$cont & ctrl_list[[i]]$poly>1) paste0("bs(c", i, ", df=", ctrl_list[[i]]$poly, ")") else paste0("as.factor(c", i, ")") }), collapse = " + ")
    form_ctrl <- paste0(" + ", form_ctrl)
  } else form_ctrl <- ""

  if(!is.null(treat)) form_treat <- paste0(" tp *") else form_treat <- ""

  if(!is.null(post)) form_post <- paste0(" + treat * (as.factor(group)", form_time, ") + post") else form_post <- ""

  # Is AME_mu and AME_sigma provided? ------------------------------------------------------------ #

  if(!notreat & is.null(AME_mu) & is.null(AME_sigma)) {

    # -------------------------------------------------------------------------------------------- #
    # Variance function regression ----
    # -------------------------------------------------------------------------------------------- #

    message("Running variance function regression ...")
    vfr <- gamlss(formula       = as.formula( paste0("y ~", form_treat, " (as.factor(group)", form_time, ")", form_post, form_ctrl) ), # treat post tp group time group#time treat#group treat#time treat#group#time tp#group tp#time tp#group#time
                  sigma.formula = as.formula( paste0( " ~", form_treat, " (as.factor(group)", form_time, ")", form_post, form_ctrl) ),
                  weights = w,
                  data = dat)
    # -------------------------------------------------------------------------------------------- #
    # Compute average marginal effects ----
    # -------------------------------------------------------------------------------------------- #

    message("Computing average marginal effects ...")
    # 2do: add SE with predict(,se.fit = T)

    AME_mu    <- suppressWarnings( calcAME("tp", "group", "time", what="mu", vfr, dat) )
    AME_sigma <- suppressWarnings( calcAME("tp", "group", "time", what="sigma", vfr, dat) )

  } else if(notreat) {

    # If notreat, AME must be empty

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

    # manual AME_mu + AME_sigma

    AME_mu[[1]]    <- AME_mu[[1]] %>% dplyr::rename(time = {{ time }}, group = {{ group }})
    AME_mu[[2]]    <- AME_mu[[2]] %>% dplyr::rename(time = {{ time }})
    AME_sigma[[1]] <- AME_sigma[[1]] %>% dplyr::rename(time = {{ time }}, group = {{ group }})
    AME_sigma[[2]] <- AME_sigma[[2]] %>% dplyr::rename(time = {{ time }})

  }

  # ---------------------------------------------------------------------------------------------- #
  # Decomposition ----
  # ---------------------------------------------------------------------------------------------- #

  message("Performing decomposition ...")

  ## Calculate within/between variance (pre-treatment, post-treatment and their difference = effect) ----

  if(!is.null(treat)) {

    # Filter to differentiate treat=t/post=NULL and treat=t/post=p
    if(nopost) {
      f.pre  <- "treat==0" # untreated
      f.post <- "treat>0"  # treated
    } else {
      f.pre  <- "treat>0&tp==0" # treated, pre-treatment
      f.post <- "treat>0&tp>0"  # treated, post-treatment
    }

    wibe.pre <- wibe(y="y", group="group", time="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f.pre))))
    wibe.post <- wibe(y="y", group="group", time="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f.post))))
    wibe.add  <- wibe.post[[2]] # will be added to tables in dW, dW etc (=wibe.post[[2]])

  } else {

    f.post <- TRUE # no filter

    wibe.pre  <- wibe(y="y", group="group", time="time", weights=w, dat=dat)
    wibe.post <- wibe(y="y", group="group", time="time", weights=w, dat=dat)
    wibe.add  <- wibe.post[[2]]

  }

  # If effect decomposition is desired, the difference of post and pre must be calculated
  if(!notreat & decomp=="effect") {

    if(ystat=="Var") {

      wibe.add <-
        wibe.pre[[1]] %>% dplyr::select(time, group, n, mu, sigma2) %>%
        inner_join(wibe.post[[1]] %>% dplyr::select(time, group, n, mu, sigma2), by=c("group","time"), suffix = c(".pre", ".post")) %>%
        group_by(time) %>%
        dplyr::mutate(pi.post=n.post/sum(n.post)) %>%
        dplyr::summarize(
          VarW=sum(pi.post*(sigma2.post-sigma2.pre)),
          VarB=sum(pi.post*((mu.post-sum(pi.post*mu.post))^2-(mu.pre-sum(pi.post*mu.pre))^2))
          ) %>%
        ungroup() %>%
        dplyr::mutate(VarT=VarW+VarB) %>%
        dplyr::select(time, VarW, VarB, VarT)
      # correct!

    }
    if(ystat=="CV2") {

      wibe.add <-
        wibe.pre[[1]] %>% dplyr::select(time, group, n, mu, sigma2) %>%
        inner_join(wibe.post[[1]] %>% dplyr::select(time, group, n, mu, sigma2), by=c("group","time"), suffix = c(".pre", ".post")) %>%
        group_by(time) %>%
        dplyr::mutate(pi.post=n.post/sum(n.post), pi.pre=n.pre/sum(n.pre)) %>%
        dplyr::summarize(
          CV2W=sum(pi.post*(sigma2.post-sigma2.pre))/sum(pi.post^2*(mu.post^2-mu.pre^2)),
          CV2B=sum(pi.post*((mu.post-sum(pi.post*mu.post))^2-(mu.pre-sum(pi.post*mu.pre))^2))/sum(pi.post^2*(mu.post^2-mu.pre)^2),
        ) %>%
        ungroup() %>%
        dplyr::mutate(CV2T=CV2W+CV2B) %>%
        dplyr::select(time, CV2W, CV2B, CV2T)

      # incorrect! NEEDS TO BE FIXED

    }

  }

  ## Gather information on t0 and t1 ---------------------------------------------------------------

  wibe.pre_postn <-
    wibe.pre[[1]] %>%
    dplyr::select(-n) %>%
    inner_join(
      dat %>% dplyr::filter(eval(parse(text=f.post))) %>% group_by(time, group) %>% dplyr::summarise(n=n()), # post-treatment n
      by=c("time", "group"))

  # > pre-treatment mu + sigma because AME_mu and AME_sigma calculate the treatment effect on pre-treatment mu + beta
  # > post-treatment n because only those who get the treatment will affect inequality

  dat01 <- calc01(group="group", time="time", ref=ref, wibe.pre=wibe.pre_postn, AME_mu=AME_mu, AME_sigma=AME_sigma, notime=notime, dat=dat)

  # If manual references are provided, the implied inequality must be calculated
  if(refm) {

    ref <- 0
    dat0 <- dat01 %>% dplyr::filter(time==0) %>% dplyr::select(time, group, ends_with("1"))

    if(ystat=="Var") {
      wibe.add <-
        wibe.add %>%
        add_row(
          tibble(
            time=0,
            VarW=VarW(dat0$n1, dat0$sigma1+dat0$lambda1),
            VarB=VarB(dat0$n1, dat0$mu1+dat0$beta1),
            VarT=VarW+VarB
          )) %>%
        arrange(time)
    }
    if(ystat=="CV2") {
      wibe.add <-
        wibe.add %>%
        add_row(
          tibble(
            time=0,
            CV2W=CV2W(dat0$n1, dat0$mu1+dat0$beta1, dat0$sigma1+dat0$lambda1),
            CV2B=CV2B(dat0$n1, dat0$mu1+dat0$beta1),
            CV2T=CV2W+CV2B
          )) %>%
        arrange(time)
    }
  }

  ## Calculate impact ------------------------------------------------------------------------------

  # Within
  deltaW <-
    deltaWB(ystat=ystat, wb="w", decomp=decomp, notreat=notreat, dat01=dat01) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.add %>% dplyr::select(time, paste0(ystat, "W"), paste0(ystat, "T")), by=c("time")) %>%
          ungroup()
      })

  # Between
  deltaB <-
    deltaWB(ystat=ystat, wb="b", decomp=decomp, notreat=notreat, dat01=dat01) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.add %>% dplyr::select(time, paste0(ystat, "B"), paste0(ystat, "T")), by=c("time")) %>%
          ungroup()
      })

  # Compositional
  deltaC <- deltaC(notreat, deltaW, deltaB)

  # Pre-treatment
  if(!is.null(treat)) {
    deltaP <- deltaP(notreat, deltaW, deltaB)
  } else {
    deltaP <- NULL
  }

  # Total
  deltaT <- deltaT(notreat, deltaW, deltaB)
  if(ystat=="Var") deltaT[[1]] <- deltaT[[1]] %>% dplyr::mutate(dVarT=VarT-VarT[time == ref[1]])
  if(ystat=="CV2") deltaT[[1]] <- deltaT[[1]] %>% dplyr::mutate(dCV2T=CV2T-CV2T[time == ref[1]])

  # ---------------------------------------------------------------------------------------------- #
  # Rename variables again and return output ----
  # ---------------------------------------------------------------------------------------------- #

  # Function to rename variables within list
  rnm <- function(input, time, group) {

    if(!is.null(input)) {

      if(!is.null(group)) input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(group) := group)
      if(!is.null(time))  input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(time) := time)

      if(!is.null(time))  input[[2]] <- input[[2]] %>% dplyr::rename(!!enquo(time) := time)

    } else {

      input <- NULL

    }

    return(input)

  }

  AME_mu    <- rnm(AME_mu, {{ time }}, {{ group }})
  AME_sigma <- rnm(AME_sigma, {{ time }}, {{ group }})
  deltaW    <- rnm(deltaW, {{ time }}, {{ group }})
  deltaB    <- rnm(deltaB, {{ time }}, {{ group }})
  deltaC    <- rnm(deltaC, {{ time }}, {{ group }})
  deltaP    <- rnm(deltaP, {{ time }}, {{ group }})
  deltaT    <- rnm(deltaT, {{ time }}, NULL)

  # Return

  vars <-  rlang::enexprs(treat, post, y, group, time, ystat, decomp, ref) %>% as.character()
  names(vars) <- c("treat", "post", "y", "group", "time", "ystat", "decomp", "ref")

  out <- list("vars"=vars, "dMu"=AME_mu, "dSigma"=AME_sigma, "dW"=deltaW, "dB"=deltaB, "dC"=deltaC, "dP"=deltaP, "dT"=deltaT)
  class(out) <- "ineqx"

  message("Done.")

  return(out)

}
