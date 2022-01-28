#' @title Descriptive and causal variance decomposition
#'
#' @description The ineqx package implements Rosche (202X). [...]
#'
#' @param treat Character string. Treatment variable. Values must be binary (0/1)
#' @param post Character string. Before/after variable. Values must be binary (0/1)
#' @param y Character string. Dependent variable. Can only be c.x
#' @param ystat Character string. Either "Var" or "CV2". Choose to analyze the effect of x on the variance or squared coefficient of variation
#' @param group Character string. Must be i.x. Grouping variable to decompose variance into within- and between-group components.
#' @param time Character string. c.x with specify penalized splines, i.x will specify dummies. Time variable to analyze change over time.
#' @param decomp Character string. Either "pre", "post", or "effect".
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

  # dat = incdat; treat="i.x"; post="i.t"; y="inc"; ystat="Var"; group="group"; time="i.year"; ref=1; decomp="post"; controls=NULL; weights=NULL; AME_mu=NULL; AME_sigma=NULL

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
    dplyr::select({{ treat }}, {{ post }}, {{ y }}, {{ group }}, {{ time }}, {{ weights }}, all_of(controls)) %>%
    rename_with(~paste0("c",1:length(controls)), .cols = all_of(controls)) %>%
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

  } else if(notreat & is.null(AME_mu) & is.null(AME_sigma)) {

    # If notreat, AME is empty

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

    AME_mu[[1]]    <- AME_mu[[1]] %>% dplyr::rename(time = {{ time }}, group = {{ group }})
    AME_mu[[2]]    <- AME_mu[[2]] %>% dplyr::rename(time = {{ time }})
    AME_sigma[[1]] <- AME_sigma[[1]] %>% dplyr::rename(time = {{ time }}, group = {{ group }})
    AME_sigma[[2]] <- AME_sigma[[2]] %>% dplyr::rename(time = {{ time }})

  }

  # ---------------------------------------------------------------------------------------------- #
  # Decomposition ----
  # ---------------------------------------------------------------------------------------------- #

  message("Performing decomposition ...")

  # Decompose data pre- and post treatment
  if(!is.null(treat)) {

    # Filter expression
    if(nopost) {
      f0 <- "treat==0" # untreated
      f1 <- "treat>0"  # treated
    } else {
      f0 <- "treat>0&tp==0" # treated, pre-treatment
      f1 <- "treat>0&tp>0"  # treated, post-treatment
    }

    n.pt <- dat %>% dplyr::filter(eval(parse(text=f1))) %>% group_by(time, group) %>% dplyr::summarise(n=n()) # post-treatment n

    wibe.xpv0 <-
      wibe(y="y", group="group", time="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f0))))[[1]] %>%
      dplyr::select(-n) %>%
      inner_join(n.pt, by=c("time", "group")) %>%
      dplyr::relocate(n, .after = "group")

    # > Pre-treatment mu + sigma but post-treatment n.
    # > This is because AME_mu and AME_sigma calculate the treatment effect as pre-treatment mu + beta
    # > However, only those who get the treatment, will affect inequality, which is why we need post-treatment n
    # > mutate(n=...) is important for x=x, t=NULL. When t is also specified, the mutation of n does
    #   not make a difference as f0/f1 both lead to the same n.

    # ------ #

    wibe.xpv1 <- wibe(y="y", group="group", time="time", weights=w, dat=dat %>% dplyr::filter(eval(parse(text=f1))))[[2]]

  } else {

    wibe.xpv0 <- wibe(y="y", group="group", time="time", weights=w, dat=dat)[[1]]
    wibe.xpv1 <- wibe(y="y", group="group", time="time", weights=w, dat=dat)[[2]]

  }

  ## Gather information t (=0) and t+1 (=1) --------------------------------------------------------

  dat01 <- calc01(group="group", time="time", ref=ref, wibe.xpv0=wibe.xpv0, AME_mu=AME_mu, AME_sigma=AME_sigma, dat=dat)

  # Add time = 0 if manual references are provided
  if(refm) {

    ref <- 0

    dat0 <- dat01 %>% dplyr::filter(time==0) %>% dplyr::select(time, group, ends_with("1"))

    # Calculate inequality at time 0
    wibe.xpv1 <-
      wibe.xpv1 %>%
      add_row(
        tibble(time=0,
               N=sum(dat0$n1),
               gmu=sum(dat0$n1/sum(dat0$n1)*(dat0$mu1+dat0$beta1)),
               VarW=VarW(dat0$n1, dat0$sigma1+dat0$lambda1),
               VarB=VarB(dat0$n1, dat0$mu1+dat0$beta1),
               CV2W=CV2W(dat0$n1, dat0$mu1+dat0$beta1, dat0$sigma1+dat0$lambda1),
               CV2B=CV2B(dat0$n1, dat0$mu1+dat0$beta1),
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
    dWB(ystat=ystat, wb="w", notreat=notreat, dat01=dat01) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.xpv1 %>% dplyr::select(time, paste0(ystat, c("W","T"))), by=c("time")) %>%
          ungroup()
      })


  # Between
  dB.out <-
    dWB(ystat=ystat, wb="b", notreat=notreat, dat01=dat01) %>%
    purrr::map(
      .f = function(x) {
        x %>%
          left_join(wibe.xpv1 %>% dplyr::select(time, paste0(ystat, c("B","T"))), by=c("time")) %>%
          ungroup()
      })

  # Compositional + Demographic
  dCD.out <- dCD(notreat, dW.out, dB.out)

  # Total
  dT.out <- dT(notreat, dW.out, dB.out)
  dT.out[[1]] <- if(ystat == "CV2") dT.out[[1]] %>% dplyr::mutate(dCV2T=CV2T-CV2T[time == ref[1]]) else dT.out[[1]] %>% dplyr::mutate(dVarT=VarT-VarT[time == ref[1]]) # add actual change in inequality

  # ---------------------------------------------------------------------------------------------- #
  # Rename variables again and return output ----
  # ---------------------------------------------------------------------------------------------- #

  # Rename

  rnm <- function(input, time, group) {

    if(!is.null(group)) input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(group) := group)
    if(!is.null(time))  input[[1]] <- input[[1]] %>% dplyr::rename(!!enquo(time) := time)

    if(!is.null(time))  input[[2]] <- input[[2]] %>% dplyr::rename(!!enquo(time) := time)

    return(input)

  }

  AME_mu    <- rnm(AME_mu, {{ time }}, {{ group }})
  AME_sigma <- rnm(AME_sigma, {{ time }}, {{ group }})
  dW.out    <- rnm(dW.out, {{ time }}, {{ group }})
  dB.out    <- rnm(dB.out, {{ time }}, {{ group }})
  dCD.out   <- rnm(dCD.out, {{ time }}, {{ group }})
  dT.out    <- rnm(dT.out, {{ time }}, NULL)

  # Return

  vars <-  rlang::enexprs(treat, post, y, group, time, ystat, ref) %>% as.character()
  names(vars) <- c("treat", "post", "y", "group", "time", "ystat", "ref")

  out <- list("vars"=vars, "dMu"=AME_mu, "dSigma"=AME_sigma, "dW"=dW.out, "dB"=dB.out, "dCD"=dCD.out, "dT"=dT.out)
  class(out) <- "ineqx"

  message("Done.")

  return(out)

}
