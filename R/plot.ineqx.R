#' @title plot function
#'
#' @description [...]
#'
#' @param ineqx.out ineqx.out object from ineqx()
#' @param type Character string. Plot type. Choose from: dMuP, dMuT, dSigmaP, dSigmaT, dWP, dWT, dBP, dBT, dPP, dPT, dT, dPA
#' @param yscale Either 1 or 2. Choose value on y-axis. 1: Effect in units of ystat. 2: Effect as \% of value at reference time (ref)
#'
#' @return Returns a ggplot2 object
#'
#' @examples data(incdat)
#' i1 <- ineqx(y=inc, group=SES, time=year, ref=1, dat=incdat)
#' plot(i1, type="dMuSigma2")
#'
#' @export plot.ineqx
#' @export
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @details
#'
#' The y-axis can be scaled in two ways:
#'
#' \code{yscale=1} (default): absolute points of Var or CV2.\cr
#' \code{yscale=2}: \% of the value of Var or CV2 at reference time. \cr
#' (Total Var or CV2 with \code{type="dT"} and within-group & between-group Var or CV2 with \code{type="dW"} & \code{type="dB"}, respectively.)

plot.ineqx <- function(ineqx.out, type, yscale=1) {

  # ---------------------------------------------------------------------------------------------- #
  # Checks
  # ---------------------------------------------------------------------------------------------- #

  if(class(ineqx.out)!="ineqx") stop("Must be given an ineqx.out object.")
  if(is.null(ineqx.out$vars["time"])) stop("For plotting, ineqx() must include a time variable.")
  if(!stringr::str_detect(type, "dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dCP|dCT|dPP|dPT|dT|dTS|dPA")) stop("type must be dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dCP|dCT|dPP|dPT|dT|dTS|dPA.")

  # ---------------------------------------------------------------------------------------------- #
  # Variables
  # ---------------------------------------------------------------------------------------------- #

  treatvar <- ineqx.out$vars[["treat"]]
  groupvar <- ineqx.out$vars[["group"]]
  timevar <- ineqx.out$vars[["time"]]
  ystat <- ineqx.out$vars[["ystat"]]
  ystatvar <- rlang::sym(paste0(ystat, substr(type,2,2)))
  typevar <-  rlang::sym(substr(type, 1,2))

  # ---------------------------------------------------------------------------------------------- #
  # ggplot theme
  # ---------------------------------------------------------------------------------------------- #

  update_geom_defaults("line", list(size = 1))

  # ---------------------------------------------------------------------------------------------- #
  # Plots
  # ---------------------------------------------------------------------------------------------- #

  if(type %in% c("dMuP", "dMuT", "dSigmaP", "dSigmaT")) {

    # -------------------------------------------------------------------------------------------- #
    # dMu, dSigma ----
    # -------------------------------------------------------------------------------------------- #

    # Data
    if(startsWith(type, "dMu")) {
      dat <- ineqx.out$dMu
      ylab <- "Mu"
    }
    if(startsWith(type, "dSigma")) {
      dat <- ineqx.out$dSigma
      ylab <- "Sigma"
    }

    # Plot
    if(endsWith(type, "P")) {
      out <-
        ggplot(data=dat[[1]], aes(x=get(timevar), y=effect, color=as.factor(get(groupvar)), linetype=as.factor(get(groupvar)))) +
        geom_line() +
        geom_point() +
        theme_ineqx +
        scale_linetype_manual(values=c(1,3,2,4,6)) +
        labs(x="", y=ylab, color=as.character(groupvar), linetype=as.character(groupvar)) +
        NULL
    }
    if(endsWith(type, "T")) {
      out <-
        ggplot(data=dat[[2]], aes(x=get(timevar), y=effect)) +
        geom_line() +
        geom_point() +
        theme_ineqx +
        scale_linetype_manual(values=c(1,3,2,4,6)) +
        labs(x="", y=ylab) +
        NULL
    }

  } else if(type %in% c("dWP", "dWT", "dBP", "dBT", "dCP", "dCT", "dPP", "dPT")) {

    # -------------------------------------------------------------------------------------------- #
    # dW, dB, dC, dP ----
    # -------------------------------------------------------------------------------------------- #

    if(type %in% c("dWP", "dBP", "dCP", "dPP")) {

      # Data
      dat <- ineqx.out[[typevar]][[1]]
      if(yscale==1) dat <- dat %>% dplyr::mutate(y={{ typevar }})
      if(yscale==2) dat <- dat %>% dplyr::mutate(y=( 1 + {{ typevar }}/({{ ystatvar }}-{{ typevar }}) ) * 100)

      # Plot
      out <-
        ggplot(data=dat, aes(x=get(timevar), y=y, color=factor(get(groupvar)), linetype=factor(get(groupvar)))) +
        geom_line() +
        geom_point() +
        theme_ineqx +
        scale_linetype_manual(values=c(1,3,2,4,6)) +
        labs(x="", y=as.character(typevar), color=as.character(groupvar), linetype=as.character(groupvar))

    } else if(type %in% c("dWT", "dBT", "dCT", "dPT")) {

      # Data
      dat <- ineqx.out[[typevar]][[2]]
      if(yscale==1) dat <- dat %>% dplyr::mutate(y={{ typevar }})
      if(yscale==2) dat <- dat %>% dplyr::mutate(y=( 1 + {{ typevar }}/({{ ystatvar }}-{{ typevar }}) ) * 100)

      # Plot
      out <-
        ggplot(data=dat, aes(x=get(timevar), y=y)) +
        geom_line() +
        geom_point() +
        theme_ineqx +
        scale_linetype_manual(values=c(1,3,2,4,6)) +
        labs(x="", y=as.character(typevar))

    }

  } else if(type %in% c("dT")) {

    # -------------------------------------------------------------------------------------------- #
    # dT ----
    # -------------------------------------------------------------------------------------------- #

    # Data
    if(treatvar=="NULL") {
      dat <-
        ineqx.out[["dT"]][[1]] %>%
        tidyr::pivot_longer(c("dW", "dB", "dC", "dT"), names_to = "delta", values_to = "value") %>%
        dplyr::mutate(delta=factor(delta, levels=c("dT", "dW", "dB", "dC"), labels = c("total", "within", "between", "compositional"), ordered = T))
    }
    if(treatvar!="NULL") {
      dat <-
        ineqx.out[["dT"]][[1]] %>%
        tidyr::pivot_longer(c("dW", "dB", "dC", "dP", "dT"), names_to = "delta", values_to = "value") %>%
        dplyr::mutate(delta=factor(delta, levels=c("dT", "dW", "dB", "dC", "dP"), labels = c("total", "within", "between", "compositional", "pre-treatment"), ordered = T))
    }

    # y-axis scaling
    if(yscale==1) dat <- dat %>% dplyr::mutate(y=value)
    if(yscale==2) dat <- dat %>% dplyr::mutate(y=( 1 + value/(.data[[paste0(ystat,"T")]]-value ) ) * 100)

    # Plot
    out <-
      ggplot(
        data=dat,
        aes(
          x=get(timevar),
          y=y,
          color=delta,
          linetype=delta
        )) +
      geom_line() +
      geom_point() +
      theme_ineqx +
      scale_linetype_manual(values=c(1,3,2,4,6)) +
      labs(x="", y="", color="", linetype="") +
      NULL

  } else if(type=="dTS") {

    # -------------------------------------------------------------------------------------------- #
    # dTS ----
    # -------------------------------------------------------------------------------------------- #

    # Data
    if(treatvar=="NULL") {
      dat <-
        ineqx.out[["dT"]][[2]] %>%
        dplyr::filter(!is.nan(share)) %>%
        dplyr::mutate(delta=factor(delta, levels=c("dW", "dB", "dC"), labels = c("within", "between", "compositional"), ordered = T))
    }
    if(treatvar!="NULL")  {
      dat <-
        ineqx.out[["dT"]][[2]] %>%
        dplyr::filter(!is.nan(share)) %>%
        dplyr::mutate(delta=factor(delta, levels=c("dW", "dB", "dC", "dP"), labels = c("within", "between", "compositional", "pre-treatment"), ordered = T))
    }

    # Plot
    out <-
      ggplot(
        data=dat,
        aes(
          x=get(timevar),
          y=share*100,
          color=delta,
          linetype=delta
        )) +
      geom_line() +
      geom_point() +
      theme_ineqx +
      scale_linetype_manual(values=c(1,3,2,4,6)) +
      labs(x="", y="% of total effect", color="", linetype="") +
      NULL

  } else if(type=="dPA") {

    # -------------------------------------------------------------------------------------------- #
    # dPA ----
    # -------------------------------------------------------------------------------------------- #

    # Data
    if(treatvar!="NULL") {
      delta <- c("dT", "dW", "dB", "dC", "dP", paste0("d", ystat, "T"))
    } else {
      delta <- c("dT", "dW", "dB", "dC", paste0("d", ystat, "T"))
    }

    dat <-
      ineqx.out$dT[[1]] %>%
      tidyr::pivot_longer(cols=all_of(delta), names_to="delta", values_to = "value") %>%
      dplyr::mutate(delta=factor(delta, levels = !!delta, ordered = T))

    # Plot
    out <-
      ggplot(dat, aes(x=get(timevar), y=value, color=delta, linetype=delta)) +
      geom_line() +
      theme_ineqx +
      scale_linetype_manual(values=c(1,3,2,4,6,5)) +
      labs(x="", y=paste0("Change in ", ystat, "T"), color="", linetype="")

  } else {

    out <- NULL

  }

  return(out)

}
