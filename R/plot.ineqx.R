#' @title plot function
#'
#' @description [...]
#'
#' @param ineqx.out ineqx.out object from ineqx()
#' @param type Character string. Plot type. Choose from: dMuP, dMuT, dSigmaP, dSigmaT, dWP, dWT, dBP, dBT, dDP, dDT, dP, dT, dPA
#' @param yint Either 1 or 2. Choose reported effect value. 1: Effect in units of ystat. 2: Effect as change in percent since 'cf' (reference time)
#' @param xlab Character string. x-axis label.
#' @param ylab Vector of length 2 with character strings. y-axis label. c("y1lab", "y2lab")
#' @param llab Vector of length equal to the number of levels of groupvar. Group labels. c("grouplbl1", "grouplbl2", "grouplbl3", ...)
#' @param titl Character string of length 2. Plot titles. c("titl1", "titl2")
#' @param xlim Sequence of numbers.breaks at x-axis. seq(2000,2020, 5)
#' @param hline Number or FALSE. Horizontal line at specified value.
#' @param se Logical. Display standard errors?
#' @param lowess Logical. Display lowess smoother?
#' @param legend Logical. Display legend?
#' @param sav Logical. Save as .png in working directory?
#'
#' @return Returns a ggplot2 object
#'
#' @examples data(incdat)
#' i1 <- ineqx(incdat, by, inc, groupvar = SES, timevar = year)
#' plot(i1, type="dMuSigma2")
#'
#' @export plot.ineqx
#' @export
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @details
#'
#' \bold{Interpretation of the y-axis}
#'
#' The interpretation of \code{yint=1} (default) is in absolute points of CV2\{W/B/T\}.
#'
#' The interpretation of \code{yint=2}
#' \itemize{
#'   \item For \code{plot(ineqx.out, type="dBP")} or \code{plot(ineqx.out, type="dBT")}, it is the change in the proportion of dB/CV2B, dD/CV2B relative to the reference value given in \code{ref}. The interpretation of \code{dB} is analogous.
#'   \item For \code{plot(ineqx.out, type="dT")}, it is the change in the proportion of dW/CV2T, dB/CV2T, dD/CV2T, dT/CV2T relative to the reference value given in ref.
#'   }
#' In the example, the CV2 is used as inequality statistic. But the interpretation for the variance as inequality statistic is analogous.
#'

plot.ineqx <- function(ineqx.out, type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = f1; type="dT"; xlab="x"; ylab=c("y1", "y2"); llab=c("a", "b", "c"); titl=c("T1", "T2"); xlim=seq(1,3,1); hline=0; se=F; lowess=T; legend=T; sav=F

  # Checks
  if(class(ineqx.out)!="ineqx") stop("Must be given an ineqx.out object.")
  if(is.null(ineqx.out$vars["timevar"])) stop("For plotting, ineqx() must include a timevar.")
  if(!stringr::str_detect(type, "dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dCP|dCT|dDP|dDT|dT|dTS|dPA")) stop("type must be dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dCP|dCT|dDP|dDT|dT|dTS|dPA.")
  if(isTRUE(hline)) stop("hline must be a value, e.g. 0.")

  # Call function according to plot type
  if(type %in% c("dMuP", "dMuT", "dSigmaP", "dSigmaT")) {

    out <- plot_dMuSigma(ineqx.out, type=type, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type %in% c("dWP", "dWT", "dBP", "dBT", "dCP", "dCT", "dDP", "dDT")) {

    out <- plot_dWBCD(ineqx.out, type=type, yint=yint, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type %in% c("dT", "dTS")) {

    out <- plot_dT(ineqx.out, type=type, yint=yint, xlab=xlab, ylab=ylab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type=="dPA") {

    out <- plot_dPA(ineqx.out, xlab=xlab, xlim=xlim, sav=sav)

  } else out <- NULL

  return(out)

}

# ================================================================================================ #
# Plot types
# ================================================================================================ #

# plot_dMuSigma

plot_dMuSigma <- function(ineqx.out, type=type, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ggplot theme
  ggpubr_ineqx <-
    theme_pubr() +
    theme(panel.grid.major.x = element_blank() , panel.grid.major.y = element_line(size=.1, color="grey", linetype = "dashed"),
          legend.position="bottom", legend.margin=margin(-20,0,0,0), legend.box.margin=margin(0,0,0,0))

  theme_set(ggpubr_ineqx)

  # Variable names
  x <- ineqx.out$vars[["x"]]
  y <- ineqx.out$vars[["y"]]
  groupvar <- ineqx.out$vars[["groupvar"]]
  timevar <- ineqx.out$vars[["timevar"]]
  ystat <- ineqx.out$vars[["ystat"]]

  # ggplot variables
  if(is.null(titl)) titl <- ""
  if(is.null(xlab)) xlab <- ""
  if(startsWith(type, "dMu")) {
    p1.dat <- ineqx.out$dMu
    if(is.null(ylab)) ylab <- "Mu"
  } else if(startsWith(type, "dSigma")) {
    p1.dat <- ineqx.out$dSigma
    if(is.null(ylab)) ylab <- "Sigma"
  }

  # Plot
  if(endsWith(type, "P")) p1 <- ggplot(data=p1.dat[[1]], aes(x=get(timevar), y=effect, color=as.factor(get(groupvar))))
  if(endsWith(type, "T")) p1 <- ggplot(data=p1.dat[[2]], aes(x=get(timevar), y=effect))
  p1 <-
    p1 +
    geom_line(alpha=1) + geom_point(alpha=1) +
    ggpubr_ineqx +
    labs(title=ifelse(is.null(titl), "", titl[1]),
         color="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), "Mu", ylab[1])) +
    NULL

  # Additions
  if(!is.null(xlim))  p1 <- p1 + scale_x_continuous(breaks = xlim)
  if(!is.null(llab))  p1 <- p1 + scale_colour_discrete(labels = llab)
  if(!isFALSE(hline)) p1 <- p1 + geom_hline(yintercept = hline, col="red")
  if(se==T)           p1 <- p1 + geom_ribbon(aes(ymin=beta-se, ymax=beta+se), linetype = 0, alpha = 0.1)
  if(lowess==T)       p1 <- p1 + geom_smooth(method = "loess", se=F, size = 1.5)

  # Save?
  if(sav==TRUE) {
    ggsave(
      paste0("./", type, ".png"),
      p1,
      width = 6,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}


# plot_dWBCD

plot_dWBCD <- function(ineqx.out, type=type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = f1; type="dCP"; yint=2; xlab=NULL; ylab=""; tlab=c("By economic position", "Total"); llab=c("low", "medium", "high"); titl="Test"; hline=0; legend=T; sav=F; dim=c(11,5)

  # ---------------------------------------------------------------------------------------------- #
  # 2do:
  # - implement se
  # ---------------------------------------------------------------------------------------------- #

  # ggplot theme
  ggpubr_ineqx <-
    theme_pubr() +
    theme(panel.grid.major.x = element_blank() , panel.grid.major.y = element_line(size=.1, color="grey", linetype = "dashed"),
          legend.position="bottom", legend.margin=margin(-20,0,0,0), legend.box.margin=margin(0,0,0,0))

  theme_set(ggpubr_ineqx)

  # Variable names
  x <- ineqx.out$vars[["x"]]
  y <- ineqx.out$vars[["y"]]
  groupvar <- ineqx.out$vars[["groupvar"]]
  timevar <- ineqx.out$vars[["timevar"]]
  ystat <- ineqx.out$vars[["ystat"]]
  ystatvar <- paste0(ystat, substr(type, 2,2))
  typevar <-  rlang::sym(substr(type, 1,2))
  listvar <- if(as.character(typevar) %in% c("dC", "dD")) substitute(dCD) else typevar

  # By group or total?

  if(type %in% c("dWP", "dBP", "dCP", "dDP")) {

    # Data
    impact <- ineqx.out[[listvar]][[1]]

    # Plot
    p1 <-
      ggplot(data=
               impact,
             aes(x=get(timevar), color=factor(get(groupvar)))) +
      ggpubr_ineqx + theme(legend.position="bottom") +
      labs(title=ifelse(is.null(titl), "", titl[1]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), as.character(typevar), ylab[1])) +
      NULL

  } else if(type %in% c("dWT", "dBT", "dCT", "dDT")) {

    # Data
    impact.total <- ineqx.out[[listvar]][[2]]

    # Plot
    p1 <-
      ggplot(data=
               impact.total,
             aes(x=get(timevar))) +
      ggpubr_ineqx + theme(legend.position="bottom") +
      labs(title=ifelse(is.null(titl), "", titl[2]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), as.character(typevar), ylab[2])) +
      NULL

  }

  # y-axis
  if(yint==1) {
    p1 <- p1 + geom_line(aes(y={{ typevar }}))
  } else if(yint==2) {
    p1 <- p1 + geom_line(aes(y= (1 + {{ typevar }} / ( {{ ystatvar }} - {{ typevar }} )) * 100))
  }

  if(!is.null(xlim))  p1 <- p1 + scale_x_continuous(breaks = xlim)
  if(!is.null(llab))  p1 <- p1 + scale_colour_discrete(labels = llab)
  if(!isFALSE(hline)) p1 <- p1 + geom_hline(yintercept = hline, col="red")

  # Save?
  if(sav==TRUE) {

    ggsave(
      paste0("./", type, ".png"),
      p1,
      width = 6,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}

# plot_dT

plot_dT <- function(ineqx.out, type=type, yint=1, xlab=NULL, ylab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ggplot theme
  ggpubr_ineqx <-
    theme_pubr() +
    theme(panel.grid.major.x = element_blank() , panel.grid.major.y = element_line(size=.1, color="grey", linetype = "dashed"),
          legend.position="bottom", legend.margin=margin(-20,0,0,0), legend.box.margin=margin(0,0,0,0))

  theme_set(ggpubr_ineqx)

  # Variable names
  x <- ineqx.out$vars[["x"]]
  y <- ineqx.out$vars[["y"]]
  groupvar <- ineqx.out$vars[["groupvar"]]
  timevar <- ineqx.out$vars[["timevar"]]
  ystat <- ineqx.out$vars[["ystat"]]

  # Plot total ----------------------------------------------------------------------------------- #

  if(type=="dT") {

    # Data
    total  <- ineqx.out[["dT"]][[1]]

    p1 <-
      ggplot(data=
               total %>%
               pivot_longer(c("dW", "dB", "dC", "dD", "dT"), names_to = "delta", values_to = "value"),
             aes(
               x=get(timevar),
               color=factor(delta,
                            levels=c("dW", "dB", "dC", "dD", "dT"),
                            labels = c("within", "between", "compositional", "demographic", "total")),
               size=factor(delta,
                           levels=c("dW", "dB", "dC", "dD", "dT"),
                           labels = c("within", "between", "compositional", "demographic", "total"))
             )) +
      ggpubr_ineqx +
      scale_color_manual(values = c("#00BA38", "#619CFF", "#FF61C3", "#F8766D", "#000000")) +
      scale_size_manual(values = c(0.75,0.75,0.75,0.75,1.2)) +
      labs(title=ifelse(is.null(titl), "", titl[1]),
           color="",
           size="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "dT", ylab[1])) +
      NULL

    if(yint==1) {
      p1 <- p1 + geom_line(aes(y=value))
    } else if(yint==2) {
      ystatvar <- as.symbol(paste0(ystat, substr(type, 2,2)))
      p1 <- p1 + geom_line(aes(y= (1 + value / ( {{ ystatvar }} - value )) * 100))

    }

    if(!is.null(xlim))  p1 <- p1 + scale_x_continuous(breaks = xlim)
    if(!isFALSE(hline)) p1 <- p1 + geom_hline(yintercept = hline, col="red")

  }

  # Plot shares ---------------------------------------------------------------------------------- #

  if(type=="dTS") {

    shares <- ineqx.out[["dT"]][[2]] %>% dplyr::filter(!is.nan(share))

    # Effect of X?
    xnotNULL <- ("dP" %in% names(ineqx.out$dT[[1]]))

    if(xnotNULL) {

      p1 <-
        ggplot(data=shares,
               aes(
                 x=get(timevar),
                 y=share*100,
                 color=factor(d,
                              levels=c("dW", "dB", "dC", "dD", "dP"),
                              labels = c("within", "between", "compositional", "demographic", "pre-x"))
               ))

    } else {

      p1 <-
        ggplot(data=shares,
               aes(
                 x=get(timevar),
                 y=share*100,
                 color=factor(d,
                              levels=c("dW", "dB", "dC", "dD"),
                              labels = c("within", "between", "compositional", "demographic"))
               ))

    }

    p1 <-
      p1 +
      geom_line(size=1) +
      ggpubr_ineqx +
      labs(title=ifelse(is.null(titl), "", titl[2]),
           color="",
           size="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "%", ylab[2])) +
      NULL

    if(!is.null(xlim))  p1 <- p1 + scale_x_continuous(breaks = xlim)
    if(!isFALSE(hline)) p1 <- p1 + geom_hline(yintercept = hline, col="red")
    if(lowess==T)       p1 <- p1 + geom_smooth(method = "loess", se=F, size = 1.5)

  }

  # Save? ---------------------------------------------------------------------------------------- #

  if(sav==T) {
    ggsave(
      "./output/graphs/dT.png",
      p1,
      width = 12,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}

plot_dPA <- function(ineqx.out, xlab, xlim, sav=F) {

  # ggplot theme
  ggpubr_ineqx <-
    theme_pubr() +
    theme(panel.grid.major.x = element_blank() , panel.grid.major.y = element_line(size=.1, color="grey", linetype = "dashed"),
          legend.position="bottom", legend.margin=margin(-20,0,0,0), legend.box.margin=margin(0,0,0,0))

  theme_set(ggpubr_ineqx)

  # Variable names
  ystat   <- ineqx.out$vars[["ystat"]]
  timevar <- ineqx.out$vars[["timevar"]]

  # Effect of X?
  xnotNULL <- ("dP" %in% names(ineqx.out$dT[[1]]))

  # ggplot data
  if(xnotNULL) {
    dat.plt <- ineqx.out$dT[[1]] %>% dplyr::mutate("dT+dP"=dT+dP) %>% pivot_longer(cols=c(paste0("d", ystat, "T"), "dW", "dB", "dC", "dD", "dP", "dT", "dT+dP"), names_to="name", values_to = "value")
  } else {
    dat.plt <- ineqx.out$dT[[1]] %>% pivot_longer(cols=c(paste0("d", ystat, "T"), "dW", "dB", "dC", "dD", "dT"), names_to="name", values_to = "value")
  }

  # Create Plot
  p1 <-
    ggplot(dat.plt, aes(x=get(timevar), y=value, linetype=name, color=name)) +
    ggpubr_ineqx +
    geom_line() +
    labs(title="Predicted change",
         x=ifelse(is.null(xlab), "", xlab),
         y=paste0("Change in ", ystat, "T"),
         color="",
         linetype="") +
    theme(legend.position="bottom")

  # Legend
  if(xnotNULL)  { c("#00BA38", "#619CFF", "#FF61C3", "#F8766D", "#000000")
    p1 <-
      p1 +
      scale_color_manual(values = c("dW" = "#00BA38", "dB" = "#619CFF", "dC" = "#FF61C3", "dD" = "#F8766D", "dT" = "#000000", "dP" = "grey", "dT+dP" = "#D39200", "dVarT" = "black", "dCV2T" = "black")) +
      scale_linetype_manual(values = c("dW" = "solid", "dB" = "solid", "dC" = "solid", "dD" = "solid", "dT" = "solid", "dP" = "solid", "dT+dP" = "solid", "dVarT" = "dashed", "dCV2T" = "dashed"))
  } else {
    p1 <-
      p1 +
      scale_color_manual(values = c("dW" = "#00BA38", "dB" = "#619CFF", "dC" = "#FF61C3", "dD" = "#F8766D", "dT" = "#000000", "dVarT" = "black", "dCV2T" = "black")) +
      scale_linetype_manual(values = c("dW" = "solid", "dB" = "solid", "dC" = "solid", "dD" = "solid", "dT" = "solid", "dVarT" = "dashed", "dCV2T" = "dashed"))
  }

  if(!is.null(xlim))  p1 <- p1 + scale_x_continuous(breaks = xlim)

  # Save?
  if(sav==T) {
    ggsave(
      "./output/graphs/dPA.png",
      p1,
      width = 6,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}
