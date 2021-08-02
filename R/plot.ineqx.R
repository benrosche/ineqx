#' @title plot function
#'
#' @description [...]
#'
#' @param ineqx.out ineqx.out object from ineqx()
#' @param type Character string. Plot type. Choose from: dMuP, dMuT, dSigmaP, dSigmaT, dWP, dWT, dBP, dBT, dDP, dDT, dO, dT, dPA
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
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
#'
#' @export plot.ineqx
#' @export

plot.ineqx <- function(ineqx.out, type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = test; type="dT"; xlab="x"; ylab=c("y1", "y2"); llab=c("a", "b", "c"); titl=c("T1", "T2"); xlim=seq(1,3,1); hline=0; se=F; lowess=T; legend=T; sav=F

  # Checks
  if(class(ineqx.out)!="ineqx") stop("Must be given an ineqx.out object.")
  if(!stringr::str_detect(type, "dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dDP|dDT|dT|dPA")) stop("type must be dMuP|dMuT|dSigmaP|dSigmaT|dWP|dWT|dBP|dBT|dDP|dDT|dT|dPA")
  if(isTRUE(hline)) stop("hline must be a value, e.g. 0.")

  # Call function according to plot type
  if(type %in% c("dMuP", "dMuT", "dSigmaP", "dSigmaT")) {

    out <- plot_dMuSigma(ineqx.out, type=type, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type %in% c("dWP", "dWT", "dBP", "dBT", "dDP", "dDT", "dCP", "dCT")) {

    out <- plot_dWBD(ineqx.out, type=type, yint=yint, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type=="dT") {

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
  x <- ineqx.out$vars["x"]
  y <- ineqx.out$vars["y"]
  groupvar <- ineqx.out$vars["groupvar"]
  timevar <- ineqx.out$vars["timevar"]
  ystat <- ineqx.out$vars["ystat"]

  if(type == "dMuP") {

    # Plot dMu by group and time

    p1 <-
      ggplot(data=ineqx.out$dMu[[1]], aes(x=get(timevar), y=beta, color=as.factor(get(groupvar)))) +
      geom_line(alpha=1) +
      geom_point(alpha=1) +
      ggpubr_ineqx +
      labs(title=ifelse(is.null(titl), "", titl[1]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "Mu", ylab[1])) +
      NULL

  } else if (type == "dMuT") {

    # Plot dMu by time

    p1 <-
      ggplot(data=ineqx.out$dMu[[2]], aes(x=get(timevar), y=beta)) +
      geom_line(alpha=1) +
      geom_point(alpha=1) +
      ggpubr_ineqx +
      labs(title=ifelse(is.null(titl), "", titl[1]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "Mu", ylab[1])) +
      NULL

    } else if (type == "dSigmaP") {

    # Plot dSigma2 by group and time

    p1 <-
      ggplot(data=ineqx.out$dSigma[[1]], aes(x=get(timevar), y=lambda, color=as.factor(get(groupvar)))) +
      geom_line(alpha=1) +
      geom_point(alpha=1) +
      ggpubr_ineqx +
      labs(title=ifelse(is.null(titl), "", titl[2]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "Sigma", ylab[2])) +
      NULL

  } else if (type == "dSigmaT") {

    # Plot dSigma2 by time

    p1 <-
      ggplot(data=ineqx.out$dSigma[[2]], aes(x=get(timevar), y=lambda)) +
      geom_line(alpha=1) +
      geom_point(alpha=1) +
      ggpubr_ineqx +
      labs(title=ifelse(is.null(titl), "", titl[2]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), "Sigma", ylab[2])) +
      NULL

  } else p1 <- NULL

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


# plot_dWBD

plot_dWBD <- function(ineqx.out, type=type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = test; type=substitute(dDP); yint=2; timeseq=seq(1985,2020,5); ylab=""; tlab=c("By economic position", "Total"); llab=c("low", "medium", "high"); hline=0; legend=T; sav=F; dim=c(11,5)

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
  x <- ineqx.out$vars["x"]
  y <- ineqx.out$vars["y"]
  groupvar <- ineqx.out$vars["groupvar"]
  timevar <- ineqx.out$vars["timevar"]
  ystat <- ineqx.out$vars["ystat"]
  ystatvar <- as.symbol(paste0(ystat, substr(type, 2,2)))
  typevar <-  as.symbol(substr(type, 1,2))

  # By group or total?

  if(type %in% c("dWP", "dBP", "dDP")) {

    # Data
    impact <- ineqx.out[[typevar]][[1]] %>% dplyr::rename(time=!!as.symbol(timevar), group=!!as.symbol(groupvar))

    # Plot
    p1 <-
      ggplot(data=
               impact,
             aes(x=time, color=factor(group))) +
      ggpubr_ineqx + theme(legend.position="bottom") +
      labs(title=ifelse(is.null(titl), "", titl[1]),
           color="",
           x=ifelse(is.null(xlab), "", xlab),
           y=ifelse(is.null(ylab), as.character(typevar), ylab[1])) +
      NULL

  } else if(type %in% c("dWT", "dBT", "dDT")) {

    # Data
    impact.total <- ineqx.out[[typevar]][[2]] %>% dplyr::rename(time=!!as.symbol(timevar))

    # Plot
    p1 <-
      ggplot(data=
               impact.total,
             aes(x=time)) +
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
    p1 <- p1 + geom_line(aes(y= ( 1 + {{ typevar }} / ( {{ ystatvar }} - {{ typevar }} )) * 100))
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
  x <- ineqx.out$vars["x"]
  y <- ineqx.out$vars["y"]
  groupvar <- ineqx.out$vars["groupvar"]
  timevar <- ineqx.out$vars["timevar"]
  ystat <- ineqx.out$vars["ystat"]

  # Separate total and shares and rename variables
  total  <- ineqx.out[[type]][[1]] %>% dplyr::rename(time=!!as.symbol(timevar))
  shares <- ineqx.out[[type]][[2]] %>% dplyr::rename(time=!!as.symbol(timevar)) %>% dplyr::filter(!is.nan(share))

  # Plot total ----------------------------------------------------------------------------------- #

  p1a <-
    ggplot(data=
             total %>%
             pivot_longer(c("dW", "dB", "dD", "dT"), names_to = "delta", values_to = "value"),
           aes(
             x=time,
             color=factor(delta,
                          levels=c("dW", "dB", "dD", "dT"),
                          labels = c("within", "between", "demographic", "total")),
             size=factor(delta,
                         levels=c("dW", "dB", "dD", "dT"),
                         labels = c("within", "between", "demographic", "total"))
           )) +
    ggpubr_ineqx +
    scale_color_manual(values = c("#F8766D", "#00BA38", "#619CFF", "#000000")) +
    scale_size_manual(values = c(1,1,1,1.3)) +
    labs(title=ifelse(is.null(titl), "", titl[1]),
         color="",
         size="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), "dT", ylab[1])) +
    NULL

  if(yint==1) {
    p1a <- p1a + geom_line(aes(y=value))
  } else if(yint==2) {
    ystatvar <- as.symbol(paste0(ystat, substr(type, 2,2)))
    p1a <- p1a + geom_line(aes(y= ( 1 + value / ( {{ ystatvar }} - value )) * 100))
  }

  if(!is.null(xlim))  p1a <- p1a + scale_x_continuous(breaks = xlim)
  if(!isFALSE(hline)) p1a <- p1a + geom_hline(yintercept = hline, col="red")

  # Plot shares ---------------------------------------------------------------------------------- #

  # Effect of X?
  xnotNULL <- ("dO" %in% names(ineqx.out$dT[[1]]))

  if(xnotNULL) {

    p1b <-
      ggplot(data=shares,
             aes(
               x=time,
               y=share*100,
               color=factor(d,
                            levels=c("dW", "dB", "dD", "dO"),
                            labels = c("within", "between", "demographic", "unrelated to x"))
             ))

  } else {

    p1b <-
      ggplot(data=shares,
             aes(
               x=time,
               y=share*100,
               color=factor(d,
                            levels=c("dW", "dB", "dD"),
                            labels = c("within", "between", "demographic"))
             ))

  }

  p1b <-
    p1b +
    geom_line(size=1) +
    ggpubr_ineqx +
    labs(title=ifelse(is.null(titl), "", titl[2]),
         color="",
         size="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), "%", ylab[2])) +
    NULL

  if(!is.null(xlim))  p1b <- p1b + scale_x_continuous(breaks = xlim)
  if(!isFALSE(hline)) p1b <- p1b + geom_hline(yintercept = hline, col="red")
  if(lowess==T)       p1b <- p1b + geom_smooth(method = "loess", se=F, size = 1.5)

  # Combine plots -------------------------------------------------------------------------------- #

  p1 <- plot_grid(p1a, p1b, align = 'v', nrow = 1, labels = "AUTO")

  # Save?
  if(sav==T) {
    ggsave(
      "./output/graphs/dT.png",
      p1a,
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
  ystat   <- ineqx.out$vars["ystat"]
  timevar <- ineqx.out$vars["timevar"]

  # Effect of X?
  xnotNULL <- ("dO" %in% names(ineqx.out$dT[[1]]))

  # ggplot data
  if(xnotNULL) {
    dat.plt <- ineqx.out$dT[[1]] %>% dplyr::mutate("dT+dO"=dT+dO) %>% pivot_longer(cols=c(paste0("d", ystat, "T"), "dW", "dB", "dD", "dO", "dT", "dT+dO"), names_to="name", values_to = "value")
  } else {
    dat.plt <- ineqx.out$dT[[1]] %>% pivot_longer(cols=c(paste0("d", ystat, "T"), "dW", "dB", "dD", "dT"), names_to="name", values_to = "value")
  }

  # Create Plot
  p1 <-
    ggplot(dat.plt, aes(x=eval(parse(text = timevar)), y=value, linetype=name, color=name)) +
    ggpubr_ineqx +
    geom_line() +
    labs(title="Predicted change",
         x=ifelse(is.null(xlab), "", xlab),
         y=paste0("Change in ", ystat, "T"),
         color="",
         linetype="") +
    theme(legend.position="bottom")

  # Legend
  if(xnotNULL)  {
    p1 <-
      p1 +
      scale_color_manual(values = c("dVarT" = "black", "dCV2T" = "black", "dB" = "blue", "dD" = "yellow", "dO" = "grey", "dT" = "red", "dT+dO" = "orange", "dW" = "green")) +
      scale_linetype_manual(values = c("dVarT" = "dashed", "dCV2T" = "dashed", "dB" = "solid", "dD" = "solid", "dO" = "solid", "dT" = "solid", "dT+dO" = "solid", "dW" = "solid"))
  } else {
    p1 <-
      p1 +
      scale_color_manual(values = c("dVarT" = "black", "dCV2T" = "black", "dB" = "blue", "dD" = "yellow", "dT" = "red", "dW" = "green")) +
      scale_linetype_manual(values = c("dVarT" = "dashed", "dCV2T" = "dashed", "dB" = "solid", "dD" = "solid", "dT" = "solid", "dW" = "solid"))
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
