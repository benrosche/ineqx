#' @title plot function
#'
#' @description [...]
#'
#' @param ineqx.out ineqx.out object from ineqx()
#' @param type dMuSigma2|dB|dW|dD|dO|dT
#' @param yint 1 or 2. 1: ..., 2: ...
#' @param xlab
#' @param ylab
#' @param llab
#' @param titl
#' @param xlim
#' @param hline
#' @param se
#' @param lowess
#' @param legend
#' @param sav
#'
#' @return Returns a ggplot2 object
#'
#' @examples data(incdat)
#' i1 <- ineqx(incdat, by, inc, groupvar = SES, timevar = year)
#' plot(i1, type="dMuSigma2")
#'
#' @export plot.ineqx
#' @method plot ineqx
#'
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>

plot.ineqx <- function(ineqx.out, type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = ineqx.out1; type="dT"; xlab="x"; ylab=c("y1", "y2"); llab=c("a", "b", "c"); titl=c("T1", "T2"); xlim=seq(1,3,1); hline=0; se=F; lowess=F; legend=T; sav=F

  # Checks
  if(!all(names(ineqx.out)==c("vars", "dMu", "dSigma", "dW", "dB", "dD", "dT"))) stop("Must be given an ineqx.out object.")
  if(!stringr::str_detect(type, "dMuSigma|dB|dW|dD|dO|dT")) stop("type must be dMuSigma|dB|dW|dD|dO|dT")
  if(isTRUE(hline)) stop("hline must be a value, e.g. 0.")

  if(type=="dMuSigma") {

    out <- plot_dMuSigma(ineqx.out, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type %in% c("dW", "dB", "dD")) {

    out <- plot_dWBD(ineqx.out, type=type, yint=yint, xlab=xlab, ylab=ylab, llab=llab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  } else if(type=="dT") {

    out <- plot_dT(ineqx.out, type=type, yint=yint, xlab=xlab, ylab=ylab, titl=titl, xlim=xlim, hline=hline, se=se, lowess=lowess, legend=legend, sav=sav)

  }

  return(out)

}

# ================================================================================================ #
# Plot types
# ================================================================================================ #

# plot_dMuSigma2

plot_dMuSigma <- function(ineqx.out, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

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

  # Plot dMu
  p1a <-
    ggplot(data=ineqx.out$dMu, aes(x=get(timevar), y=beta, color=as.factor(get(groupvar)))) +
    geom_line(alpha=0.3) +
    geom_point(alpha=1) +
    ggpubr_ineqx + theme(legend.position="none") +
    labs(title=ifelse(is.null(titl), "", titl[1]),
         color="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), "Mu", ylab[1])) +
    NULL

  if(!is.null(xlim))  p1a <- p1a + scale_x_continuous(breaks = xlim)
  if(!is.null(llab))  p1a <- p1a + scale_colour_discrete(labels = llab)
  if(!isFALSE(hline)) p1a <- p1a + geom_hline(yintercept = hline, col="red")
  if(se==T)           p1a <- p1a + geom_ribbon(aes(ymin=beta-se, ymax=beta+se), linetype = 0, alpha = 0.1)
  if(lowess==T)       p1a <- p1a + geom_smooth(method = "loess", se=F, size = 1.5)

  # Plot dSigma2
  p1b <-
    ggplot(data=ineqx.out$dSigma, aes(x=get(timevar), y=lambda, color=as.factor(get(groupvar)))) +
    geom_line(alpha=0.3) +
    geom_point(alpha=1) +
    ggpubr_ineqx + theme(legend.position="none") +
    labs(title=ifelse(is.null(titl), "", titl[2]),
         color="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), "Sigma2", ylab[2])) +
    NULL

  if(!is.null(xlim))  p1b <- p1b + scale_x_continuous(breaks = xlim)
  if(!is.null(llab))  p1b <- p1b + scale_colour_discrete(labels = llab)
  if(!isFALSE(hline)) p1b <- p1b + geom_hline(yintercept = hline, col="red")
  if(se==T)           p1b <- p1b + geom_ribbon(aes(ymin=lambda-se, ymax=lambda+se), linetype = 0, alpha = 0.1)
  if(lowess==T)       p1b <- p1b + geom_smooth(method = "loess", se=F, size = 1.5)

  # Combine Plots
  p1 <- if(legend==F) plot_grid(p1a, p1b, align = 'h', labels="AUTO", nrow = 1) else plot_grid(plot_grid(p1a, p1b, align = 'h', labels="AUTO", nrow = 1), get_legend(p1a + theme(legend.position="bottom")) , ncol=1, rel_heights = c(2, .1))

  # Save?
  if(sav==TRUE) {
    ggsave(
      paste0("./dMuSigma.png"),
      p1,
      width = 12,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}

# plot_dWBD

plot_dWBD <- function(ineqx.out, type=type, yint=1, xlab=NULL, ylab=NULL, llab=NULL, titl=NULL, xlim=NULL, hline=F, se=F, lowess=F, legend=T, sav=F) {

  # ineqx.out = ineqx.out; type=substitute(dB); int=2; timeseq=seq(1985,2020,5); ylab=""; tlab=c("By economic position", "Total"); llab=c("low", "medium", "high"); hline=0; legend=T; sav=F; dim=c(11,5)

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

  # Separate impact and impact.total and rename variables
  impact <- ineqx.out[[type]][[1]] %>% dplyr::rename(time=!!as.symbol(timevar), group=!!as.symbol(groupvar))
  impact.total <- ineqx.out[[type]][[2]] %>% dplyr::rename(time=!!as.symbol(timevar))

  # plot impact

  p1a <-
    ggplot(data=
             impact,
           aes(x=time, color=factor(group))) +
    ggpubr_ineqx + theme(legend.position="none") +
    labs(title=ifelse(is.null(titl), "", titl[1]),
         color="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), type, ylab[1])) +
    NULL

  # Need to be defined separately
  typevar <- as.symbol(type)
  ystatvar <- as.symbol(paste0(ystat, substr(type, 2,2)))

  if(yint==1) {
    p1a <- p1a + geom_line(aes(y={{ typevar }}))
  } else if(int==2) {
    p1a <- p1a + geom_line(aes(y= ( 1 + {{ typevar }} / ( {{ ystatvar }} - {{ typevar }} )) * 100))
  }

  if(!is.null(xlim))  p1a <- p1a + scale_x_continuous(breaks = xlim)
  if(!is.null(llab))  p1a <- p1a + scale_colour_discrete(labels = llab)
  if(!isFALSE(hline)) p1a <- p1a + geom_hline(yintercept = hline, col="red")
  if(lowess==T)       p1a <- p1a + geom_smooth(method = "loess", se=F, size = 1.5)

  # plot impact total

  p1b <-
    ggplot(data=
             impact.total,
           aes(x=time)) +
    ggpubr_ineqx + theme(legend.position="none") +
    labs(title=ifelse(is.null(titl), "", titl[2]),
         color="",
         x=ifelse(is.null(xlab), "", xlab),
         y=ifelse(is.null(ylab), type, ylab[2])) +
    NULL

  if(yint==1) {
    p1b <- p1b + geom_line(aes(y={{ typevar }}))
  } else if(int==2) {
    p1b <- p1b + geom_line(aes(y= ( 1 + {{ typevar }} / ( {{ ystatvar }} - {{ typevar }} )) * 100))
  }

  if(!is.null(xlim))  p1b <- p1b + scale_x_continuous(breaks = xlim)
  if(!isFALSE(hline)) p1b <- p1b + geom_hline(yintercept = hline, col="red")
  if(lowess==T)       p1b <- p1b + geom_smooth(method = "loess", se=F, size = 1.5)

  # Combine plots

  p1 <- if(legend==F) plot_grid(p1a, p1b, align = 'hv', labels="AUTO", nrow = 1) else plot_grid(plot_grid(p1a, p1b, align = 'hv', labels="AUTO", nrow = 1), get_legend(p1a + theme(legend.position="bottom")), ncol=1, rel_heights = c(2, .1))

  # Save?
  if(sav==TRUE) {

    ggsave(
      paste0("./", type, ".png"),
      p1,
      width = 12,
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
  shares <- ineqx.out[[type]][[2]] %>% dplyr::rename(time=!!as.symbol(timevar))

  # plot total
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
  } else if(int==2) {
    ystatvar <- as.symbol(paste0(ystat, substr(type, 2,2)))
    p1a <- p1a + geom_line(aes(y= ( 1 + value / ( {{ ystatvar }} - value )) * 100))
  }

  if(!is.null(xlim))  p1a <- p1a + scale_x_continuous(breaks = xlim)
  if(!isFALSE(hline)) p1a <- p1a + geom_hline(yintercept = hline, col="red")
  if(lowess==T)       p1a <- p1a + geom_smooth(method = "loess", se=F, size = 1.5)

  # plot shares
  p1b <-
    ggplot(data=shares,
           aes(
             x=time,
             y=share*100,
             color=factor(d,
                          levels=c("dW", "dB", "dD", "dO"),
                          labels = c("within", "between", "demographic", "unrelated to x"))
           )) +
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

  # Combine plots
  p1 <- plot_grid(p1a, p1b, align = 'v', nrow = 1, labels = "AUTO")

  # Save?
  if(sav==T) {
    ggsave(
      "./output/graphs/dT-shares.png",
      p1a,
      width = 12,
      height = 4,
      units='in',
      dpi=300
    )
  }

  return(p1)

}
