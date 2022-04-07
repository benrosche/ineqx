#' @title theme_ineqx
#' @description ggplot theme used for the plotting. Note that `theme_ineqx` rather than `theme_ineqx()` must be used.
#' @examples theme_ineqx
#' @export theme_ineqx
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
theme_ineqx <-
  theme_bw() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(face="bold", hjust = 0.5),
    axis.title.x = element_text(face="bold"),
    axis.title.y = element_text(face="bold"),
    axis.line = element_line(color = "black"),
    axis.text=element_text(color="black", size = 15),
    panel.grid = element_line(size=0.25),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#A9A9A9", linetype = 2),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    legend.title = element_text(face="bold"),
    legend.text = element_text(face="bold"),
    legend.position="bottom",
    legend.margin=margin(-20,0,0,0),
    legend.box.margin=margin(0,0,0,0),
    legend.key.width=unit(1.4,"cm"),
    strip.background=element_rect(fill="white", color="white"),
    strip.text=element_text(face="bold", colour = "black", size=rel(1.2)))
theme_set(theme_ineqx)
