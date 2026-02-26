#' ineqx ggplot2 theme
#'
#' A custom ggplot2 theme for consistent styling of ineqx plots.
#' Use as \code{+ theme_ineqx} in ggplot2 pipelines.
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_ineqx
#'
#' @export
theme_ineqx <-
  ggplot2::theme_bw() +
  ggplot2::theme(
    text = ggplot2::element_text(size = 15),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
    axis.title.x = ggplot2::element_text(face = "bold"),
    axis.title.y = ggplot2::element_text(face = "bold"),
    axis.line = ggplot2::element_line(color = "black"),
    axis.text = ggplot2::element_text(color = "black", size = 15),
    panel.grid = ggplot2::element_line(linewidth = 0.25),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(color = "#A9A9A9", linetype = 2),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(face = "bold"),
    legend.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.key.width = ggplot2::unit(1.4, "cm"),
    strip.background = ggplot2::element_rect(fill = "white", color = "white"),
    strip.text = ggplot2::element_text(face = "bold", colour = "black",
                                       size = ggplot2::rel(1.2)))
