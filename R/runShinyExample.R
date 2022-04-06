#' @title runShinyExample
#' @description Run to view shiny app demonstrating the ineqx package.
#' @examples runShinyExample()
#' @export runShinyExample
#' @author Benjamin Rosche <benjamin.rosche@@gmail.com>
runShinyExample <- function() {
  appDir <- system.file("ineqx-app", "app.R", package = "ineqx")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
