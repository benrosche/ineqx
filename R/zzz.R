# ============================================================================ #
# Misc utility functions
# ============================================================================ #

#' Run Shiny example app
#'
#' Launches an interactive Shiny application demonstrating
#' the variance decomposition.
#'
#' @export
runShinyExample <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the example app. ",
         "Install it with install.packages('shiny')")
  }
  app_dir <- system.file("ineqx-app", package = "ineqx")
  if (app_dir == "") {
    stop("Could not find Shiny app directory. Try reinstalling the package.")
  }
  shiny::runApp(app_dir, display.mode = "normal")
}
