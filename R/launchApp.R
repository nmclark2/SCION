#' Launch the SCION Shiny app
#'
#' A 3-tab app covering the full run -> permute/FDR -> visualize workflow.
#' Requires `shiny`, `shinydashboard`, `shinydashboardPlus`, `shinyjs`, and
#' `plotly`, which are `Suggests` (not hard dependencies) so that scripted/HPC
#' use of the package doesn't require the Shiny stack to be installed.
#'
#' @param ... additional arguments passed to `shiny::shinyApp()`.
#' @return a `shiny.appobj`, as returned by `shiny::shinyApp()`.
#' @export
launchApp <- function(...) {
  required <- c("shiny", "shinydashboard", "shinydashboardPlus", "shinyjs", "plotly")
  missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("launchApp() requires: ", paste(missing_pkgs, collapse = ", "),
         ". Install with install.packages(c(", paste(sprintf("\"%s\"", missing_pkgs), collapse = ", "), ")).")
  }
  shiny::shinyApp(ui = app_UI, server = app_server, ...)
}
