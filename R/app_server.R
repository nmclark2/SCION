################################################################################
# SERVER
################################################################################

#' @keywords internal
app_server <- function(input, output, session) {
  # shinydashboardPlus boxSidebars finish their client-side layout after Shiny's
  # first renderPlot() measurement pass, which can otherwise briefly report a
  # near-zero plot area ("figure margins too large"). Forcing a resize event once
  # the DOM has settled makes Shiny re-measure before any real plot is requested.
  session$onFlushed(function() {
    shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
  })

  network_result <- runSidebarServer("run", parent_session = session)
  thresholded_network <- diagnosticsTabServer("diagnostics", network_result = network_result)
  visualizeTabServer("visualize", thresholded_network = thresholded_network)
}
