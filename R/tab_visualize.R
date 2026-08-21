################################################################################
# Module: Visualize
################################################################################

#' @keywords internal
visualizeTabUI <- function(id = "visualize") {
  ns <- shiny::NS(id)
  shinydashboardPlus::box(
    title = "Network", width = 12, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
    shiny::fluidRow(
      shiny::column(6, shiny::checkboxInput(ns("interactive"), "Interactive (visNetwork)", value = TRUE)),
      shiny::column(6, shiny::downloadButton(ns("download"), "Download network file"))
    ),
    shiny::uiOutput(ns("plot_container"))
  )
}

#' Always shows the thresholded network (real network filtered by FDR and/or a manual
#' weight cutoff, whichever is currently active on the Network Diagnostics tab, or the
#' unfiltered real network if neither has been applied there) -- adjust the cutoff on
#' that tab rather than switching views here.
#'
#' @param thresholded_network a reactive holding the currently displayed edge table (or
#'   `NULL`), as returned by the Network Diagnostics tab's server.
#' @noRd
visualizeTabServer <- function(id = "visualize", thresholded_network) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_network <- shiny::reactive({
      net <- thresholded_network()
      shiny::validate(shiny::need(net, "Run a network first (see the sidebar)."))
      net
    })

    output$plot_container <- shiny::renderUI({
      net <- current_network()
      if (nrow(net) == 0) {
        return(shiny::helpText("No edges to display."))
      }
      if (isTRUE(input$interactive) && requireNamespace("visNetwork", quietly = TRUE)) {
        visNetwork::visNetworkOutput(ns("network_vis"), height = "700px")
      } else {
        shiny::plotOutput(ns("network_plot"), height = "700px")
      }
    })

    output$network_vis <- visNetwork::renderVisNetwork({
      plot_network(current_network(), interactive = TRUE)
    })

    output$network_plot <- shiny::renderPlot({
      plot_network(current_network(), interactive = FALSE)
    })

    output$download <- shiny::downloadHandler(
      filename = function() "network.tsv",
      content = function(file) write_scion_network(current_network(), file)
    )
  })
}
