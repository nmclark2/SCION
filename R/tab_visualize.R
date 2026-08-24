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
      # a visualization, not the underlying table -- that's already downloadable
      # from the sidebar (full network) and the Network Diagnostics tab
      # (thresholded network); no need for a third, duplicate copy here
      shiny::column(6, shiny::downloadButton(ns("download"), "Download network image (PDF)"))
    ),
    shiny::uiOutput(ns("legend")),
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

    # a separate, plain Shiny UI element rather than something baked into the
    # visNetwork widget itself -- Shiny's renderVisNetwork()/visNetworkOutput()
    # binding doesn't transmit htmlwidgets::prependContent() to the client, so
    # this is the only reliable way to show a legend alongside the interactive
    # view. The static (igraph) view already bakes its own legend into the
    # plotted image, so this only needs to show for the interactive case.
    output$legend <- shiny::renderUI({
      net <- current_network()
      if (nrow(net) == 0 || !isTRUE(input$interactive)) {
        return(NULL)
      }
      network_legend_html()
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
      # legend = FALSE: the legend is rendered separately, above, as its own
      # UI element (output$legend) -- see plot_network()'s "legend" argument.
      plot_network(current_network(), interactive = TRUE, legend = FALSE)
    })

    output$network_plot <- shiny::renderPlot({
      plot_network(current_network(), interactive = FALSE)
    })

    # always the static (igraph) rendering, regardless of the on-screen
    # interactive/static toggle -- a vector PDF is the highest-resolution,
    # most broadly usable format, and rendering the live interactive widget
    # to an image would need extra screenshot tooling (e.g. webshot2) for no
    # real benefit over the static layout.
    output$download <- shiny::downloadHandler(
      filename = function() "network.pdf",
      content = function(file) {
        net <- current_network()
        grDevices::pdf(file, width = 9, height = 7)
        on.exit(grDevices::dev.off(), add = TRUE)
        plot_network(net, interactive = FALSE)
      }
    )
  })
}
