################################################################################
# Module: Visualize
################################################################################

#' @keywords internal
visualizeTabUI <- function(id = "visualize") {
  ns <- shiny::NS(id)
  shinydashboardPlus::box(
    title = "Network", width = 12, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
    shiny::uiOutput(ns("plot_container")),
    sidebar = shinydashboardPlus::boxSidebar(
      shiny::selectInput(ns("which_network"), "Network to display",
                          choices = c("Real" = "real", "FDR-thresholded" = "thresholded")),
      shiny::checkboxInput(ns("interactive"), "Interactive (visNetwork)", value = TRUE),
      shiny::downloadButton(ns("download"), "Download network file"),
      id = ns("viz_sidebar"), icon = shiny::icon("gears", class = "fa-2xl"), width = 25,
      background = "rgba(91, 98, 104, 0.9)"
    )
  )
}

#' @param network_result a reactiveVal/reactive holding [run_scion()]'s result (or `NULL`).
#' @param fdr_result a reactive holding [compute_fdr_threshold()]'s result (or `NULL`), as
#'   returned by the Permutation tab's server.
#' @noRd
visualizeTabServer <- function(id = "visualize", network_result, fdr_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    current_network <- shiny::reactive({
      if (identical(input$which_network, "thresholded")) {
        fdr <- fdr_result()
        shiny::validate(shiny::need(fdr, "Run permutations first, or switch to the real network."))
        fdr$thresholded_network
      } else {
        res <- network_result()
        shiny::validate(shiny::need(res, "Run a network first (see the sidebar)."))
        res$network
      }
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
      filename = function() paste0(input$which_network, "_network.tsv"),
      content = function(file) write_scion_network(current_network(), file)
    )
  })
}
