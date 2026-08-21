################################################################################
# Module: Network Diagnostics
#
# Pure display -- permutations are triggered from the sidebar's "Run
# permutations" checkbox (run_scion(permute = TRUE, ...)), not from this tab.
# Weight/out-degree distributions show as soon as a network exists; the FDR
# section only appears once permutation results are present.
################################################################################

#' @keywords internal
diagnosticsTabUI <- function(id = "diagnostics") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::fluidRow(
      shinydashboardPlus::box(
        title = "Edge weight distribution", width = 6, status = "primary", solidHeader = TRUE,
        headerBorder = TRUE,
        shiny::plotOutput(ns("weight_distribution"))
      ),
      shinydashboardPlus::box(
        title = "Regulator out-degree distribution", width = 6, status = "primary", solidHeader = TRUE,
        headerBorder = TRUE,
        shiny::plotOutput(ns("outdegree_distribution"))
      )
    ),
    shiny::uiOutput(ns("fdr_section"))
  )
}

#' @param network_result a reactiveVal/reactive holding [run_scion()]'s result (or `NULL`),
#'   as returned by the Run sidebar's server.
#' @return a reactive holding the current [compute_fdr_threshold()] result (or `NULL`),
#'   i.e. `network_result()$fdr_result`.
#' @noRd
diagnosticsTabServer <- function(id = "diagnostics", network_result) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$weight_distribution <- shiny::renderPlot({
      res <- network_result()
      shiny::validate(shiny::need(res, "Run a network first (see the sidebar)."))
      plot_weight_distribution(res$network)
    })

    output$outdegree_distribution <- shiny::renderPlot({
      res <- network_result()
      shiny::validate(shiny::need(res, "Run a network first (see the sidebar)."))
      plot_outdegree_distribution(res$network)
    })

    output$fdr_section <- shiny::renderUI({
      res <- network_result()
      if (is.null(res) || is.null(res$fdr_result)) {
        return(shiny::fluidRow(
          shinydashboardPlus::box(
            title = "FDR distribution", width = 12, status = "primary", solidHeader = TRUE,
            headerBorder = TRUE,
            shiny::helpText("Check \"Run permutations\" in the sidebar to see the FDR-based",
                             "edge weight cutoff here.")
          )
        ))
      }
      shiny::fluidRow(
        shinydashboardPlus::box(
          title = "FDR curve", width = 6, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
          shiny::plotOutput(ns("fdr_curve"))
        ),
        shinydashboardPlus::box(
          title = "Real vs. permuted weight distribution", width = 6, status = "primary",
          solidHeader = TRUE, headerBorder = TRUE,
          shiny::plotOutput(ns("weight_comparison"))
        ),
        shinydashboardPlus::box(
          title = "Result", width = 12, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
          shiny::uiOutput(ns("summary")),
          shiny::downloadButton(ns("download_thresholded"), "Download thresholded network")
        )
      )
    })

    output$fdr_curve <- shiny::renderPlot({
      fdr <- network_result()$fdr_result
      shiny::validate(shiny::need(fdr, "No permutation results yet."))
      plot_fdr_curve(fdr, type = "curve")
    })

    output$weight_comparison <- shiny::renderPlot({
      fdr <- network_result()$fdr_result
      shiny::validate(shiny::need(fdr, "No permutation results yet."))
      plot_fdr_curve(fdr, type = "weight_comparison")
    })

    output$summary <- shiny::renderUI({
      res <- network_result()
      fdr <- res$fdr_result
      if (is.null(fdr)) {
        return(NULL)
      }
      shiny::tagList(
        shiny::strong(sprintf("Threshold: %s", if (is.na(fdr$threshold)) "not reached" else signif(fdr$threshold, 4))),
        shiny::p(sprintf("%d of %d edges kept", nrow(fdr$thresholded_network), nrow(res$network)))
      )
    })

    output$download_thresholded <- shiny::downloadHandler(
      filename = function() "thresholded_network.tsv",
      content = function(file) {
        fdr <- network_result()$fdr_result
        shiny::req(fdr)
        write_scion_network(fdr$thresholded_network, file)
      }
    )

    shiny::reactive({
      res <- network_result()
      if (!is.null(res)) res$fdr_result else NULL
    })
  })
}
