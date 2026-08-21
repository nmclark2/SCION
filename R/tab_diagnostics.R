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
          shiny::fluidRow(
            shiny::column(
              5,
              shiny::numericInput(ns("target_fdr"), "Target FDR", value = 0.05, min = 0, max = 1, step = 0.01),
              shiny::actionButton(ns("recompute_fdr"), "Recompute at this FDR", icon = shiny::icon("rotate"))
            ),
            shiny::column(
              5,
              shiny::numericInput(ns("manual_threshold"), "...or apply a flat weight cutoff instead",
                                   value = NA, min = 0, step = 0.01),
              shiny::actionButton(ns("apply_manual_threshold"), "Apply cutoff", icon = shiny::icon("filter"))
            ),
            shiny::column(
              2,
              shiny::br(),
              shiny::actionLink(ns("clear_override"), "Reset to default FDR")
            )
          ),
          shiny::hr(),
          shiny::uiOutput(ns("summary")),
          shiny::downloadButton(ns("download_thresholded"), "Download thresholded network")
        )
      )
    })

    # The permutation results (network_result()$permuted_networks) never need
    # to be recomputed here -- compute_fdr_threshold() is pure rank/p-value math
    # over already-inferred weights, and a flat cutoff is just a data frame
    # filter. Both are effectively free, so "recompute" never re-runs any
    # random forest inference.
    fdr_override <- shiny::reactiveVal(NULL)
    manual_threshold <- shiny::reactiveVal(NULL)

    shiny::observeEvent(network_result(), {
      fdr_override(NULL)
      manual_threshold(NULL)
      res <- network_result()
      if (!is.null(res$fdr_result)) {
        shiny::updateNumericInput(session, "target_fdr", value = res$fdr_result$target_fdr)
      }
    })

    shiny::observeEvent(input$recompute_fdr, {
      res <- network_result()
      shiny::req(res$fdr_result)
      fdr_override(compute_fdr_threshold(res$network, res$permuted_networks, target_fdr = input$target_fdr))
      manual_threshold(NULL)
    })

    shiny::observeEvent(input$apply_manual_threshold, {
      shiny::req(!is.na(input$manual_threshold))
      manual_threshold(input$manual_threshold)
    })

    shiny::observeEvent(input$clear_override, {
      fdr_override(NULL)
      manual_threshold(NULL)
    })

    # single source of truth for the plots/summary/download -- a flat cutoff
    # (if applied) always takes precedence over whichever FDR result (default
    # or recomputed) is currently active, but keeps that FDR result's curve/
    # permuted weights so the FDR plots still show the full picture.
    display_fdr <- shiny::reactive({
      base <- if (!is.null(fdr_override())) fdr_override() else network_result()$fdr_result
      if (is.null(base)) {
        return(NULL)
      }
      mt <- manual_threshold()
      if (!is.null(mt)) {
        base$threshold <- mt
        base$thresholded_network <- network_result()$network[network_result()$network$Weight >= mt, , drop = FALSE]
      }
      base
    })

    output$fdr_curve <- shiny::renderPlot({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      plot_fdr_curve(display_fdr(), type = "curve")
    })

    output$weight_comparison <- shiny::renderPlot({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      plot_fdr_curve(display_fdr(), type = "weight_comparison")
    })

    output$summary <- shiny::renderUI({
      res <- network_result()
      if (is.null(res$fdr_result)) {
        return(NULL)
      }
      fdr <- display_fdr()
      basis <- if (!is.null(manual_threshold())) {
        "manual cutoff"
      } else if (!is.null(fdr_override())) {
        sprintf("FDR < %s", fdr$target_fdr)
      } else {
        sprintf("default FDR < %s", fdr$target_fdr)
      }
      shiny::tagList(
        shiny::strong(sprintf("Threshold: %s (%s)",
                               if (is.na(fdr$threshold)) "not reached" else signif(fdr$threshold, 4), basis)),
        shiny::p(sprintf("%d of %d edges kept", nrow(fdr$thresholded_network), nrow(res$network)))
      )
    })

    output$download_thresholded <- shiny::downloadHandler(
      filename = function() "thresholded_network.tsv",
      content = function(file) {
        fdr <- display_fdr()
        shiny::req(fdr)
        write_scion_network(fdr$thresholded_network, file)
      }
    )

    display_fdr
  })
}
