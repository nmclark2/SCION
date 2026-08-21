################################################################################
# Module: Network Diagnostics
#
# Pure display -- permutations are triggered from the sidebar's "Run
# permutations" checkbox (run_scion(permute = TRUE, ...)), not from this tab.
# Weight/out-degree distributions, and the flat weight-cutoff control, show as
# soon as a network exists; the FDR curve/plots and "recompute at this FDR"
# control only appear once permutation results are present -- but a weight
# cutoff can always be applied, permutations or not.
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
    shiny::uiOutput(ns("fdr_plots")),
    shiny::uiOutput(ns("threshold_section"))
  )
}

#' @param network_result a reactiveVal/reactive holding [run_scion()]'s result (or `NULL`),
#'   as returned by the Run sidebar's server.
#' @return a reactive holding the currently displayed edge table (or `NULL` if no network
#'   has been run yet) -- the real network, FDR-thresholded, or manually thresholded,
#'   whichever is currently selected below.
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

    # The permutation results (network_result()$permuted_networks) never need
    # to be recomputed here -- compute_fdr_threshold() is pure rank/p-value math
    # over already-inferred weights, and a flat cutoff is just a data frame
    # filter. Both are effectively free, so "recompute"/"apply" never re-run
    # any random forest inference.
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

    # the active FDR result (default or recomputed), or NULL if permutations were never run
    base_fdr <- shiny::reactive({
      if (!is.null(fdr_override())) fdr_override() else network_result()$fdr_result
    })

    # single source of truth for the summary/download/Visualize tab: a manual
    # cutoff always takes precedence when set; otherwise the active FDR result's
    # thresholded network; otherwise (no permutations run) the full real network.
    display_network <- shiny::reactive({
      res <- network_result()
      if (is.null(res)) {
        return(NULL)
      }
      mt <- manual_threshold()
      if (!is.null(mt)) {
        return(res$network[res$network$Weight >= mt, , drop = FALSE])
      }
      fdr <- base_fdr()
      if (!is.null(fdr)) {
        return(fdr$thresholded_network)
      }
      res$network
    })

    # same as base_fdr(), but with $threshold/$thresholded_network swapped to the
    # manual cutoff when one is applied -- purely so plot_fdr_curve()'s vline
    # reflects whichever cutoff is actually in effect. Only meaningful (and only
    # ever plotted) when permutations were run, since it still needs a curve/
    # permuted_weights to draw.
    display_fdr_for_plot <- shiny::reactive({
      base <- base_fdr()
      if (is.null(base)) {
        return(NULL)
      }
      mt <- manual_threshold()
      if (!is.null(mt)) {
        base$threshold <- mt
        base$thresholded_network <- display_network()
      }
      base
    })

    output$fdr_plots <- shiny::renderUI({
      res <- network_result()
      if (is.null(res) || is.null(res$fdr_result)) {
        return(NULL)
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
        )
      )
    })

    output$fdr_curve <- shiny::renderPlot({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      plot_fdr_curve(display_fdr_for_plot(), type = "curve")
    })

    output$weight_comparison <- shiny::renderPlot({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      plot_fdr_curve(display_fdr_for_plot(), type = "weight_comparison")
    })

    # always available once a network exists -- a flat weight cutoff doesn't
    # require permutations; the "recompute at this FDR" half only shows up
    # when there's an FDR result to recompute.
    output$threshold_section <- shiny::renderUI({
      res <- network_result()
      if (is.null(res)) {
        return(NULL)
      }
      has_fdr <- !is.null(res$fdr_result)
      shiny::fluidRow(
        shinydashboardPlus::box(
          title = "Result", width = 12, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
          shiny::fluidRow(
            if (has_fdr) {
              shiny::column(
                5,
                shiny::numericInput(ns("target_fdr"), "Target FDR", value = 0.05, min = 0, max = 1, step = 0.01),
                shiny::actionButton(ns("recompute_fdr"), "Recompute at this FDR", icon = shiny::icon("rotate"))
              )
            },
            shiny::column(
              5,
              shiny::numericInput(ns("manual_threshold"),
                                   if (has_fdr) "...or apply a flat weight cutoff instead" else "Apply a flat weight cutoff",
                                   value = NA, min = 0, step = 0.01),
              shiny::actionButton(ns("apply_manual_threshold"), "Apply cutoff", icon = shiny::icon("filter"))
            ),
            shiny::column(
              2,
              shiny::br(),
              shiny::actionLink(ns("clear_override"), if (has_fdr) "Reset to default FDR" else "Reset")
            )
          ),
          shiny::hr(),
          shiny::uiOutput(ns("summary")),
          shiny::downloadButton(ns("download_thresholded"), "Download thresholded network")
        )
      )
    })

    output$summary <- shiny::renderUI({
      res <- network_result()
      if (is.null(res)) {
        return(NULL)
      }
      net <- display_network()
      mt <- manual_threshold()
      if (!is.null(mt)) {
        threshold_label <- signif(mt, 4)
        basis <- "manual cutoff"
      } else if (!is.null(res$fdr_result)) {
        fdr <- base_fdr()
        threshold_label <- if (is.na(fdr$threshold)) "not reached" else signif(fdr$threshold, 4)
        basis <- sprintf("%s FDR < %s", if (!is.null(fdr_override())) "recomputed" else "default", fdr$target_fdr)
      } else {
        threshold_label <- 0
        basis <- "no additional cutoff beyond the sidebar's edge weight cutoff"
      }
      shiny::tagList(
        shiny::strong(sprintf("Threshold: %s (%s)", threshold_label, basis)),
        shiny::p(sprintf("%d of %d edges kept", nrow(net), nrow(res$network)))
      )
    })

    output$download_thresholded <- shiny::downloadHandler(
      filename = function() "thresholded_network.tsv",
      content = function(file) {
        net <- display_network()
        shiny::req(net)
        write_scion_network(net, file)
      }
    )

    display_network
  })
}
