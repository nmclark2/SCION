################################################################################
# Module: Network Diagnostics
#
# Pure display -- permutations are triggered from the sidebar's "Run
# permutations" checkbox (run_scion(permute = TRUE, ...)), not from this tab.
# Weight/out-degree distributions, and the flat weight-cutoff control, show as
# soon as a network exists; the "Target FDR" control and FDR curve/plots only
# appear once permutation results are present -- but a weight cutoff can
# always be applied, permutations or not.
################################################################################

#' @keywords internal
diagnosticsTabUI <- function(id = "diagnostics") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("threshold_section")),
    shiny::fluidRow(
      shinydashboardPlus::box(
        title = "Edge weight distribution", width = 6, status = "primary", solidHeader = TRUE,
        headerBorder = TRUE,
        plotly::plotlyOutput(ns("weight_distribution"))
      ),
      shinydashboardPlus::box(
        title = "Regulator out-degree distribution", width = 6, status = "primary", solidHeader = TRUE,
        headerBorder = TRUE,
        plotly::plotlyOutput(ns("outdegree_distribution"))
      )
    ),
    shiny::uiOutput(ns("fdr_plots"))
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

    output$weight_distribution <- plotly::renderPlotly({
      res <- network_result()
      shiny::validate(shiny::need(res, "Run a network first (see the sidebar)."))
      # always the full, unthresholded network -- the cutoff line only means
      # something if you can see where it falls in the whole distribution
      plotly::ggplotly(plot_weight_distribution(res$network, cutoff = current_threshold()))
    })

    output$outdegree_distribution <- plotly::renderPlotly({
      res <- network_result()
      shiny::validate(shiny::need(res, "Run a network first (see the sidebar)."))
      plotly::ggplotly(plot_outdegree_distribution(res$network))
    })

    # compute_fdr_threshold() is pure rank/p-value math over already-inferred
    # weights (no random forest re-inference), so there's no need to cache a
    # separate "recomputed" result behind a button -- target_fdr is the single
    # source of truth for the FDR criterion, and every plot/summary/download
    # just recomputes from its current value. This also means there's no
    # "default vs. recomputed" distinction to explain to the user: whatever
    # target_fdr currently reads is the FDR in effect, full stop.
    manual_threshold <- shiny::reactiveVal(NULL)

    shiny::observeEvent(network_result(), {
      manual_threshold(NULL)
      res <- network_result()
      if (!is.null(res$fdr_result)) {
        shiny::updateNumericInput(session, "target_fdr", value = res$fdr_result$target_fdr)
      }
    })

    shiny::observeEvent(input$apply_manual_threshold, {
      shiny::req(!is.na(input$manual_threshold))
      manual_threshold(input$manual_threshold)
    })

    # discards any manual cutoff *and* any target_fdr edit, returning to
    # exactly the FDR criterion the original run used -- this is a real reset
    # (both inputs visibly revert), not just clearing internal state the user
    # can't see.
    shiny::observeEvent(input$clear_override, {
      manual_threshold(NULL)
      shiny::updateNumericInput(session, "manual_threshold", value = NA)
      res <- network_result()
      if (!is.null(res$fdr_result)) {
        shiny::updateNumericInput(session, "target_fdr", value = res$fdr_result$target_fdr)
      }
    })

    # the FDR result for whatever target_fdr currently reads, or NULL if
    # permutations were never run. Recomputed live on every target_fdr change
    # -- while the field is mid-edit (blank/invalid), falls back to the
    # original run's result rather than erroring or going blank.
    base_fdr <- shiny::reactive({
      res <- network_result()
      if (is.null(res) || is.null(res$fdr_result)) {
        return(NULL)
      }
      target_fdr <- input$target_fdr
      if (is.null(target_fdr) || is.na(target_fdr)) {
        return(res$fdr_result)
      }
      compute_fdr_threshold(res$network, res$permuted_networks, target_fdr = target_fdr)
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

    # the single cutoff value currently in effect (manual, or FDR-based), for
    # marking on the edge weight distribution -- NULL if neither is active.
    current_threshold <- shiny::reactive({
      mt <- manual_threshold()
      if (!is.null(mt)) {
        return(mt)
      }
      fdr <- base_fdr()
      if (!is.null(fdr) && !is.na(fdr$threshold)) {
        return(fdr$threshold)
      }
      NULL
    })

    output$fdr_plots <- shiny::renderUI({
      res <- network_result()
      if (is.null(res) || is.null(res$fdr_result)) {
        return(NULL)
      }
      shiny::fluidRow(
        shinydashboardPlus::box(
          title = "FDR curve", width = 6, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
          plotly::plotlyOutput(ns("fdr_curve"))
        ),
        shinydashboardPlus::box(
          title = "Real vs. permuted weight distribution", width = 6, status = "primary",
          solidHeader = TRUE, headerBorder = TRUE,
          plotly::plotlyOutput(ns("weight_comparison"))
        )
      )
    })

    output$fdr_curve <- plotly::renderPlotly({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      # threshold/show_target_fdr_line reflect whichever cutoff is actually in
      # effect -- when a manual cutoff overrides the FDR criterion, the vline
      # moves to it and the FDR-target hline disappears, since the FDR target
      # is no longer what determined the cutoff.
      plotly::ggplotly(plot_fdr_curve(base_fdr(), type = "curve", threshold = current_threshold(),
                                       show_target_fdr_line = is.null(manual_threshold())))
    })

    output$weight_comparison <- plotly::renderPlotly({
      shiny::validate(shiny::need(network_result()$fdr_result, "No permutation results yet."))
      plotly::ggplotly(plot_fdr_curve(base_fdr(), type = "weight_comparison"))
    })

    # always available once a network exists -- a flat weight cutoff doesn't
    # require permutations; the "Target FDR" half only shows up when there's
    # an FDR result to compute one from, and updates live as it's edited (no
    # separate "recompute" action -- see base_fdr() above).
    output$threshold_section <- shiny::renderUI({
      res <- network_result()
      if (is.null(res)) {
        return(NULL)
      }
      has_fdr <- !is.null(res$fdr_result)
      shiny::fluidRow(
        shinydashboardPlus::box(
          title = "Threshold network", width = 12, status = "primary", solidHeader = TRUE, headerBorder = TRUE,
          shiny::fluidRow(
            if (has_fdr) {
              shiny::column(
                5,
                shiny::numericInput(ns("target_fdr"), "Target FDR", value = 0.05, min = 0, max = 1, step = 0.01)
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
              shiny::actionLink(ns("clear_override"), "Reset")
            )
          ),
          shiny::hr(),
          shiny::uiOutput(ns("summary")),
          shiny::downloadButton(ns("download_thresholded"), "Download thresholded network"),
          shiny::downloadButton(ns("download_plots"), "Download all plots (PDF)")
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
        basis <- sprintf("FDR < %s", fdr$target_fdr)
      } else {
        threshold_label <- 0
        basis <- "no cutoff applied -- showing the full network"
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

    # a single multi-page, vector PDF (one plot per page) via the same ggplot2
    # objects behind the interactive plotly views above -- plotly's own image
    # export rasterizes (and needs an external "kaleido" install), so this
    # goes straight to ggplot2 -> grDevices::pdf() for genuinely unbounded
    # resolution, matching what save_diagnostic_plots() does for CLI runs.
    output$download_plots <- shiny::downloadHandler(
      filename = function() "diagnostic_plots.pdf",
      content = function(file) {
        res <- network_result()
        shiny::req(res)
        grDevices::pdf(file, width = 7, height = 5)
        on.exit(grDevices::dev.off(), add = TRUE)
        # label_line(s) = TRUE -- a static PDF has no hover, so the cutoff/
        # target-FDR lines need their value written on the plot rather than
        # left to a tooltip (the interactive plotly views above rely on hover
        # instead, via the same functions' default label_line(s) = FALSE).
        print(plot_weight_distribution(res$network, cutoff = current_threshold(),
                                        title = "Edge weight distribution", label_line = TRUE))
        print(plot_outdegree_distribution(res$network, title = "Regulator out-degree distribution"))
        if (!is.null(res$fdr_result)) {
          print(plot_fdr_curve(base_fdr(), type = "curve", title = "FDR curve",
                                threshold = current_threshold(), show_target_fdr_line = is.null(manual_threshold()),
                                label_lines = TRUE))
          print(plot_fdr_curve(base_fdr(), type = "weight_comparison",
                                title = "Real vs. permuted weight distribution"))
        }
      }
    )

    display_network
  })
}
