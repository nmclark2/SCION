################################################################################
# Module: Run (lives in the persistent dashboard sidebar, not a navbar tab --
# every other tab needs this run's result, so it stays visible regardless of
# which tab is active, matching protigy-v2's setup-sidebar convention).
################################################################################

#' @keywords internal
runSidebarUI <- function(id = "run") {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::selectInput(ns("example_regulator_type"), "Example regulator type",
                        choices = c("Protein" = "protein", "Phospho" = "phospho"),
                        selected = "protein"),
    shiny::actionButton(ns("load_example"), "Load example data (Arabidopsis)",
                         icon = shiny::icon("flask"), class = "btn-default"),
    shiny::hr(),
    shiny::fileInput(ns("target_data_file"), "Target matrix", accept = c(".csv", ".gct")),
    shiny::fileInput(ns("reg_data_file"), "Regulator matrix", accept = c(".csv", ".gct")),
    shiny::fileInput(ns("target_genes_file"), "Target gene list (optional)"),
    shiny::fileInput(ns("reg_genes_file"), "Regulator gene list (optional)"),
    shiny::checkboxInput(ns("gene_list_header"), "Gene lists have a header row", value = TRUE),
    shiny::selectInput(ns("clustering_method"), "Clustering",
                        choices = c("None" = "none", "Temporal (DTW)" = "dtw",
                                    "Non-temporal (ICA)" = "ica", "Non-temporal (k-means)" = "kmeans",
                                    "Upload clusters" = "upload"),
                        selected = "none"),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] != 'none' && input['%s'] != 'upload'",
                           ns("clustering_method"), ns("clustering_method")),
      shiny::fileInput(ns("clustering_data_file"),
                        "Clustering matrix (optional -- defaults to the combined target + regulator data)"),
      shiny::numericInput(ns("clustering_threshold"), "Clustering threshold", value = 0.5, min = 0, step = 0.1)
    ),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("clustering_method")),
      shiny::fileInput(ns("clusters_file"), "Pre-computed clusters file")
    ),
    shiny::checkboxInput(ns("connect_hubs"), "Connect cluster hubs", value = TRUE),
    shiny::numericInput(ns("weightthreshold"), "Edge weight cutoff", value = 0, min = 0, step = 0.1),
    shiny::checkboxInput(ns("normalize"), "Normalize edge weights", value = TRUE),
    shiny::selectInput(ns("engine"), "Random forest engine", choices = c("randomForest", "ranger"),
                        selected = "randomForest"),
    shiny::textInput(ns("ptm_sep"), "PTM site separator", value = "."),
    shiny::numericInput(ns("num_cores"), "Number of cores", value = 1, min = 1, step = 1),
    shiny::numericInput(ns("seed"), "Random seed", value = 2020, step = 1),
    shiny::checkboxInput(ns("run_permutations"), "Run permutations", value = FALSE),
    shiny::conditionalPanel(
      condition = sprintf("input['%s']", ns("run_permutations")),
      shiny::numericInput(ns("n_permutations"), "Number of permutations", value = 100, min = 1, step = 1),
      shiny::selectInput(ns("permute_dim"), "Permute dimension", choices = c("col", "row"), selected = "col"),
      shiny::numericInput(ns("base_seed"), "Base seed", value = 0, step = 1),
      shiny::numericInput(ns("target_fdr"), "Target FDR", value = 0.05, min = 0, max = 1, step = 0.01)
    ),
    shiny::actionButton(ns("run"), "Run network", icon = shiny::icon("play"), class = "btn-primary"),
    shiny::hr(),
    shiny::uiOutput(ns("status"))
  )
}

#' @param parent_session the top-level app session (from `app_server()`), used to
#'   switch the navbar to the Network Diagnostics tab once a run completes --
#'   `updateNavbarPage()` must target the app's actual navbar input, which lives
#'   outside this module's own namespace.
#' @noRd
runSidebarServer <- function(id = "run", parent_session) {
  shiny::moduleServer(id, function(input, output, session) {
    network_result <- shiny::reactiveVal(NULL)

    # shared by both the "Run network" button and "Load example data": runs
    # run_scion() with whatever args it's given, reports the result, and
    # auto-navigates to Network Diagnostics.
    run_and_report <- function(args) {
      shiny::withProgress(message = "Running SCION...", value = 0.3, {
        result <- tryCatch(do.call(run_scion, args), error = function(e) {
          shiny::showNotification(paste("SCION run failed:", conditionMessage(e)),
                                   type = "error", duration = NULL)
          NULL
        })
        if (!is.null(result)) {
          network_result(result)
          msg <- sprintf("Network inferred: %d edges.", nrow(result$network))
          if (isTRUE(args$permute)) {
            msg <- paste(msg, sprintf("FDR threshold: %s.",
                                       if (is.na(result$fdr_result$threshold)) "not reached"
                                       else signif(result$fdr_result$threshold, 4)))
          }
          shiny::showNotification(msg, type = "message", duration = 8)
          shiny::updateNavbarPage(parent_session, "navbar-tabs", selected = "diagnostics")
        }
      })
    }

    shiny::observeEvent(input$run, {
      shiny::req(input$target_data_file, input$reg_data_file)
      run_and_report(list(
        target_data_file = input$target_data_file$datapath,
        reg_data_file = input$reg_data_file$datapath,
        target_genes_file = if (!is.null(input$target_genes_file)) input$target_genes_file$datapath,
        reg_genes_file = if (!is.null(input$reg_genes_file)) input$reg_genes_file$datapath,
        gene_list_header = input$gene_list_header,
        clustering_data_file = if (!is.null(input$clustering_data_file)) input$clustering_data_file$datapath,
        clustering_method = input$clustering_method,
        clustering_threshold = input$clustering_threshold,
        clusters_file = if (!is.null(input$clusters_file)) input$clusters_file$datapath,
        connect_hubs = input$connect_hubs,
        weightthreshold = input$weightthreshold,
        normalize = input$normalize,
        num.cores = input$num_cores,
        engine = input$engine,
        ptm_sep = input$ptm_sep,
        seed = input$seed,
        permute = input$run_permutations,
        n_permutations = input$n_permutations,
        permute_dim = input$permute_dim,
        base_seed = input$base_seed,
        target_fdr = input$target_fdr
      ))
    })

    shiny::observeEvent(input$load_example, {
      example_dir <- system.file("extdata", "arabidopsis", package = "SCION")
      shiny::validate(shiny::need(nzchar(example_dir), "Example data not found in the installed package."))
      reg_type <- input$example_regulator_type
      run_and_report(list(
        target_data_file = file.path(example_dir, "target_mat_RNA.csv"),
        reg_data_file = file.path(example_dir, sprintf("reg_mat_%s.csv", reg_type)),
        target_genes_file = file.path(example_dir, "target_list_RNA.csv"),
        reg_genes_file = file.path(example_dir, sprintf("reg_list_%s.csv", reg_type)),
        gene_list_header = TRUE,
        clustering_method = "none", # kept fast for a demo; DTW in particular is O(n^2) on ~1100 genes
        connect_hubs = TRUE,
        weightthreshold = 0.33, # matches the published tutorial's documented default
        normalize = TRUE,
        num.cores = input$num_cores,
        engine = "randomForest",
        seed = 2020,
        permute = FALSE
      ))
    })

    output$status <- shiny::renderUI({
      res <- network_result()
      if (is.null(res)) {
        return(shiny::helpText("No network yet -- upload files and click Run network,",
                                "or click \"Load example data\"."))
      }
      shiny::tagList(
        shiny::strong(sprintf("%d edges", nrow(res$network))),
        if (!is.null(res$cluster_assignment)) {
          shiny::p(sprintf("%d clusters", max(res$cluster_assignment$clusters)))
        },
        if (!is.null(res$fdr_result)) {
          shiny::p(sprintf(
            "FDR threshold: %s (%d edges kept) -- see the Network Diagnostics tab.",
            if (is.na(res$fdr_result$threshold)) "not reached" else signif(res$fdr_result$threshold, 4),
            nrow(res$network_thresholded)
          ))
        }
      )
    })

    network_result
  })
}
