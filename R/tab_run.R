################################################################################
# Module: Run (lives in the persistent dashboard sidebar, not a navbar tab --
# every other tab needs this run's result, so it stays visible regardless of
# which tab is active, matching protigy-v2's setup-sidebar convention).
################################################################################

#' A small "(i)" icon that shows `...` (pasted together) as a hover tooltip --
#' for a label that needs a sentence of explanation without permanently taking
#' up sidebar space the way a `helpText()` paragraph would. Pure CSS (see the
#' `.scion-tooltip` rule in [app_UI()]) rather than the native `title`
#' attribute, which is inconsistent across browsers (slow to appear, or not at
#' all for inline elements in some layouts).
#' @keywords internal
info_tooltip <- function(...) {
  shiny::tags$span(
    class = "scion-tooltip", `data-tooltip` = paste0(...),
    shiny::icon("circle-info"),
    style = "color: #888; margin-left: 4px; font-size: 0.85em;"
  )
}

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
    shiny::uiOutput(ns("data_inputs")),
    shiny::selectInput(ns("clustering_method"), "Clustering",
                        choices = c("None" = "none", "Temporal (DTW)" = "dtw",
                                    "Non-temporal (ICA)" = "ica", "Non-temporal (k-means)" = "kmeans",
                                    "Upload clusters" = "upload"),
                        selected = "none"),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] != 'none' && input['%s'] != 'upload'",
                           ns("clustering_method"), ns("clustering_method")),
      shiny::uiOutput(ns("clustering_data_input")),
      shiny::numericInput(ns("clustering_threshold"), "Clustering threshold", value = 0.5, min = 0, step = 0.1)
    ),
    shiny::conditionalPanel(
      condition = sprintf("input['%s'] == 'upload'", ns("clustering_method")),
      shiny::fileInput(ns("clusters_file"), "Pre-computed clusters file")
    ),
    shiny::checkboxInput(ns("connect_hubs"), "Connect cluster hubs", value = TRUE),
    shiny::checkboxInput(ns("normalize"),
                          shiny::tagList("Normalize edge weights", info_tooltip(
                            "Locked off during permutations -- it would cap every permutation's top weight at 1."
                          )), value = TRUE),
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
    ns <- session$ns
    network_result <- shiny::reactiveVal(NULL)
    using_example_data <- shiny::reactiveVal(FALSE)

    # run_scion() (and cluster_genes()/permute_network()) emit
    # "SCION_STAGE: <stage>" messages at the start/end of each major phase.
    # Maps each to a progress-bar position/label; stages that don't run for a
    # given call (e.g. no clustering, no permutations) are simply skipped over.
    # Permutation testing additionally emits "permutation testing progress
    # <done>/<total>" as each permutation (or, when parallelized, each worker-
    # sized batch) finishes -- handled separately below so the bar actually
    # advances across that whole phase instead of sitting at 0.65 until done.
    stage_progress <- list(
      "clustering started" = list(value = 0.05, detail = "Clustering..."),
      "clustering complete" = list(value = 0.3, detail = "Clustering complete"),
      "network inference started" = list(value = 0.35, detail = "Inferring network..."),
      "network inference complete" = list(value = 0.6, detail = "Network inference complete"),
      "permutation testing started" = list(value = 0.65, detail = "Running permutations..."),
      "permutation testing complete" = list(value = 0.95, detail = "Permutations complete")
    )
    permutation_progress_pattern <- "^permutation testing progress (\\d+)/(\\d+)$"

    # shared by both the "Run network" button and the example-data path: runs
    # run_scion() with whatever args it's given, reports the result, and
    # auto-navigates to Network Diagnostics. trace = FALSE always -- otherwise
    # RS.Get.Weight.Matrix() prints one message per target gene, which floods
    # the console/logs for anything beyond a handful of genes. Other,
    # non-stage messages (e.g. "No clustering_data supplied...") are left
    # alone and just print normally.
    run_and_report <- function(args) {
      args$trace <- FALSE
      shiny::withProgress(message = "Running SCION...", value = 0, {
        handle_stage_message <- function(m) {
          txt <- trimws(conditionMessage(m))
          if (startsWith(txt, "SCION_STAGE: ")) {
            stage <- sub("^SCION_STAGE: ", "", txt)
            progress_match <- regmatches(stage, regexec(permutation_progress_pattern, stage))[[1]]
            if (length(progress_match) == 3) {
              done <- as.integer(progress_match[2])
              total <- as.integer(progress_match[3])
              shiny::setProgress(value = 0.65 + 0.3 * (done / total),
                                  detail = sprintf("Running permutations... (%d/%d)", done, total))
            } else {
              info <- stage_progress[[stage]]
              if (!is.null(info)) {
                shiny::setProgress(value = info$value, detail = info$detail)
              }
            }
            invokeRestart("muffleMessage")
          }
        }

        result <- withCallingHandlers(
          tryCatch(do.call(run_scion, args), error = function(e) {
            shiny::showNotification(paste("SCION run failed:", conditionMessage(e)),
                                     type = "error", duration = NULL)
            NULL
          }),
          message = handle_stage_message
        )

        if (!is.null(result)) {
          shiny::setProgress(value = 1, detail = "Done")
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

    # "Load example data" only pre-fills parameters -- it does NOT run
    # anything. The user reviews/adjusts them, then clicks "Run network"
    # themselves, same as with their own uploaded files. Mirrors the settings
    # documented in the README/tutorial for this dataset: temporal (DTW)
    # clustering with the bundled clustering matrix. The tutorial's edge
    # cutoff (0.33) is applied afterward, on the Network Diagnostics tab --
    # the app always runs the full, unthresholded network (see below).
    shiny::observeEvent(input$load_example, {
      using_example_data(TRUE)
      shiny::updateCheckboxInput(session, "gene_list_header", value = TRUE)
      shiny::updateSelectInput(session, "clustering_method", selected = "dtw")
      # skipped when permutations are already on -- that lock keeps this at FALSE
      if (!isTRUE(input$run_permutations)) {
        shiny::updateCheckboxInput(session, "normalize", value = TRUE)
      }
      shiny::updateCheckboxInput(session, "connect_hubs", value = TRUE)
      shiny::updateSelectInput(session, "engine", selected = "randomForest")
      shiny::updateTextInput(session, "ptm_sep", value = ".")
      shiny::updateNumericInput(session, "seed", value = 2020)
      shiny::showNotification(
        "Example data selected -- review the parameters below, then click \"Run network\".",
        type = "message", duration = 8
      )
    })

    shiny::observeEvent(input$clear_example, using_example_data(FALSE))

    # Permutation testing needs the full, unnormalized network to compare
    # against -- per-network normalization would force every permutation's
    # top edge weight to 1 regardless of signal (see compute_fdr_threshold()).
    # Lock it off for the duration, restoring whatever the user had set when
    # they turn permutations back off. (No separate weight-cutoff lock needed
    # here -- the app never applies one at run time at all; see input$run below.)
    normalize_before_permute <- shiny::reactiveVal(TRUE)
    shiny::observeEvent(input$run_permutations, {
      if (isTRUE(input$run_permutations)) {
        normalize_before_permute(input$normalize)
        shiny::updateCheckboxInput(session, "normalize", value = FALSE)
        shinyjs::disable("normalize")
      } else {
        shinyjs::enable("normalize")
        shiny::updateCheckboxInput(session, "normalize", value = normalize_before_permute())
      }
    }, ignoreInit = TRUE)

    # uploading your own file stops using the example data
    shiny::observeEvent(input$target_data_file, using_example_data(FALSE), ignoreInit = TRUE)
    shiny::observeEvent(input$reg_data_file, using_example_data(FALSE), ignoreInit = TRUE)

    output$data_inputs <- shiny::renderUI({
      if (isTRUE(using_example_data())) {
        reg_type <- input$example_regulator_type
        shiny::div(
          class = "well", style = "padding: 10px 12px; color: #333;",
          shiny::tags$strong(shiny::icon("circle-check"), " Using bundled Arabidopsis example data"),
          shiny::tags$ul(
            style = "padding-left: 18px; margin: 6px 0; color: #333;",
            shiny::tags$li("Target matrix: target_mat_RNA.csv"),
            shiny::tags$li(sprintf("Regulator matrix: reg_mat_%s.csv", reg_type)),
            shiny::tags$li("Target gene list: target_list_RNA.csv"),
            shiny::tags$li(sprintf("Regulator gene list: reg_list_%s.csv", reg_type))
          ),
          shiny::actionButton(ns("clear_example"), "Use my own data instead",
                               icon = shiny::icon("rotate-left"), class = "btn-xs btn-default",
                               style = "margin-top: 4px;")
        )
      } else {
        shiny::tagList(
          shiny::fileInput(ns("target_data_file"), "Target matrix", accept = c(".csv", ".gct")),
          shiny::fileInput(ns("reg_data_file"), "Regulator matrix", accept = c(".csv", ".gct")),
          shiny::fileInput(ns("target_genes_file"), "Target gene list (optional)"),
          shiny::fileInput(ns("reg_genes_file"), "Regulator gene list (optional)"),
          shiny::checkboxInput(ns("gene_list_header"), "Gene lists have a header row", value = TRUE)
        )
      }
    })

    output$clustering_data_input <- shiny::renderUI({
      if (isTRUE(using_example_data())) {
        shiny::helpText(sprintf("Using bundled clustering matrix: cluster_mat_%s.csv",
                                 input$example_regulator_type))
      } else {
        shiny::fileInput(ns("clustering_data_file"),
                          "Clustering matrix (optional -- defaults to the combined target + regulator data)")
      }
    })

    shiny::observeEvent(input$run, {
      if (isTRUE(using_example_data())) {
        example_dir <- system.file("extdata", "arabidopsis", package = "SCION")
        shiny::validate(shiny::need(nzchar(example_dir), "Example data not found in the installed package."))
        reg_type <- input$example_regulator_type
        target_data_file <- file.path(example_dir, "target_mat_RNA.csv")
        reg_data_file <- file.path(example_dir, sprintf("reg_mat_%s.csv", reg_type))
        target_genes_file <- file.path(example_dir, "target_list_RNA.csv")
        reg_genes_file <- file.path(example_dir, sprintf("reg_list_%s.csv", reg_type))
        clustering_data_file <- if (input$clustering_method %in% c("dtw", "ica", "kmeans")) {
          file.path(example_dir, sprintf("cluster_mat_%s.csv", reg_type))
        }
        gene_list_header <- TRUE
      } else {
        shiny::req(input$target_data_file, input$reg_data_file)
        target_data_file <- input$target_data_file$datapath
        reg_data_file <- input$reg_data_file$datapath
        target_genes_file <- if (!is.null(input$target_genes_file)) input$target_genes_file$datapath
        reg_genes_file <- if (!is.null(input$reg_genes_file)) input$reg_genes_file$datapath
        clustering_data_file <- if (!is.null(input$clustering_data_file)) input$clustering_data_file$datapath
        gene_list_header <- input$gene_list_header
      }

      run_and_report(list(
        target_data_file = target_data_file,
        reg_data_file = reg_data_file,
        target_genes_file = target_genes_file,
        reg_genes_file = reg_genes_file,
        gene_list_header = gene_list_header,
        clustering_data_file = clustering_data_file,
        clustering_method = input$clustering_method,
        clustering_threshold = input$clustering_threshold,
        clusters_file = if (!is.null(input$clusters_file)) input$clusters_file$datapath,
        connect_hubs = input$connect_hubs,
        # always the full, unthresholded network -- apply a cutoff afterward,
        # on the Network Diagnostics tab, where it can be adjusted without
        # re-running inference (see "Threshold network" there)
        weightthreshold = 0,
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

    output$status <- shiny::renderUI({
      res <- network_result()
      if (is.null(res)) {
        return(NULL)
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
        },
        shiny::downloadButton(ns("download_full"), "Download full network")
      )
    })

    output$download_full <- shiny::downloadHandler(
      filename = function() "full_network.tsv",
      content = function(file) {
        res <- network_result()
        shiny::req(res)
        write_scion_network(res$network, file)
      }
    )

    network_result
  })
}
