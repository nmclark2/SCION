test_that("Shiny app: running a network shows weight/outdegree diagnostics (no permutations)", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- tempfile("scion-app-")
  dir.create(app_dir)
  writeLines("SCION::launchApp()", file.path(app_dir, "app.R"))

  mats <- make_test_matrices(n_targets = 15, n_regs = 6, n_samples = 10)
  target_file <- file.path(app_dir, "target.csv")
  reg_file <- file.path(app_dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)

  app <- shinytest2::AppDriver$new(app_dir, name = "scion-smoke", height = 900, width = 1400,
                                    load_timeout = 30000)
  on.exit(app$stop(), add = TRUE)

  # default tab on load is Help, not Network Diagnostics
  expect_equal(app$get_value(input = "navbar-tabs"), "help")

  # Run tab (sidebar): upload the small test matrices and run the real network
  app$upload_file(`run-target_data_file` = target_file)
  app$upload_file(`run-reg_data_file` = reg_file)
  app$set_inputs(`run-num_cores` = 1, wait_ = FALSE)
  app$click("run-run")
  app$wait_for_idle(timeout = 60000)

  status_html <- app$get_html("#run-status")
  expect_match(status_html, "edges")
  expect_match(status_html, "Download full network")

  # a completed run should auto-navigate away from the default Help tab
  expect_equal(app$get_value(input = "navbar-tabs"), "diagnostics")

  # weight/outdegree distributions show without ever running permutations
  expect_match(app$get_html("#diagnostics-weight_distribution"), "plotly")
  expect_match(app$get_html("#diagnostics-outdegree_distribution"), "plotly")

  # no FDR result yet, so no FDR curve/weight-comparison plots -- but a flat
  # weight cutoff is still available (this is the point: thresholding by
  # weight shouldn't require having run permutations)
  expect_no_match(app$get_html("#diagnostics-fdr_plots"), "FDR curve")
  summary_html <- app$get_html("#diagnostics-summary")
  expect_match(summary_html, "no cutoff applied")

  app$set_inputs(`diagnostics-manual_threshold` = 0.2, wait_ = FALSE)
  app$click("diagnostics-apply_manual_threshold")
  app$wait_for_idle(timeout = 10000)
  summary_html <- app$get_html("#diagnostics-summary")
  expect_match(summary_html, "Threshold: 0.2 \\(manual cutoff\\)")

  # Visualize tab -- static plot (no visNetwork widget in a headless snapshot)
  app$set_inputs(`visualize-interactive` = FALSE)
  app$wait_for_idle(timeout = 30000)
})

test_that("Shiny app: sidebar's 'Run permutations' checkbox populates the Network Diagnostics FDR section", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- tempfile("scion-app-")
  dir.create(app_dir)
  writeLines("SCION::launchApp()", file.path(app_dir, "app.R"))

  mats <- make_test_matrices(n_targets = 15, n_regs = 6, n_samples = 10)
  target_file <- file.path(app_dir, "target.csv")
  reg_file <- file.path(app_dir, "reg.csv")
  utils::write.csv(mats$target, target_file)
  utils::write.csv(mats$reg, reg_file)

  app <- shinytest2::AppDriver$new(app_dir, name = "scion-smoke-sidebar-permute", height = 900,
                                    width = 1400, load_timeout = 30000)
  on.exit(app$stop(), add = TRUE)

  app$upload_file(`run-target_data_file` = target_file)
  app$upload_file(`run-reg_data_file` = reg_file)
  app$set_inputs(`run-num_cores` = 1, wait_ = FALSE)
  app$set_inputs(`run-run_permutations` = TRUE)
  app$set_inputs(`run-n_permutations` = 3, wait_ = FALSE)
  app$click("run-run")
  app$wait_for_idle(timeout = 60000)

  # the sidebar itself should report the FDR threshold
  run_status_html <- app$get_html("#run-status")
  expect_match(run_status_html, "FDR threshold")

  # ... and the Network Diagnostics tab (auto-navigated to) should show it too,
  # with no click on any of its own controls needed -- it no longer has any
  summary_html <- app$get_html("#diagnostics-summary")
  expect_match(summary_html, "Threshold")
})

test_that("Shiny app: 'Load example data' only pre-fills parameters, does not auto-run", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- tempfile("scion-app-")
  dir.create(app_dir)
  writeLines("SCION::launchApp()", file.path(app_dir, "app.R"))

  app <- shinytest2::AppDriver$new(app_dir, name = "scion-smoke-example-data-noautorun", height = 900,
                                    width = 1400, load_timeout = 30000)
  on.exit(app$stop(), add = TRUE)

  app$click("run-load_example", wait_ = FALSE)
  app$wait_for_idle(timeout = 10000)

  # parameters got pre-filled, including clustering (a clustering file is bundled)...
  expect_equal(app$get_value(input = "run-clustering_method"), "dtw")
  # ...but nothing has actually run, and we're still on the default Help tab
  expect_equal(app$get_value(input = "navbar-tabs"), "help")
  expect_no_match(app$get_html("#run-status"), "edges")
})

test_that("Shiny app: 'Use my own data instead' reverts the example-data summary back to file inputs", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- tempfile("scion-app-")
  dir.create(app_dir)
  writeLines("SCION::launchApp()", file.path(app_dir, "app.R"))

  app <- shinytest2::AppDriver$new(app_dir, name = "scion-smoke-clear-example", height = 900,
                                    width = 1400, load_timeout = 30000)
  on.exit(app$stop(), add = TRUE)

  app$click("run-load_example", wait_ = FALSE)
  app$wait_for_idle(timeout = 10000)
  expect_match(app$get_html("#run-data_inputs"), "Using bundled Arabidopsis example data")

  app$click("run-clear_example", wait_ = FALSE)
  app$wait_for_idle(timeout = 10000)
  data_inputs_html <- app$get_html("#run-data_inputs")
  expect_no_match(data_inputs_html, "Using bundled Arabidopsis example data")
  expect_match(data_inputs_html, "Target matrix")
})

test_that("Shiny app: 'Load example data' + 'Run network' runs the bundled Arabidopsis network", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")

  app_dir <- tempfile("scion-app-")
  dir.create(app_dir)
  writeLines("SCION::launchApp()", file.path(app_dir, "app.R"))

  app <- shinytest2::AppDriver$new(app_dir, name = "scion-smoke-example-data", height = 900,
                                    width = 1400, load_timeout = 30000)
  on.exit(app$stop(), add = TRUE)

  app$set_inputs(`run-num_cores` = 6, wait_ = FALSE)
  app$click("run-load_example", wait_ = FALSE)
  app$wait_for_idle(timeout = 10000)

  # example-data summary shows (with the bundled files listed), and the user still
  # has to click "Run network" themselves
  data_inputs_html <- app$get_html("#run-data_inputs")
  expect_match(data_inputs_html, "Using bundled Arabidopsis example data")
  expect_match(data_inputs_html, "reg_mat_protein.csv")
  # the bundled clustering matrix note shows too, since clustering is on by default
  expect_match(app$get_html("#run-clustering_data_input"), "cluster_mat_protein.csv")

  app$click("run-run", wait_ = FALSE)
  app$wait_for_idle(timeout = 180000)

  status_html <- app$get_html("#run-status")
  expect_match(status_html, "edges")
  expect_equal(app$get_value(input = "navbar-tabs"), "diagnostics")
  expect_match(app$get_html("#diagnostics-weight_distribution"), "plotly")
})
