#' Run SCION: read inputs, cluster once, infer a network
#'
#' User-facing wrapper around [read_scion_inputs()], [cluster_genes()], and
#' [infer_network()]. Unlike the legacy `SCION()` function, this does not
#' require a working directory or hardcode file output -- it returns
#' everything in memory, and optionally writes a Cytoscape-importable edge
#' table if `output_file` is given. The returned list (particularly `target`,
#' `reg`, and `cluster_assignment`) has everything [permute_network()] needs to
#' run permutations against this same network without recomputing clustering.
#'
#' @param target_data_file,reg_data_file,target_genes_file,reg_genes_file,
#'   gene_list_header,format passed to [read_scion_inputs()].
#' @param clustering_method passed to [cluster_genes()] as `method`: `"none"`
#'   (default), `"dtw"`, `"ica"`, `"kmeans"`, or `"upload"`.
#' @param clustering_data_file path to a clustering matrix, required unless
#'   `clustering_method` is `"none"` or `"upload"`.
#' @param clustering_threshold passed to [cluster_genes()] as `threshold`.
#' @param clusters_file passed to [cluster_genes()] as `clusters_file`, required
#'   when `clustering_method = "upload"`.
#' @param connect_hubs,num.cores,engine,ptm_sep passed to [infer_network()].
#' @param weightthreshold passed to [infer_network()]. Forced to `0` (with a warning) whenever
#'   `permute = TRUE`, regardless of what's passed -- the FDR comparison needs the full,
#'   unthresholded network on both sides; apply a cutoff to the result afterward instead (see
#'   [compute_fdr_threshold()]).
#' @param normalize passed to [infer_network()]. Forced to `FALSE` (with a warning) whenever
#'   `permute = TRUE`, regardless of what's passed -- normalizing rescales each network (real
#'   and every permutation) independently to `[0, 1]`, which would force every permutation's
#'   top edge weight to exactly 1 and invalidate the rank-based FDR comparison.
#' @param seed RNG seed set once, before clustering and inference, so the same
#'   inputs produce the same network every run. Default matches the legacy
#'   `SCION()` behavior. Set to `NULL` to skip seeding.
#' @param permute if `TRUE`, also run [permute_network()] and
#'   [compute_fdr_threshold()] against the just-inferred network, in this same
#'   call. Runs `n_permutations` permutations using the SAME `weightthreshold`,
#'   `normalize`, `connect_hubs`, `engine`, `ptm_sep`, `num.cores`, and `...`
#'   used for the real network above -- there is no separate way to set these
#'   for the permutations, since the FDR calculation requires them to match.
#'   For sharding permutations across an HPC job array instead of running them
#'   all in this one call, use `permute = FALSE` here and call
#'   [permute_network()] / [save_permutation()] / [load_permutations()]
#'   directly against this call's `target`/`reg`/`cluster_assignment`.
#' @param n_permutations,permute_dim,base_seed passed to [permute_network()]
#'   when `permute = TRUE`.
#' @param target_fdr passed to [compute_fdr_threshold()] when `permute = TRUE`.
#' @param output_file optional path to write the final edge table to (tab-
#'   separated, Cytoscape-importable). `NULL` (default) writes nothing. When
#'   `permute = TRUE`, writes the FDR-thresholded network, not the raw one.
#'   Also triggers [save_diagnostic_plots()], saved alongside it (same
#'   directory, named from `output_file`'s base name) -- there's no Shiny app
#'   to view/download them from interactively in this CLI/scripted path, so
#'   they're written automatically instead of requiring a separate call.
#' @param ... additional arguments passed to [infer_network()] /
#'   [RS.Get.Weight.Matrix()] (and, when `permute = TRUE`, to the internal
#'   [permute_network()] call as well, so e.g. `nb.trees` stays consistent
#'   between the real network and its permutations).
#' @return a list with `network` (the real edge table), `target`, `reg` (the
#'   processed input matrices), `cluster_assignment` (or `NULL`), `params`
#'   (the arguments used, for reference / for feeding into [permute_network()]
#'   yourself), and -- only when `permute = TRUE` -- `permuted_networks` (the
#'   [permute_network()] result), `fdr_result` (the [compute_fdr_threshold()]
#'   result), and `network_thresholded` (shorthand for
#'   `fdr_result$thresholded_network`).
#' @export
run_scion <- function(target_data_file, reg_data_file, target_genes_file = NULL,
                       reg_genes_file = NULL, gene_list_header = TRUE, format = c("auto", "csv", "gct"),
                       clustering_method = c("none", "dtw", "ica", "kmeans", "upload"),
                       clustering_data_file = NULL, clustering_threshold = 0.5,
                       clusters_file = NULL, connect_hubs = TRUE, weightthreshold = 0,
                       normalize = TRUE, num.cores = 1,
                       engine = c("randomForest", "ranger"), ptm_sep = ".", seed = 2020,
                       permute = FALSE, n_permutations = 100,
                       permute_dim = c("col", "row"), base_seed = 0, target_fdr = 0.05,
                       output_file = NULL, ...) {
  format <- match.arg(format)
  clustering_method <- match.arg(clustering_method)
  engine <- match.arg(engine)
  permute_dim <- match.arg(permute_dim)

  if (permute && isTRUE(normalize)) {
    warning("normalize = TRUE rescales each network (real and every permutation) to its own ",
            "[0, 1] range, which would force every permutation's top edge weight to 1 ",
            "regardless of its actual signal and invalidate the rank-based FDR comparison in ",
            "compute_fdr_threshold(). Running with normalize = FALSE for both the real network ",
            "and its permutations instead.", call. = FALSE)
    normalize <- FALSE
  }

  if (permute && !identical(weightthreshold, 0)) {
    warning("weightthreshold != 0 would apply the same manual cutoff to the real network and ",
            "every permutation BEFORE the FDR comparison, biasing which edges compute_fdr_",
            "threshold() ever gets to compare. Running with weightthreshold = 0 for both the ",
            "real network and its permutations instead -- apply a cutoff to the result ",
            "afterward via compute_fdr_threshold() or a manual filter.", call. = FALSE)
    weightthreshold <- 0
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  inputs <- read_scion_inputs(target_data_file, reg_data_file, target_genes_file = target_genes_file,
                               reg_genes_file = reg_genes_file, gene_list_header = gene_list_header,
                               clustering_data_file = clustering_data_file, format = format)

  cluster_assignment <- cluster_genes(inputs$cluster_data, method = clustering_method,
                                       threshold = clustering_threshold, clusters_file = clusters_file,
                                       target_data = inputs$target, reg_data = inputs$reg)

  message("SCION_STAGE: network inference started")
  network <- infer_network(inputs$target, inputs$reg, cluster_assignment = cluster_assignment,
                            weightthreshold = weightthreshold, normalize = normalize,
                            connect_hubs = connect_hubs, num.cores = num.cores, engine = engine,
                            ptm_sep = ptm_sep, ...)
  message("SCION_STAGE: network inference complete")

  result <- list(network = network, target = inputs$target, reg = inputs$reg,
                  cluster_assignment = cluster_assignment,
                  params = list(weightthreshold = weightthreshold, normalize = normalize,
                                 connect_hubs = connect_hubs, num.cores = num.cores,
                                 engine = engine, ptm_sep = ptm_sep, seed = seed))

  if (permute) {
    permuted_networks <- permute_network(inputs$target, inputs$reg,
                                          cluster_assignment = cluster_assignment,
                                          n_permutations = n_permutations,
                                          permute_dim = permute_dim, base_seed = base_seed,
                                          num.cores = num.cores,
                                          weightthreshold = weightthreshold, normalize = normalize,
                                          connect_hubs = connect_hubs, engine = engine,
                                          ptm_sep = ptm_sep, ...)
    fdr_result <- compute_fdr_threshold(network, permuted_networks, target_fdr = target_fdr)

    result$permuted_networks <- permuted_networks
    result$fdr_result <- fdr_result
    result$network_thresholded <- fdr_result$thresholded_network
  }

  if (!is.null(output_file)) {
    write_scion_network(if (permute) result$network_thresholded else network, output_file)
    save_diagnostic_plots(result, dirname(output_file), prefix = tools::file_path_sans_ext(basename(output_file)))
  }

  result
}

#' Write a SCION network to a Cytoscape-importable file
#'
#' @param network an edge table as returned in the `network` element of
#'   [run_scion()]'s result.
#' @param output_file path to write to (tab-separated).
#' @export
write_scion_network <- function(network, output_file) {
  utils::write.table(network, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
  invisible(output_file)
}
