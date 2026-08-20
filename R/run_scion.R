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
#' @param target_data_file,reg_data_file,target_genes_file,reg_genes_file,format
#'   passed to [read_scion_inputs()].
#' @param clustering_method passed to [cluster_genes()] as `method`: `"none"`
#'   (default), `"dtw"`, `"ica"`, `"kmeans"`, or `"upload"`.
#' @param clustering_data_file path to a clustering matrix, required unless
#'   `clustering_method` is `"none"` or `"upload"`.
#' @param clustering_threshold passed to [cluster_genes()] as `threshold`.
#' @param clusters_file passed to [cluster_genes()] as `clusters_file`, required
#'   when `clustering_method = "upload"`.
#' @param connect_hubs,weightthreshold,normalize,num.cores,engine,ptm_sep passed
#'   to [infer_network()].
#' @param seed RNG seed set once, before clustering and inference, so the same
#'   inputs produce the same network every run. Default matches the legacy
#'   `SCION()` behavior. Set to `NULL` to skip seeding.
#' @param output_file optional path to write the final edge table to (tab-
#'   separated, Cytoscape-importable). `NULL` (default) writes nothing.
#' @param ... additional arguments passed to [infer_network()] /
#'   [RS.Get.Weight.Matrix()].
#' @return a list with `network` (the edge table), `target`, `reg` (the
#'   processed input matrices), `cluster_assignment` (or `NULL`), and `params`
#'   (the arguments used, for reference / for feeding into [permute_network()]).
#' @export
run_scion <- function(target_data_file, reg_data_file, target_genes_file = NULL,
                       reg_genes_file = NULL, format = c("csv", "gct"),
                       clustering_method = c("none", "dtw", "ica", "kmeans", "upload"),
                       clustering_data_file = NULL, clustering_threshold = 0.5,
                       clusters_file = NULL, connect_hubs = TRUE, weightthreshold = 0,
                       normalize = TRUE, num.cores = 1,
                       engine = c("randomForest", "ranger"), ptm_sep = ".", seed = 2020,
                       output_file = NULL, ...) {
  format <- match.arg(format)
  clustering_method <- match.arg(clustering_method)
  engine <- match.arg(engine)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  inputs <- read_scion_inputs(target_data_file, reg_data_file, target_genes_file,
                               reg_genes_file, clustering_data_file, format)

  cluster_assignment <- cluster_genes(inputs$cluster_data, method = clustering_method,
                                       threshold = clustering_threshold,
                                       clusters_file = clusters_file, reg_data = inputs$reg)

  network <- infer_network(inputs$target, inputs$reg, cluster_assignment = cluster_assignment,
                            weightthreshold = weightthreshold, normalize = normalize,
                            connect_hubs = connect_hubs, num.cores = num.cores, engine = engine,
                            ptm_sep = ptm_sep, ...)

  if (!is.null(output_file)) {
    write_scion_network(network, output_file)
  }

  result <- list(network = network, target = inputs$target, reg = inputs$reg,
                  cluster_assignment = cluster_assignment,
                  params = list(weightthreshold = weightthreshold, normalize = normalize,
                                 connect_hubs = connect_hubs, num.cores = num.cores,
                                 engine = engine, ptm_sep = ptm_sep, seed = seed))
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
