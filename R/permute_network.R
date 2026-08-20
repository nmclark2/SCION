#' Shuffle a gene x sample matrix for permutation testing
#'
#' @param mat matrix/data frame, genes as rows, samples as columns.
#' @param dim `"col"`: independently shuffle each sample's values across genes.
#'   `"row"`: independently shuffle each gene's values across samples.
#' @return a matrix of the same dimensions and dimnames as `mat`, with values
#'   shuffled per the validated scheme -- row/column *labels* are restored to
#'   their original positions; only the underlying values move.
#' @keywords internal
shuffle_matrix <- function(mat, dim = c("col", "row")) {
  dim <- match.arg(dim)
  mat <- as.matrix(mat)
  orig_rownames <- rownames(mat)
  orig_colnames <- colnames(mat)
  shuffled <- if (dim == "col") apply(mat, 2, sample) else t(apply(mat, 1, sample))
  rownames(shuffled) <- orig_rownames
  colnames(shuffled) <- orig_colnames
  shuffled
}

#' Generate permuted null networks for FDR thresholding
#'
#' Shuffles the target/regulator matrices and reruns [infer_network()] `n_permutations`
#' times, reusing the SAME fixed `cluster_assignment` every time (clustering is never
#' recomputed per permutation -- see [cluster_genes()]). Ported from the validated
#' lab procedure (`scion_data_processing.R`): each permutation `i` is seeded with
#' `base_seed + i` before shuffling, which makes results independent of `num.cores`
#' for the same reason the per-target-gene seeding in [RS.Get.Weight.Matrix()] does.
#'
#' @param target,reg the processed target/regulator matrices used for the real
#'   network (the `target`/`reg` elements of [run_scion()]'s result).
#' @param cluster_assignment the fixed cluster assignment used for the real network
#'   (the `cluster_assignment` element of [run_scion()]'s result, or `NULL`).
#' @param n_permutations number of permutations (typically 100).
#' @param permute_dim passed to [shuffle_matrix()] as `dim`; `"col"` (default)
#'   matches the validated lab scheme.
#' @param base_seed permutation `i` uses seed `base_seed + i`. Default 0
#'   reproduces the lab scheme (`set.seed(i)`).
#' @param num.cores when `> 2`, permutations run across a FORK cluster (one
#'   permutation per worker); each permutation's own [infer_network()] call is
#'   then forced to `num.cores = 1` internally to avoid nesting parallelism.
#'   When `num.cores <= 2`, permutations run serially and their own
#'   `num.cores` is passed straight through to each [infer_network()] call.
#' @param weightthreshold,normalize,connect_hubs,ptm_sep passed to [infer_network()].
#' @param engine passed to [infer_network()]. **Must be the same engine used to
#'   produce `target`/`reg`'s real network** -- see [RS.Get.Weight.Matrix()] and
#'   [compute_fdr_threshold()].
#' @param ... additional arguments passed to [infer_network()].
#' @return a list of length `n_permutations`, each element an edge table as
#'   returned by [infer_network()].
#' @export
permute_network <- function(target, reg, cluster_assignment = NULL, n_permutations = 100,
                             permute_dim = c("col", "row"), base_seed = 0, num.cores = 1,
                             weightthreshold = 0, normalize = TRUE, connect_hubs = TRUE,
                             engine = c("randomForest", "ranger"), ptm_sep = ".", ...) {
  permute_dim <- match.arg(permute_dim)
  engine <- match.arg(engine)

  outer_parallel <- num.cores > 2
  inner_num_cores <- if (outer_parallel) 1 else num.cores

  run_one <- function(i) {
    set.seed(base_seed + i)
    shuffled_target <- shuffle_matrix(target, permute_dim)
    shuffled_reg <- shuffle_matrix(reg, permute_dim)
    infer_network(shuffled_target, shuffled_reg, cluster_assignment = cluster_assignment,
                   weightthreshold = weightthreshold, normalize = normalize,
                   connect_hubs = connect_hubs, num.cores = inner_num_cores, engine = engine,
                   ptm_sep = ptm_sep, ...)
  }

  if (outer_parallel) {
    clst <- parallel::makeCluster(num.cores - 1, type = "FORK")
    on.exit(parallel::stopCluster(clst), add = TRUE)
    parallel::parLapply(clst, seq_len(n_permutations), run_one)
  } else {
    lapply(seq_len(n_permutations), run_one)
  }
}
