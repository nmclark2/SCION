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
#' @param n_permutations number of permutations (typically 100). Ignored if
#'   `indices` is given explicitly.
#' @param indices which permutation indices to run; each index `i` is seeded
#'   with `base_seed + i`, independent of every other index and of `num.cores`
#'   (see Details). Defaults to `seq_len(n_permutations)`, i.e. all of them in
#'   one call. Pass a single index -- e.g. an HPC job array's task ID -- to run
#'   just one permutation per job; see [save_permutation()]/[load_permutations()]
#'   for a file-based workflow built around exactly that pattern.
#' @param permute_dim passed to [shuffle_matrix()] as `dim`; `"col"` (default)
#'   matches the validated lab scheme.
#' @param base_seed permutation `i` uses seed `base_seed + i`. Default 0
#'   reproduces the lab scheme (`set.seed(i)`).
#' @param num.cores when `> 2`, the requested indices run across a FORK cluster
#'   (one permutation per worker at a time); each permutation's own
#'   [infer_network()] call is then forced to `num.cores = 1` internally to
#'   avoid nesting parallelism. When `num.cores <= 2`, indices run serially (in
#'   this process or, for a single index, as whatever single HPC job called
#'   this) and their own `num.cores` is passed straight through to each
#'   [infer_network()] call.
#' @param weightthreshold passed to [infer_network()]. **Must match whatever was used for the
#'   real network** (see `normalize` below), and should generally be `0`: the FDR comparison in
#'   [compute_fdr_threshold()] needs the full, unthresholded network on both sides -- a nonzero
#'   value applies the same manual cutoff before that comparison ever happens, biasing which
#'   edges get compared (see [run_scion()], which enforces this automatically when calling this
#'   function with `permute = TRUE`).
#' @param connect_hubs,ptm_sep passed to [infer_network()].
#' @param normalize passed to [infer_network()]. **Must match whatever was used for the real
#'   network** this permutation set will be compared against in [compute_fdr_threshold()], and
#'   should generally be `FALSE` for both: normalizing rescales each network independently to
#'   `[0, 1]`, forcing every permutation's top edge weight to exactly 1 regardless of its
#'   actual signal, which invalidates the rank-based FDR comparison (see [run_scion()], which
#'   enforces this automatically when calling this function with `permute = TRUE`).
#' @param ... additional arguments passed to [infer_network()].
#' @return a list the same length as `indices`, each element an edge table as
#'   returned by [infer_network()], named by its permutation index (as a string).
#' @details
#' Emits a `"SCION_STAGE: permutation testing progress <done>/<total>"` message
#' after every completed permutation (serial path) or every completed batch of
#' `num.cores - 1` permutations (parallel path -- a single `parLapply()` call
#' blocks until it returns, so there's no way to observe progress *within* one
#' call; splitting the work into per-worker-sized batches trades a little
#' dispatch overhead for a progress update roughly every `num.cores - 1`
#' permutations). The Shiny app's Run tab listens for these to drive its
#' progress bar smoothly across the whole permutation phase instead of sitting
#' still until it's all done.
#' @export
permute_network <- function(target, reg, cluster_assignment = NULL, n_permutations = 100,
                             indices = seq_len(n_permutations), permute_dim = c("col", "row"),
                             base_seed = 0, num.cores = 1, weightthreshold = 0, normalize = TRUE,
                             connect_hubs = TRUE, ptm_sep = ".", ...) {
  permute_dim <- match.arg(permute_dim)

  outer_parallel <- num.cores > 2
  inner_num_cores <- if (outer_parallel) 1 else num.cores
  n_total <- length(indices)

  message("SCION_STAGE: permutation testing started")

  run_one <- function(i) {
    set.seed(base_seed + i)
    shuffled_target <- shuffle_matrix(target, permute_dim)
    shuffled_reg <- shuffle_matrix(reg, permute_dim)
    infer_network(shuffled_target, shuffled_reg, cluster_assignment = cluster_assignment,
                   weightthreshold = weightthreshold, normalize = normalize,
                   connect_hubs = connect_hubs, num.cores = inner_num_cores,
                   ptm_sep = ptm_sep, ...)
  }

  results <- if (outer_parallel) {
    clst <- parallel::makeCluster(num.cores - 1, type = "FORK")
    on.exit(parallel::stopCluster(clst), add = TRUE)

    # one parLapply() call per worker-sized batch, not one call for everything --
    # a single call blocks until ALL of its indices finish, so batching is what
    # makes a progress message possible at all in the parallel case.
    batch_size <- num.cores - 1
    batch_starts <- seq(1, n_total, by = batch_size)
    n_done <- 0
    batches <- lapply(batch_starts, function(start) {
      batch <- indices[start:min(start + batch_size - 1, n_total)]
      res <- parallel::parLapply(clst, batch, run_one)
      n_done <<- n_done + length(batch)
      message(sprintf("SCION_STAGE: permutation testing progress %d/%d", n_done, n_total))
      res
    })
    do.call(c, batches)
  } else {
    lapply(seq_along(indices), function(pos) {
      res <- run_one(indices[pos])
      message(sprintf("SCION_STAGE: permutation testing progress %d/%d", pos, n_total))
      res
    })
  }

  message("SCION_STAGE: permutation testing complete")
  stats::setNames(results, as.character(indices))
}
