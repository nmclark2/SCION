#' Convert a weight matrix into a Cytoscape-style edge table
#'
#' @param network numeric matrix of edge weights, targets as rows, regulators as columns.
#' @param weightthreshold edges with weight strictly below this value are dropped.
#' @return a data frame with columns `Regulator`, `Interaction`, `Target`, `Weight`, in
#'   target-major, then regulator-major order (matching a row-by-row, column-by-column scan
#'   of `network`).
#' @keywords internal
weight_matrix_to_edges <- function(network, weightthreshold) {
  trimmed <- network
  trimmed[trimmed < weightthreshold] <- NA
  keep <- which(!is.na(trimmed), arr.ind = TRUE)
  keep <- keep[order(keep[, "row"], keep[, "col"]), , drop = FALSE]
  data.frame(
    Regulator = colnames(trimmed)[keep[, "col"]],
    Interaction = rep("regulates", nrow(keep)),
    Target = rownames(trimmed)[keep[, "row"]],
    Weight = trimmed[keep],
    stringsAsFactors = FALSE
  )
}

#' Pick the hub gene(s) (highest out-degree regulator) from an edge table
#'
#' @param edge_table a data frame as returned by [weight_matrix_to_edges()].
#' @return a character vector of regulator names with the maximum out-degree
#'   (length > 1 if there is a tie), or `character(0)` if `edge_table` has no rows.
#' @keywords internal
pick_hub_genes <- function(edge_table) {
  if (nrow(edge_table) == 0) {
    return(character(0))
  }
  edge_counts <- table(edge_table$Regulator)
  names(edge_counts)[edge_counts == max(edge_counts)]
}

#' Infer a single (non-clustered) network
#' @keywords internal
infer_network_single <- function(target_data, reg_data, weightthreshold, normalize,
                                  num.cores, ...) {
  if (dim(target_data)[1] < 1 || dim(reg_data)[1] < 1) {
    message("Need at least one target and at least one regulator to infer a network. SCION will not run")
    return(NULL)
  }
  network <- RS.Get.Weight.Matrix(t(target_data), t(reg_data), normalize = normalize,
                                  num.cores = num.cores, ...)
  if (is.null(network)) {
    return(NULL)
  }
  weight_matrix_to_edges(network, weightthreshold)
}

#' Infer one network per cluster, optionally connecting cluster hubs
#' @keywords internal
infer_network_clustered <- function(target_data, reg_data, cluster_assignment, weightthreshold,
                                      normalize, connect_hubs, num.cores, ptm_sep, ...) {
  finalnetwork <- data.frame(Regulator = character(), Interaction = character(),
                              Target = character(), Weight = double(), stringsAsFactors = FALSE)
  myhubs <- character(0)

  n_clusters <- max(cluster_assignment$clusters)
  # Draw one seed per cluster (plus one for the hub network) up front, before any
  # RS.Get.Weight.Matrix() call. A cluster's own output is already num.cores-invariant
  # given its seed (see RS.Get.Weight.Matrix()'s per-target seeding) -- but without this,
  # cluster i+1's *own* seed draw would depend on whatever ambient RNG state cluster i's
  # call happened to leave behind, which differs between forked (num.cores > 2) and serial
  # execution. Fixing that here makes the whole multi-cluster pipeline num.cores-invariant,
  # not just each individual cluster.
  cluster_seeds <- sample.int(.Machine$integer.max, n_clusters + 1)

  for (i in seq_len(n_clusters)) {
    mygenes <- row.names(cluster_assignment)[cluster_assignment$clusters == i]
    clustertargetdata <- target_data[row.names(target_data) %in% mygenes, , drop = FALSE]
    clusterregdata <- reg_data[row.names(reg_data) %in% mygenes, , drop = FALSE]

    # randomForest's %IncMSE importance is deterministically NA with exactly one
    # regulator (confirmed empirically -- not data-dependent), so a 1-regulator
    # cluster can never produce a usable edge; skip it outright rather than
    # silently getting an all-NA weight matrix from RS.Get.Weight.Matrix().
    if (dim(clusterregdata)[1] <= 1) {
      next
    }
    # GENIE3 cannot infer autoregulation on one TF, and errors with exactly two
    # (since one gets removed) -- so at least 3 TFs are needed when TFs are also
    # targets. Only relevant when TFs overlap with targets.
    if (sum(row.names(clusterregdata) %in% row.names(clustertargetdata)) > 0 &&
          dim(clusterregdata)[1] <= 2) {
      next
    }
    if (dim(clustertargetdata)[1] < 1 || dim(clusterregdata)[1] < 1) {
      next
    }

    network <- RS.Get.Weight.Matrix(t(clustertargetdata), t(clusterregdata), normalize = normalize,
                                     num.cores = num.cores, seed = cluster_seeds[i], ...)
    if (is.null(network)) {
      next
    }

    networktable <- weight_matrix_to_edges(network, weightthreshold)
    finalnetwork <- rbind(finalnetwork, networktable)
    myhubs <- c(myhubs, pick_hub_genes(networktable))
  }

  if (connect_hubs && length(myhubs) > 2) {
    hub_network <- infer_hub_network(target_data, reg_data, myhubs, weightthreshold, normalize,
                                      num.cores, ptm_sep, cluster_seeds[n_clusters + 1], ...)
    finalnetwork <- rbind(finalnetwork, hub_network)
  }

  finalnetwork
}

#' Infer the network connecting cluster hub genes
#' @keywords internal
infer_hub_network <- function(target_data, reg_data, myhubs, weightthreshold, normalize,
                               num.cores, ptm_sep, seed, ...) {
  hubtargetdata <- target_data[row.names(target_data) %in% myhubs, , drop = FALSE]
  hubregdata <- reg_data[row.names(reg_data) %in% myhubs, , drop = FALSE]
  if (dim(hubtargetdata)[1] == 0) {
    # strip PTM site information (assumes exactly one separator per hub name) to find targets
    genes <- unlist(strsplit(myhubs, ptm_sep, fixed = TRUE))
    genes <- genes[seq(1, length(genes), by = 2)]
    hubtargetdata <- target_data[row.names(target_data) %in% genes, , drop = FALSE]
  }
  if (dim(hubregdata)[1] <= 1) {
    return(NULL)
  }
  network <- RS.Get.Weight.Matrix(t(hubtargetdata), t(hubregdata), normalize = normalize,
                                   num.cores = num.cores, seed = seed, ...)
  if (is.null(network)) {
    return(NULL)
  }
  weight_matrix_to_edges(network, weightthreshold)
}

#' Infer a (possibly clustered) network from target/regulator matrices
#'
#' The core network-inference step of SCION, with no file I/O: given target and
#' regulator expression matrices (genes as rows, samples as columns) and an
#' optional, already-computed cluster assignment, infers one network per
#' cluster (or a single network if `cluster_assignment` is `NULL`), optionally
#' connects cluster hubs, and returns a single edge table. This function is
#' called identically by [run_scion()] for the real network and by
#' [permute_network()] for every permutation -- it never recomputes clustering.
#'
#' @param target_data data frame/matrix, genes as rows, samples as columns.
#' @param reg_data data frame/matrix, genes as rows, samples as columns.
#' @param cluster_assignment optional data frame with a `clusters` column (row
#'   names = gene names), as returned by [cluster_genes()]. `NULL` means infer a
#'   single network across all genes.
#' @param weightthreshold edges with weight below this value are dropped.
#' @param normalize passed to [RS.Get.Weight.Matrix()].
#' @param connect_hubs if clustering, whether to additionally infer a network
#'   connecting each cluster's hub gene (highest out-degree regulator).
#' @param num.cores passed to [RS.Get.Weight.Matrix()].
#' @param ptm_sep separator used to split a PTM-site regulator name into gene
#'   symbol + site (e.g. `"."` for `SOX2.S35`, `"_"` for `MEF2C_S453s`): when
#'   connecting cluster hubs and no hub gene is directly present in
#'   `target_data`'s row names, and always in the returned edge table (see
#'   [split_ptm_sites()]) so a PTM-annotated regulator and its plain-symbol
#'   target are recognized as the same gene.
#' @param seed optional RNG seed set once, at the very start of this call
#'   (before clustering-aware seed draws or any inference). Without it, the
#'   result depends on the ambient global RNG state at call time -- pass an
#'   explicit `seed` (or call `set.seed()` yourself right before calling) any
#'   time you need two calls to be comparable/reproducible.
#' @param ... additional arguments passed to [RS.Get.Weight.Matrix()].
#' @return a data frame with columns `Regulator`, `Interaction`, `Target`, `Weight`,
#'   or `NULL` if inference could not be run (e.g. no targets/regulators).
#' @export
infer_network <- function(target_data, reg_data, cluster_assignment = NULL,
                           weightthreshold = 0, normalize = TRUE, connect_hubs = TRUE,
                           num.cores = 1, ptm_sep = ".", seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  result <- if (is.null(cluster_assignment)) {
    infer_network_single(target_data, reg_data, weightthreshold, normalize, num.cores, ...)
  } else {
    infer_network_clustered(target_data, reg_data, cluster_assignment, weightthreshold, normalize,
                             connect_hubs, num.cores, ptm_sep, ...)
  }
  split_ptm_sites(result, ptm_sep)
}

#' Split a PTM-site regulator name into gene symbol + site, if present
#'
#' Regulator names get matched against target gene names elsewhere (hub
#' connection, network visualization, out-degree counting); a PTM-annotated
#' regulator like `SOX2.S35` would otherwise show up as a node distinct from
#' the plain-gene-symbol target `SOX2`, when they're the same underlying gene.
#' This strips the site off `Regulator` into its own `Site` column so
#' `Regulator` is always just the gene symbol.
#'
#' @param network an edge table with a `Regulator` column, or `NULL`.
#' @param ptm_sep separator to split on (e.g. `"."` for `SOX2.S35`, `"_"` for
#'   `MEF2C_S453s`). Only regulators actually containing `ptm_sep` are split;
#'   others are left as-is. If none contain it, `network` is returned
#'   unchanged (no `Site` column added) -- so non-PTM data keeps exactly the
#'   `Regulator`/`Interaction`/`Target`/`Weight` schema it always has.
#' @return `network` with `Regulator` reduced to gene symbols and a `Site`
#'   column added (`NA` for regulators that had no `ptm_sep`), or `network`
#'   unchanged if it's `NULL`, empty, or no `Regulator` contains `ptm_sep`.
#' @keywords internal
split_ptm_sites <- function(network, ptm_sep) {
  if (is.null(network) || nrow(network) == 0 || !any(grepl(ptm_sep, network$Regulator, fixed = TRUE))) {
    return(network)
  }
  parts <- strsplit(network$Regulator, ptm_sep, fixed = TRUE)
  gene <- vapply(parts, `[`, character(1), 1)
  site <- vapply(parts, function(p) if (length(p) >= 2) paste(p[-1], collapse = ptm_sep) else NA_character_,
                  character(1))
  data.frame(Regulator = gene, Site = site, Interaction = network$Interaction,
             Target = network$Target, Weight = network$Weight, stringsAsFactors = FALSE)
}
