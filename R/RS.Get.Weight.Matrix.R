#' Infer a weighted regulator -> target network via random-forest importance
#'
#' For each target gene, fits a random forest predicting its expression from all
#' regulator features, and records each regulator's importance as the edge
#' weight. This is the GENIE3-style core of SCION's network inference, and is
#' called identically for the real network and for every permutation generated
#' by [permute_network()].
#'
#' @param target.matrix data frame/matrix of target expression, samples as rows,
#'   genes as columns.
#' @param input.matrix data frame/matrix of regulator (input) expression, samples
#'   as rows, genes as columns.
#' @param K number of candidate regulators considered at each tree split: `"sqrt"`,
#'   `"all"`, or an integer.
#' @param nb.trees number of trees per random forest.
#' @param importance.measure `"IncNodePurity"` or `"%IncMSE"` (randomForest engine);
#'   mapped to `"impurity"`/`"permutation"` respectively when `engine = "ranger"`.
#' @param seed optional RNG seed for the top-level per-target seed draw. Must be
#'   set consistently between a real network and its permutations for the
#'   permutation seeding scheme in [permute_network()] to be reproducible.
#' @param trace if `TRUE`, emit a progress message per target gene (randomForest
#'   engine only).
#' @param normalize if `TRUE`, rescale the weight matrix to `[0, 1]`.
#' @param num.cores number of cores to use. For `engine = "randomForest"`, a
#'   `num.cores - 1` FORK cluster is used to parallelize across target genes
#'   (disabled when `num.cores <= 2`). For `engine = "ranger"`, target genes are
#'   processed serially and `num.cores` is instead passed to each ranger fit's
#'   own internal thread parallelism, since nesting FORK-based parallelism
#'   around an already-multithreaded ranger call would oversubscribe cores.
#' @param engine `"randomForest"` (default) or `"ranger"`. **Must be the same for
#'   a real network and every one of its permutations** -- the FDR calculation in
#'   [compute_fdr_threshold()] rank-matches edge weights between the real and
#'   permuted networks, which is only valid when they come from the same engine.
#' @param ... additional arguments passed to the underlying random forest call.
#' @return a numeric matrix of edge weights (targets x regulators), or `NULL` if
#'   there are no targets or no regulators.
#' @export
RS.Get.Weight.Matrix <- function(target.matrix, input.matrix, K = "sqrt", nb.trees = 10000,
                                  importance.measure = "%IncMSE", seed = NULL, trace = TRUE,
                                  normalize = TRUE, num.cores = 1,
                                  engine = c("randomForest", "ranger"), ...) {
  engine <- match.arg(engine)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }

  # normalize expression matrix
  target.matrix <- apply(target.matrix, 2, function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE))
  input.matrix <- apply(input.matrix, 2, function(x) (x - mean(x, na.rm = TRUE)) / stats::sd(x, na.rm = TRUE))
  input.matrix <- input.matrix[, !is.na(colSums(input.matrix))]

  num.samples <- dim(target.matrix)[1]
  num.targets <- dim(target.matrix)[2]
  num.inputs <- dim(input.matrix)[2]
  target.names <- colnames(target.matrix)
  input.names <- colnames(input.matrix)

  if (is.null(num.inputs) | is.null(num.targets)) {
    return(NULL)
  }

  weight.matrix <- matrix(0.0, nrow = num.targets, ncol = num.inputs)
  rownames(weight.matrix) <- target.names
  colnames(weight.matrix) <- input.names

  if (is.numeric(K)) {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.inputs))
  } else if (K == "all") {
    mtry <- num.inputs - 1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }

  names(target.names) <- target.names

  # one seed per target, drawn in the parent, so the forest for a given target is
  # the same whichever worker fits it -- and whether or not there is a worker at all
  target.seeds <- stats::setNames(sample.int(.Machine$integer.max, length(target.names)), target.names)

  if (engine == "ranger") {
    # ranger parallelizes internally (num.threads); looping across targets in an
    # outer FORK cluster on top of that would oversubscribe cores, so we process
    # targets serially here and let ranger fan out within each fit instead.
    imList <- lapply(target.names, function(x) {
      rsgwm2_ranger(x, target.matrix, input.matrix, mtry, nb.trees, importance.measure,
                    seed = target.seeds[[x]], num.threads = num.cores, ...)
    })
  } else if (num.cores > 2) {
    clst <- parallel::makeCluster(num.cores - 1, type = "FORK", outfile = "log.txt")
    doParallel::registerDoParallel(clst)
    imList <- parallel::parLapply(cl = clst, X = target.names, function(x) {
      rsgwm2_randomforest(x, num.targets, target.names, input.matrix, target.matrix, trace,
                          mtry, nb.trees, importance.measure, seed = target.seeds[[x]], ...)
    })
    parallel::stopCluster(cl = clst)
  } else {
    imList <- lapply(target.names, function(x) {
      rsgwm2_randomforest(x, num.targets, target.names, input.matrix, target.matrix, trace,
                          mtry, nb.trees, importance.measure, seed = target.seeds[[x]], ...)
    })
  }

  for (nm in names(imList)) {
    tcols <- names(imList[[nm]])
    weight.matrix[nm, tcols] <- imList[[nm]]
  }

  mynet <- weight.matrix / num.samples
  if (normalize) {
    mynet <- (mynet - min(mynet, na.rm = TRUE)) / (max(mynet, na.rm = TRUE) - min(mynet, na.rm = TRUE))
  }
  mynet
}

#' @keywords internal
rsgwm2_randomforest <- function(target.gene.name, num.targets, target.names, input.matrix,
                                 target.matrix, trace, mtry, nb.trees, importance.measure,
                                 seed = NULL, ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (trace) {
    target.gene.idx <- which(target.names == target.gene.name)
    message(sprintf("Computing gene %d/%d", target.gene.idx, num.targets))
  }

  # NOTE: the target gene is deliberately NOT removed from the input (regulator)
  # matrix, even if present there -- removing it breaks inference when there is
  # only one regulator in a network. This also means autoregulation is not
  # excluded from the network.
  x <- input.matrix
  y <- target.matrix[, target.gene.name]

  rf <- randomForest::randomForest(x = x, y = y, mtry = mtry, ntree = nb.trees,
                                    keep.forest = FALSE, importance = TRUE, ...)
  randomForest::importance(rf)[, importance.measure]
}

#' @keywords internal
rsgwm2_ranger <- function(target.gene.name, target.matrix, input.matrix, mtry, nb.trees,
                           importance.measure, seed = NULL, num.threads = 1, ...) {
  x <- input.matrix
  y <- target.matrix[, target.gene.name]
  ranger_importance <- if (importance.measure == "IncNodePurity") "impurity" else "permutation"

  rf <- ranger::ranger(x = x, y = y, mtry = mtry, num.trees = nb.trees,
                       importance = ranger_importance, num.threads = num.threads,
                       seed = seed, ...)
  rf$variable.importance
}
