################################################################################
# Module: Help (static content)
################################################################################

#' A label/text pair matching the sidebar's own control order, so this tab
#' can double as a full parameter reference without turning back into a wall
#' of prose -- one short point per parameter.
#' @keywords internal
help_item <- function(label, text) {
  shiny::tagList(shiny::strong(label), shiny::p(text))
}

#' @keywords internal
helpTabUI <- function(id = "help") {
  shiny::tagList(
    shiny::h3("Run network (sidebar)"),
    help_item("Load example data", "Loads a real example dataset so you can try the app before
               uploading your own data."),
    help_item("Target matrix", "Your input data for the targets (e.g. transcript data). Rows are
               genes, columns are samples. CSV or GCT format."),
    help_item("Regulator matrix", "Your input data for the regulators (e.g. protein or PTM data).
               Same format as the target matrix. If a regulator name includes a PTM site
               (e.g. SOX2.S35), set the separator below."),
    help_item("Target / Regulator gene list", "Optional: restrict to specific genes (first column
               = gene names). Leave blank to use every gene in the matrix."),
    help_item("Gene lists have a header row", "Check this only if your gene list files have a
               header row -- a plain one-gene-per-line file needs it unchecked."),
    help_item("Clustering", "Optionally group genes before inference: temporal (DTW), non-temporal
               (ICA or k-means), or your own pre-computed clusters."),
    help_item("Clustering matrix", "Data used to compute clusters. Defaults to your target and
               regulator matrices combined if left blank."),
    help_item("Clustering threshold", "Cutoff used by DTW/ICA to decide cluster membership. Not
               used for k-means, which picks its own number of clusters."),
    help_item("Pre-computed clusters file", "A file assigning each gene to a cluster, if you
               already have one."),
    help_item("Connect cluster hubs", "Adds edges between each cluster's hub gene (its highest
               out-degree regulator) and the other clusters' hubs, linking the separate cluster
               networks together."),
    help_item("Normalize edge weights", "Rescales edge weights to [0, 1]. Not available when
               running permutations. This does not filter the network -- use \"Threshold network\"
               on the Network Diagnostics tab for that."),
    help_item("PTM site separator", "Character separating a gene name from its PTM site in
               regulator names (e.g. \".\" for SOX2.S35)."),
    help_item("Number of cores", "Number of CPU cores to use. Leave at 1 to disable
               parallelization."),
    help_item("Random seed", "Seed for reproducible results."),
    help_item("Run permutations", "Shuffles your data and re-runs inference many times to estimate
               a false discovery rate (FDR) for each edge weight, for a principled cutoff on the
               Network Diagnostics tab."),
    help_item("Number of permutations", "How many shuffled datasets to run (typically 100)."),
    help_item("Permute dimension", "\"col\" shuffles each sample's gene values independently;
               \"row\" shuffles each gene's values across samples."),
    help_item("Base seed", "Starting seed for the permutations -- permutation i uses base seed + i."),
    help_item("Target FDR", "The false discovery rate used to pick an edge-weight cutoff.
               Adjustable later on the Network Diagnostics tab without re-running anything."),
    help_item("Download full network", "Downloads the complete, unfiltered network as a
               tab-separated file."),
    shiny::h3("Network Diagnostics"),
    shiny::p("Shows edge weight and regulator out-degree distributions. If you ran permutations, also
              shows the FDR curve. Use \"Threshold network\" to filter by FDR or by a flat weight
              cutoff. \"Download all plots (PDF)\" saves everything as one file."),
    shiny::h3("Visualize"),
    shiny::p("Draws the network as currently thresholded on the Network Diagnostics tab. For very
              large networks or a publication figure, download the network and open it in Cytoscape
              instead.")
  )
}
