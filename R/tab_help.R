################################################################################
# Module: Help (static content)
################################################################################

#' @keywords internal
helpTabUI <- function(id = "help") {
  shiny::tagList(
    shiny::h3("Run network (sidebar)"),
    shiny::strong("Load example data"),
    shiny::p("Runs a real Arabidopsis network (from Zander, Lewsey, Clark et al., Nature Plants
              6, 290-302 (2020)) with the parameters used in the README/tutorial -- a quick way
              to see the app working end-to-end before uploading your own data. Choose \"Example
              regulator type\" for a protein-level network (transcription factors) or a
              phosphoproteome-level one (TFs with PTM sites); both pair with the same RNA
              target data and their own bundled clustering matrix, used with temporal (DTW)
              clustering (switched on for this example, since a pre-computed clustering file
              is included)."),
    shiny::strong("Target / Regulator matrix"),
    shiny::p("Data matrices for your gene targets (e.g. transcript data) and regulators (e.g. protein
              data). Rows are genes, columns are samples. CSV or GCT format. If your regulators contain
              additional information such as a PTM site, denote it with the PTM site separator below
              (e.g. SOX2.S35 or SOX2_S35)."),
    shiny::strong("Target / Regulator gene list"),
    shiny::p("Optional list of genes to use as targets/regulators (first column = gene names). Leave
              blank to use every gene present in the corresponding matrix. Check \"Gene lists have a
              header row\" only if the file has a header -- a plain one-gene-per-line file (no header)
              needs it unchecked, or its first gene will be silently dropped."),
    shiny::strong("Clustering"),
    shiny::p("Optionally cluster genes before network inference: temporal (DTW), non-temporal (ICA or
              k-means), or upload your own pre-computed clusters. If you don't supply a clustering
              matrix, it defaults to your target and regulator data combined."),
    shiny::strong("Edge weight cutoff / Normalize edge weights"),
    shiny::p("A manual cutoff for trimming low-confidence edges (0 keeps everything). Normalizing
              rescales edge weights to [0, 1] before applying the cutoff. For a principled,
              FDR-controlled cutoff instead of a manual one, leave this at 0 and check
              \"Run permutations\" below instead."),
    shiny::strong("Random forest engine"),
    shiny::p("randomForest (default) matches previously published results. ranger is a much faster
              alternative, at the cost of no longer being numerically comparable to randomForest-based
              runs."),
    shiny::strong("Run permutations"),
    shiny::p("Runs repeated permutations of your target/regulator data (typically 100) to build a null
              distribution of edge weights, then picks the largest edge-weight cutoff for which the
              false discovery rate stays below your target (default 0.05). Results show up in the
              Network Diagnostics tab alongside the weight and out-degree distributions."),
    shiny::h3("Network Diagnostics"),
    shiny::p("Always shows the real network's edge weight distribution and regulator out-degree
              distribution. Once you've run permutations, the FDR curve and the real-vs-permuted
              weight distribution appear here too, to sanity-check the FDR-based cutoff."),
    shiny::h3("Visualize"),
    shiny::p("Always shows the thresholded network from the Network Diagnostics tab (FDR-based,
              a manual weight cutoff, or the full real network if neither has been applied there) --
              adjust the cutoff on that tab rather than switching views here. Interactive mode uses
              visNetwork; for large networks or a publication figure, download the network file and
              import it into Cytoscape instead.")
  )
}
