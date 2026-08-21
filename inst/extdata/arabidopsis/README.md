# Arabidopsis example data

Used by the Shiny app's "Load example data" button and available via
`system.file("extdata", "arabidopsis", package = "SCION")`.

- `target_mat_RNA.csv` / `target_list_RNA.csv` -- RNA-seq expression, target genes.
- `reg_mat_protein.csv` / `reg_list_protein.csv` -- protein expression, regulators (TFs).
- `reg_mat_phospho.csv` / `reg_list_phospho.csv` -- phosphoproteome expression, regulators
  (TFs with PTM sites, dot notation e.g. `AT1G01260.p395`). The app's "Example regulator
  type" selector picks between this and the protein variant above; both pair with the same
  RNA target matrix/list.
- `cluster_mat_protein.csv` / `cluster_mat_phospho.csv` -- pre-computed clustering matrices,
  one per regulator type, used when the "Load example data" button's default temporal (DTW)
  clustering step runs.

From the jasmonic acid (JA) response dataset published in Zander, M., Lewsey, M.G., Clark,
N.M. et al. *Integrated multi-omics framework of the plant response to jasmonic acid*. Nature
Plants 6, 290-302 (2020). Also used as the worked example in the book chapter tutorial: Clark,
N.M., Hurgobin, B., Kelley, D.R., Lewsey, M.G., Walley, J.W. (2023). *A Practical Guide to
Inferring Multi-Omics Networks in Plant Systems*. Methods in Molecular Biology, vol 2698.

**Provenance note**: these are the already-processed `Network_Files/` outputs from the
repo's original `TEST.zip` (produced by `utilities/create_data_tables.R` from the raw
`JA_FPKM.csv`/proteomics files). Copied here as-is to unblock the example-data button;
a fully reproducible `data-raw/` script regenerating these from the raw inputs -- the
"revised test-data plan" from the restructuring effort -- is still a follow-up, not done
here.
