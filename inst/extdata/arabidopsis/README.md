# Arabidopsis example data

Used by the Shiny app's "Load example data" button and available via
`system.file("extdata", "arabidopsis", package = "SCION")`.

- `target_mat_RNA.csv` / `target_list_RNA.csv` -- RNA-seq expression, target genes.
- `reg_mat_protein.csv` / `reg_list_protein.csv` -- protein expression, regulators (TFs).
- `reg_mat_phospho.csv` / `reg_list_phospho.csv` -- phosphoproteome expression, regulators
  (TFs with PTM sites, dot notation e.g. `AT1G01260.p395`). The app's "Example regulator
  type" selector picks between this and the protein variant above; both pair with the same
  RNA target matrix/list.

From the JA/brassinosteroid dataset published in Clark, N.M., Nolan, T.M., Wang, P. et al.
*Integrated omics networks reveal the temporal signaling events of brassinosteroid response
in Arabidopsis*. Nat Commun 12, 5858 (2021).

**Provenance note**: these are the already-processed `Network_Files/` outputs from the
repo's original `TEST.zip` (produced by `utilities/create_data_tables.R` from the raw
`JA_FPKM.csv`/proteomics files). Copied here as-is to unblock the example-data button;
a fully reproducible `data-raw/` script regenerating these from the raw inputs -- the
"revised test-data plan" from the restructuring effort -- is still a follow-up, not done
here.
