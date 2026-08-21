# Running SCION on Terra

`scion.wdl` runs the full run -> permute -> FDR-threshold workflow, sharding
permutations one-per-scatter-shard instead of running them all in one task.

There is deliberately only one SCION Docker image. It's built from
[broadinstitute/PANOPLY](https://github.com/broadinstitute/PANOPLY)'s
`panoply_SCION` branch (`hydrant/tasks/panoply_scion/panoply_scion/Dockerfile`),
following PANOPLY's own hydrant task conventions, and published as
`broadcptacdev/panoply_scion:latest` -- the same image PANOPLY's own
`panoply_scion.wdl` uses. This WDL hardcodes that image directly; there's no
separate build for standalone use.

## Run the workflow

```json
{
  "scion_workflow.target_data_file": "gs://.../target_mat.csv",
  "scion_workflow.reg_data_file": "gs://.../reg_mat.csv",
  "scion_workflow.n_permutations": 100
}
```

See `scion.wdl`'s `workflow` input block for every available parameter --
they mirror `run_scion()`'s arguments directly (format, clustering_method,
weightthreshold, engine, etc.), plus `n_permutations`/`permute_dim`/
`base_seed`/`target_fdr` for the permutation/FDR step.

## Outputs

- `real_network_tsv` -- the real (unpermuted) network, Cytoscape-importable.
- `thresholded_network_tsv` -- the real network filtered to the FDR-selected
  weight cutoff.
- `fdr_curve_png`, `weight_comparison_png` -- diagnostic plots from
  `plot_fdr_curve()`.
- `network_plot_png` -- a static `plot_network()` render of the thresholded
  network (omitted if the threshold left zero edges).
- `fdr_result_rds` -- the full `compute_fdr_threshold()` result object, for
  programmatic follow-up in R.

## Updating the image

The image is rebuilt/republished from the PANOPLY repo, not from here. If
`broadcptacdev/panoply_scion:latest` moves to a different SCION package ref,
that's a PANOPLY-side change (its Dockerfile's `scion_ref` build-arg) --
nothing in this repo needs to change unless the driver scripts' CLI
interface itself changes (`scripts/scion_run_real.R`,
`scion_run_permutation.R`, `scion_aggregate_fdr.R`), which PANOPLY vendors a
copy of.
