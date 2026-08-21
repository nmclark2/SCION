version 1.0

## Runs the SCION network-inference + permutation-FDR-thresholding workflow on Terra,
## sharding the permutations across one scatter shard each.
##
## Uses the same `broadcptacdev/panoply_scion` image as PANOPLY's own
## panoply_scion.wdl (broadinstitute/PANOPLY, panoply_SCION branch) -- there is
## deliberately only one SCION Docker image, built from PANOPLY's hydrant task
## conventions, shared by both workflows. See that repo for the Dockerfile.

workflow scion_workflow {
  input {
    File target_data_file
    File reg_data_file
    File? target_genes_file
    File? reg_genes_file
    Boolean gene_list_header = true
    File? clustering_data_file
    String format = "csv"
    String clustering_method = "none"
    Float clustering_threshold = 0.5
    File? clusters_file
    Boolean connect_hubs = true
    Float weightthreshold = 0
    Boolean normalize = true
    Int num_cores = 1
    String engine = "randomForest"
    String ptm_sep = "."
    Int seed = 2020
    Int nb_trees = 10000

    Int n_permutations = 100
    String permute_dim = "col"
    Int base_seed = 0
    Float target_fdr = 0.05

    Int memory_gb = 16
    Int disk_gb = 50
    Int preemptible = 0
  }

  call run_real {
    input:
      target_data_file = target_data_file,
      reg_data_file = reg_data_file,
      target_genes_file = target_genes_file,
      reg_genes_file = reg_genes_file,
      gene_list_header = gene_list_header,
      clustering_data_file = clustering_data_file,
      format = format,
      clustering_method = clustering_method,
      clustering_threshold = clustering_threshold,
      clusters_file = clusters_file,
      connect_hubs = connect_hubs,
      weightthreshold = weightthreshold,
      normalize = normalize,
      num_cores = num_cores,
      engine = engine,
      ptm_sep = ptm_sep,
      seed = seed,
      nb_trees = nb_trees,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      preemptible = preemptible
  }

  scatter (i in range(n_permutations)) {
    call run_permutation {
      input:
        target_rds = run_real.target_rds,
        reg_rds = run_real.reg_rds,
        cluster_assignment_rds = run_real.cluster_assignment_rds,
        params_rds = run_real.params_rds,
        index = i + 1,
        base_seed = base_seed,
        permute_dim = permute_dim,
        preemptible = preemptible
    }
  }

  call aggregate_fdr {
    input:
      network_rds = run_real.network_rds,
      permutation_files = run_permutation.permutation_rds,
      target_fdr = target_fdr,
      preemptible = preemptible
  }

  output {
    File real_network_tsv = run_real.network_tsv
    File thresholded_network_tsv = aggregate_fdr.thresholded_network_tsv
    File fdr_curve_png = aggregate_fdr.fdr_curve_png
    File weight_comparison_png = aggregate_fdr.weight_comparison_png
    File? network_plot_png = aggregate_fdr.network_plot_png
    File fdr_result_rds = aggregate_fdr.fdr_result_rds
  }

  meta {
    author: "Natalie Clark"
    email: "nclark@broadinstitute.org"
  }
}

task run_real {
  input {
    File target_data_file
    File reg_data_file
    File? target_genes_file
    File? reg_genes_file
    Boolean gene_list_header
    File? clustering_data_file
    String format
    String clustering_method
    Float clustering_threshold
    File? clusters_file
    Boolean connect_hubs
    Float weightthreshold
    Boolean normalize
    Int num_cores
    String engine
    String ptm_sep
    Int seed
    Int nb_trees
    Int memory_gb
    Int disk_gb
    Int preemptible
  }

  command <<<
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/scion_run_real.R \
      --target_data_file ~{target_data_file} \
      --reg_data_file ~{reg_data_file} \
      ~{"--target_genes_file " + target_genes_file} \
      ~{"--reg_genes_file " + reg_genes_file} \
      --gene_list_header ~{gene_list_header} \
      ~{"--clustering_data_file " + clustering_data_file} \
      --format ~{format} \
      --clustering_method ~{clustering_method} \
      --clustering_threshold ~{clustering_threshold} \
      ~{"--clusters_file " + clusters_file} \
      --connect_hubs ~{connect_hubs} \
      --weightthreshold ~{weightthreshold} \
      --normalize ~{normalize} \
      --num_cores ~{num_cores} \
      --engine ~{engine} \
      --ptm_sep '~{ptm_sep}' \
      --seed ~{seed} \
      --nb_trees ~{nb_trees} \
      --out_dir out
  >>>

  output {
    File network_tsv = "out/network.tsv"
    File network_rds = "out/network.rds"
    File target_rds = "out/target.rds"
    File reg_rds = "out/reg.rds"
    File cluster_assignment_rds = "out/cluster_assignment.rds"
    File params_rds = "out/params.rds"
  }

  runtime {
    docker: "broadcptacdev/panoply_scion:latest"
    memory: memory_gb + "GB"
    disks: "local-disk " + disk_gb + " SSD"
    cpu: num_cores
    preemptible: preemptible
  }
}

task run_permutation {
  input {
    File target_rds
    File reg_rds
    File cluster_assignment_rds
    File params_rds
    Int index
    Int base_seed
    String permute_dim
    Int preemptible
  }

  command <<<
    set -euo pipefail
    Rscript /prot/proteomics/Projects/PGDAC/src/scion_run_permutation.R \
      --target_rds ~{target_rds} \
      --reg_rds ~{reg_rds} \
      --cluster_assignment_rds ~{cluster_assignment_rds} \
      --params_rds ~{params_rds} \
      --index ~{index} \
      --base_seed ~{base_seed} \
      --permute_dim ~{permute_dim} \
      --num_cores 1 \
      --out_dir out
  >>>

  output {
    File permutation_rds = glob("out/permutation_*.rds")[0]
  }

  runtime {
    docker: "broadcptacdev/panoply_scion:latest"
    memory: "8GB"
    disks: "local-disk 20 SSD"
    cpu: 1
    preemptible: preemptible
  }
}

task aggregate_fdr {
  input {
    File network_rds
    Array[File] permutation_files
    Float target_fdr
    Int preemptible
  }

  command <<<
    set -euo pipefail
    mkdir -p permutations
    for f in ~{sep=" " permutation_files}; do
      cp "$f" permutations/
    done

    Rscript /prot/proteomics/Projects/PGDAC/src/scion_aggregate_fdr.R \
      --network_rds ~{network_rds} \
      --permutation_dir permutations \
      --target_fdr ~{target_fdr} \
      --out_dir out
  >>>

  output {
    File thresholded_network_tsv = "out/thresholded_network.tsv"
    File fdr_curve_png = "out/fdr_curve.png"
    File weight_comparison_png = "out/weight_comparison.png"
    File? network_plot_png = "out/network_plot.png"
    File fdr_result_rds = "out/fdr_result.rds"
  }

  runtime {
    docker: "broadcptacdev/panoply_scion:latest"
    memory: "8GB"
    disks: "local-disk 20 SSD"
    cpu: 1
    preemptible: preemptible
  }
}
