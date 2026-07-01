library(targets)
# library(ScIGMA) # Load the package when it is available

# Source custom functions
source("R/01_download.R")
source("R/02_analysis.R")
source("R/03_figures.R")

# Pipeline definition
list(
  # ----------------------------------------------------------------------------
  # 0. CONFIGURATION
  # ----------------------------------------------------------------------------
  tar_target(
    config_file,
    "config.yml",
    format = "file"
  ),
  tar_target(
    params,
    config::get(file = config_file)
  ),
  tar_target(
    manual_gates_file,
    "data/manual_gates.rds",
    format = "file"
  ),
  tar_target(
    manual_gates,
    readRDS(manual_gates_file)
  ),

  # ----------------------------------------------------------------------------
  # 1. DATA ACQUISITION & PROCESSING
  # ----------------------------------------------------------------------------
  tar_target(
    raw_h5_file,
    download_zenodo_h5(
      url = params$zenodo$dataset_url,
      dest = "data/test_dataset.h5"
    ),
    format = "file"
  ),

  tar_target(
    ScIGMA_obj,
    run_scigma_pipeline(raw_h5_file, params, manual_gates)
  ),

  # ----------------------------------------------------------------------------
  # 2. SNV FIGURES (DNA)
  # ----------------------------------------------------------------------------
  tar_target(
    fig_snv_heatmap,
    plot_fig_snv_heatmap(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_snv_phylogeny,
    plot_fig_snv_phylogeny(ScIGMA_obj),
    format = "file"
  ),

  # ----------------------------------------------------------------------------
  # 3. CNV FIGURES
  # ----------------------------------------------------------------------------
  tar_target(
    fig_cnv_heatmap,
    plot_fig_cnv_heatmap(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_cnv_lineplot,
    plot_fig_cnv_lineplot(ScIGMA_obj),
    format = "file"
  ),

  # ----------------------------------------------------------------------------
  # 4. PROTEIN FIGURES
  # ----------------------------------------------------------------------------
  tar_target(
    fig_protein_ridgeplot,
    plot_fig_protein_ridgeplot(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_barplot,
    plot_fig_protein_barplot(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_biplot_gating,
    plot_fig_protein_biplot_gating(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_umap,
    plot_fig_protein_umap(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_umap_markers,
    plot_fig_protein_umap_markers(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_umap_gating,
    plot_fig_protein_umap_gating(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_protein_umap_clustering,
    plot_fig_protein_umap_clustering(ScIGMA_obj),
    format = "file"
  ),

  # ----------------------------------------------------------------------------
  # 5. MULTIOMICS INTEGRATION FIGURES
  # ----------------------------------------------------------------------------
  tar_target(
    fig_multi_umap_dna_clones,
    plot_fig_multi_umap_dna_clones(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_multi_umap_dna_genotype,
    plot_fig_multi_umap_dna_genotype(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_multi_barplot_clones_per_cluster,
    plot_fig_multi_barplot_clones_per_cluster(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_multi_barplot_variant_genotype,
    plot_fig_multi_barplot_variant_genotype(ScIGMA_obj),
    format = "file"
  ),
  tar_target(
    fig_multi_barplot_gating_genotype,
    plot_fig_multi_barplot_gating_genotype(ScIGMA_obj),
    format = "file"
  )
)
