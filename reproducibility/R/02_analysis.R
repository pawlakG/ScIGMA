#' Run ScIGMA pipeline on the downloaded dataset
#' @param h5_path Path to the raw .h5 file
#' @param params Parsed YAML configuration list
#' @return A processed data object
run_scigma_pipeline <- function(h5_path, params, manual_gates) {
    message("Running analytical pipeline on: ", h5_path)

    # 1. High-Fidelity Data Ingestion
    message("Loading HDF5 matrix...")
    ScIGMA_data <- ScIGMA:::loadH5_HDF5_biocond(
        filepath = h5_path,
        sample_name = "Zenodo_AML_Test"
    )

    # 2. Rigorous Variant Filtering & Annotation (VEP/ClinVar)
    message("Filtering and annotating variants...")
    ScIGMA_data <- ScIGMA:::filter_and_annotate_variants(
        obj = ScIGMA_data,
        min_dp = params$variant_filtering$min_dp,
        min_gq = params$variant_filtering$min_gq,
        vaf_ref = params$variant_filtering$vaf_ref,
        vaf_hom = params$variant_filtering$vaf_hom,
        vaf_het = params$variant_filtering$vaf_het,
        min_cell_pt = params$variant_filtering$min_cell_pt,
    )

    # 3. Select DNA variant
    message("Select DNA variant...")
    selected_variants <- SummarizedExperiment::rowData(ScIGMA_data$mae[[
        "dna_variants"
    ]]) |>
        as.data.frame() |>
        dplyr::filter(cdna %in% c("c.700T>C", "c.672+1G>A", "c.638G>A"))
    ScIGMA_data$variants.filtered <- selected_variants

    # 4. Run COMPASS
    message("Run COMPASS...")
    ScIGMA_data <- ScIGMA:::run_compass_inference(
        obj = ScIGMA_data,
        chains = params$compass$chains,
        chain_length = params$compass$chain_length,
        patient_sex = params$compass$patient_sex
    )

    # 5. DNA Variant Heatmap Visualization
    message("Generating DNA variant heatmap...")
    ht_res <- ScIGMA:::generate_dna_variant_heatmap(
        obj = ScIGMA_data,
        selected_variants_df = ScIGMA_data$variants.filtered,
        heatmap_include_all_samples = FALSE,
        use_imputed = TRUE
    )

    message("ht_res computeted")
    active_clones <- ht_res$clones
    ScIGMA_data$cnv.active.clones <- active_clones

    # 6. Process protein data
    message("Processing protein data...")
    ScIGMA_data$mae <- ScIGMA:::sanitize_mae_strings(ScIGMA_data$mae)
    safe_rownames <- ScIGMA:::sanitize_protein_markers(rownames(ScIGMA_data$mae[[
        "proteins"
    ]]))
    rownames(ScIGMA_data$mae[["proteins"]]) <- safe_rownames
    ScIGMA_data$seurat_object <- ScIGMA:::protein_run_pca(ScIGMA_data)

    # 7. Filter CNV
    message("Filtering CNV...")
    filtered_data <- ScIGMA:::filter_cnv_profile(
        ScIGMA_data,
        active_clones,
        amp_completeness = params$cnv_filtering$cnv_ampCompleteness,
        amp_readDepth = params$cnv_filtering$cnv_ampReadDepth,
        amp_meanCellRead = params$cnv_filtering$cnv_meanCellReadDepth
    )
    ScIGMA_data$cnv_dp_filtered <- filtered_data
    ScIGMA_data$is_cnv_filtered <- TRUE

    ploidy_data <- ScIGMA:::process_cnv_to_clonal_profile(
        ScIGMA_data$cnv_dp_filtered,
        active_clones,
        diploid_ref = params$cnv_clonal_inference$ref_clone,
        exclude_clone = "small"
    )

    # Change R6 object
    ScIGMA_data$ploidy.mtx <- ploidy_data

    # 8. Immunophentype gating
    message("Immunophentype gating...")
    if (!is.null(manual_gates)) {
        message("Injecting manual immunophenotype gates...")
        ScIGMA_data$protein_gating_tree <- manual_gates
    }

    # 9. Run UMAP
    message("Running UMAP...")
    ScIGMA_data$seurat_object <- Seurat::RunUMAP(
        ScIGMA_data$seurat_object,
        features = params$protein_umap$umap_features,
        min.dist = params$protein_umap$umap_min_dist,
        n.neighbors = params$protein_umap$umap_n_neighbors,
        seed.use = 42
    )

    # 10. Unsupervised clustering
    message("Unsupervised clustering...")
    ScIGMA_data <- ScIGMA:::run_protein_clustering(
        obj = ScIGMA_data,
        resolution = params$unsupervised_custering$resolution,
        seed = 42
    )

    return(ScIGMA_data)
}
