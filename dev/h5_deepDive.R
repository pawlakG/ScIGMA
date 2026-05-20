devtools::load_all()

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/Jas_Rec/Jas_rec.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/2_PBMC_Mix_KG-1_Spike_In/Sample_Data_Set_2_PBMC_Mix_KG-1_Spike_In.labeled.dna+protein.h5"
directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"
directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"
directory <- "../inputs/bPodvinDatasets/dew_diag.h5"

rhdf5::h5ls(directory, all = TRUE, recursive = TRUE)

h5f <- H5Fopen(directory, flags = "H5F_ACC_RDONLY")

h5f$assays$dna_variants$layers$NGT |> table()

h5f$assays$dna_read_counts$ca

dim(h5f$assays$protein_read_counts$layers$read_counts)

h5f$assays$protein_read_counts$ra$sample_name


h5closeAll()


# [ NODE_ACCESS : Test MAE pipeline ]
# ----------------------------------------------------- _
# ========================================================================= #
# Pipeline d'Inférence Phylogénétique COMPASS (In-Memory Rcpp)
# ========================================================================= #

devtools::load_all()

directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"
# directory <- "../inputs/bPodvinDatasets/Dut_Ev/Dut_Ev_Rec.dna+protein.h5"
# directory <- "../inputs/bPodvinDatasets/LMC/M35.dna+protein.h5"
# directory <- "../inputs/tapestriDatasets/KG-1-Raji-50-50-Myeloid/KG-1-Raji-50-50-Myeloid.dna.h5"

# 1. Chargement Out-of-Core HDF5
ScIGMA_data <- loadH5_HDF5_biocond(
    filepath = directory, # Assure-toi que la variable directory est bien définie
    sample_name = "aml_4_lines",
    omic_type = "DNA+protein"
    # omic_type = "DNA"
)

ScIGMA_data$mae <- sanitize_mae_strings(ScIGMA_data$mae)

print(head(rownames(ScIGMA_data$mae[["dna_variants"]])))
print(length(rownames(ScIGMA_data$mae[["dna_variants"]])))

# 2. Filtrage & Annotation (Un seul appel robuste)
tryCatch({
    # Assure-toi que cfg$paths est disponible dans l'environnement
    ScIGMA_data <- filter_and_annotate_variants(
        obj = ScIGMA_data,
        paths = cfg$paths,
        min_dp = 10,
        min_gq = 30,
        min_cell_pt = 10,
        min_mut_cell_pt = 1
    )
}, error = function(e) {
    stop(sprintf("Pipeline failed during filtering/annotation: %s", e$message))
})
# saveRDS(ScIGMA_data, "dev/ScIGMA_data_after_filter_and_annotate_variants.rds")
# ScIGMA_data <- readRDS("dev/ScIGMA_data_after_filter_and_annotate_variants.rds")
# 3. Sélection des cibles (CNA / SNV)
annotation_df <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]))



# Filtrage sur les variants Pathogènes purs
pathogenic_variants <- annotation_df[grepl("Pathogenic", annotation_df$clinvar, ignore.case = TRUE), ]
pathogenic_variants <- pathogenic_variants[order(pathogenic_variants$cell_proportion, decreasing = TRUE), ]

if (nrow(pathogenic_variants) == 0){
    pathogenic_variants <- annotation_df |> arrange(desc(clinvar))
    pathogenic_variants <- pathogenic_variants[1:10,]
}

target_variants <- rownames(pathogenic_variants)[1:min(20, nrow(pathogenic_variants))]

# >> CNV _

sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
    as.data.frame() |>
    dplyr::select(variant_id, gene, transcript_id, protein, cdna, variant_type, gene_function, impact, clinvar, cell_proportion) |>
    dplyr::arrange(desc(cell_proportion), desc(impact))

# Extraction sécurisée des variants sélectionnés
selected_df <- sorted_annotation[target_variants, , drop = FALSE]

# Mise à jour de l'objet global (au cas où d'autres modules l'utilisent)
ScIGMA_data$variants.filtered <- selected_df

# Génération de la Heatmap
ht_res <- generate_dna_variant_heatmap(
    obj = ScIGMA_data,
    selected_variants_df = selected_df,
    heatmap_include_all_samples = TRUE
)

ht <- ComplexHeatmap::draw(ht_res$heatmap)
print("New heatmap rendered")



# Gestion des clones
if (is.null(ScIGMA_data$dna_clones_renamed)) {
    ScIGMA_data$dna.clones <- ht_res$clones
}

filtered_data <- filter_cnv_profile(ScIGMA_data,
                                    ScIGMA_data$dna.clones,
                                    amp_completeness = 50,
                                    amp_readDepth = 10,
                                    amp_meanCellRead = 10)




ScIGMA_data <- infer_clonal_architecture(scigma_data = ScIGMA_data,
                                         target_variants = target_variants,
                                         chain_length = 200)


ScIGMA_data$seurat_object <- protein_run_pca(ScIGMA_data)



ScIGMA_data$seurat_object <- RunUMAP(ScIGMA_data$seurat_object,
                                     dims = 1:(nrow(ScIGMA_data$seurat_object)-2),
                                     min.dist = 0.15,
                                     n.neighbors = 30,
                                     future.seed=TRUE)

ScIGMA_data$seurat_object <- FindNeighbors(
    ScIGMA_data$seurat_object,
    features = rownames(ScIGMA_data$seurat_object),
    dims = 1:ncol(ScIGMA_data$seurat_object@reductions$pca@cell.embeddings)
)

ScIGMA_data$seurat_object <- FindClusters(
    ScIGMA_data$seurat_object,
    resolution = 0.15
)

# Get dna clones from compass mtx
extracted_variant_genotypes <- extract_variant_genotypes(mae_data = ScIGMA_data$mae,
                                                         variant_id = rownames(ScIGMA_data$variants.filtered)[1], use_compass = T)

extracted_variant_genotypes_df_plot <- extracted_variant_genotypes |> t() |> as.data.frame() |> select(1) |> rownames_to_column("cell_barcode")
ref_umap <- as.data.frame(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings) |> rownames_to_column("cell_barcode")
extracted_variant_genotypes_df_plot <- merge(extracted_variant_genotypes_df_plot, ref_umap) |> column_to_rownames("cell_barcode")
