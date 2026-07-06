# ========================================================================= #
# Pipeline d'Exploration ScIGMA (h5_deepDive.R)
# Reproduction de l'analyse (reproducibility) + Intégration Normalisation DSB
# ========================================================================= #

devtools::load_all()

# Chemins de fichiers (ajuster selon le besoin)
directory <- "../inputs/tapestriDatasets/4-cell-lines-AML-multiomics/4-cell-lines-AML-multiomics.dna+protein.h5"

# 1. High-Fidelity Data Ingestion
message("Loading HDF5 matrix...")
ScIGMA_data <- loadH5_HDF5_biocond(
    filepath = directory,
    sample_name = "aml_4_lines",
    omic_type = "DNA+protein"
)

# 2. Rigorous Variant Filtering & Annotation (VEP/ClinVar)
message("Filtering and annotating variants...")
ScIGMA_data <- filter_and_annotate_variants(
    obj = ScIGMA_data,
    min_dp = 10,
    min_gq = 30,
    vaf_ref = 5,
    vaf_hom = 95,
    vaf_het = 20,
    min_cell_pt = 10,
    min_mut_cell_pt = 1
)

# 3. Select DNA variant
message("Select DNA variant...")
# Extraction de l'annotation native sans altérer la colonne variant_id d'origine
# qui contient le préfixe du gène (ex: NRAS:chr1...).
sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[[
    "dna_variants"
]]) |>
    as.data.frame() |>
    dplyr::select(
        variant_id,
        gene,
        transcript_id,
        protein,
        cdna,
        variant_type,
        gene_function,
        clinvar,
        cell_proportion
    ) |>
    dplyr::arrange(desc(cell_proportion), clinvar)

sorted_annotation <- sorted_annotation %>% filter(grepl("pathogenic", clinvar))

# Sélection des 20 variants les plus prévalents
target_variants <- sorted_annotation$variant_id[
    1:min(20, nrow(sorted_annotation))
]
selected_df <- sorted_annotation |>
    dplyr::filter(variant_id %in% target_variants)
# On s'assure que les rownames matchent le format d'origine de la matrice (l'index)
rownames(
    selected_df
) <- rownames(SummarizedExperiment::rowData(ScIGMA_data$mae[[
    "dna_variants"
]]))[match(
    selected_df$variant_id,
    SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]])$variant_id
)]

ScIGMA_data$variants.filtered <- selected_df

# 4. Run COMPASS
message("Run COMPASS...")
ScIGMA_data <- run_compass_inference(
    obj = ScIGMA_data,
    chains = 4,
    chain_length = 300,
    patient_sex = "M"
)

# 5. DNA Variant Heatmap Visualization
message("Generating DNA variant heatmap...")
ht_res <- generate_dna_variant_heatmap(
    obj = ScIGMA_data,
    selected_variants_df = ScIGMA_data$variants.filtered,
    heatmap_include_all_samples = FALSE,
    use_imputed = TRUE
)

# Rendu de la heatmap pour vérification visuelle (Optionnel)
ht <- ComplexHeatmap::draw(ht_res$heatmap)

active_clones <- ht_res$clones
ScIGMA_data$cnv.active.clones <- active_clones
ScIGMA_data$dna.clones <- active_clones

# 6. Process protein data (AVEC LE COMMUTATEUR DSB)
message("Processing protein data...")
ScIGMA_data$mae <- sanitize_mae_strings(ScIGMA_data$mae)
safe_rownames <- sanitize_protein_markers(rownames(ScIGMA_data$mae[[
    "proteins"
]]))
rownames(ScIGMA_data$mae[["proteins"]]) <- safe_rownames

# --- Choix de la Méthode de Normalisation ---
normalization_method <- "dsb" # Options: "dsb" ou "linear_regression"
message(sprintf("Running Normalization using: %s", normalization_method))

if (normalization_method == "dsb") {
    # Normalisation par soustraction du bruit de fond (Gouttes vides)
    ScIGMA_data <- normalize_protein_dsb(
        ScIGMA_data,
        min_bg_counts = 10,
        max_bg_counts = 500
    )
    norm_mat <- SummarizedExperiment::assay(
        ScIGMA_data$mae[["proteins"]],
        "normalized"
    )
} else {
    # Normalisation standard CLR / régression linéaire
    prot_counts <- SummarizedExperiment::assay(
        ScIGMA_data$mae[["proteins"]],
        "counts"
    )
    norm_mat <- normalize_linear_regression(
        as.matrix(prot_counts),
        jitter = 0.5
    )
    SummarizedExperiment::assay(
        ScIGMA_data$mae[["proteins"]],
        "normalized"
    ) <- norm_mat
    S4Vectors::metadata(
        ScIGMA_data$mae
    )$protein_normalize_method <- "Linear Regression"
}

ScIGMA_data$protein.filtered <- TRUE

# Création de l'objet Seurat
min_features_threshold <- floor(sqrt(nrow(norm_mat)))
seurat_obj <- Seurat::CreateSeuratObject(
    counts = norm_mat,
    project = "ScIGMA_proteins",
    min.cells = 3,
    min.features = min_features_threshold
)

# Restauration des codes-barres (Seurat remplace '_' par '-')
original_cells <- colnames(norm_mat)
seurat_cells <- gsub("_", "-", original_cells)
cell_mapping <- setNames(original_cells, seurat_cells)
current_seurat_cells <- colnames(seurat_obj)
seurat_obj <- Seurat::RenameCells(
    seurat_obj,
    new.names = unname(cell_mapping[current_seurat_cells])
)

seurat_obj <- Seurat::SetAssayData(
    seurat_obj,
    layer = "data",
    new.data = norm_mat
)

# Identification des protéines les plus variables
seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = "vst",
    nfeatures = nrow(seurat_obj)
)

# Fix ScaleData : Ne pas centrer/scaler si les valeurs sont déjà normalisées par dsb
if (normalization_method == "dsb") {
    seurat_obj <- Seurat::ScaleData(
        seurat_obj,
        features = rownames(seurat_obj),
        do.scale = FALSE,
        do.center = FALSE
    )
} else {
    seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))
}

# ACP Protéique
max_npcs <- min(nrow(seurat_obj) - 1, ncol(seurat_obj) - 1, 50)
seurat_obj <- Seurat::RunPCA(
    seurat_obj,
    features = Seurat::VariableFeatures(object = seurat_obj),
    npcs = max_npcs,
    verbose = FALSE,
    seed.use = 42
)

ScIGMA_data$seurat_object <- seurat_obj

# 7. Filter CNV
message("Filtering CNV...")
filtered_data <- filter_cnv_profile(
    ScIGMA_data,
    active_clones,
    amp_completeness = 50,
    amp_readDepth = 10,
    amp_meanCellRead = 10
)
ScIGMA_data$cnv_dp_filtered <- filtered_data
ScIGMA_data$is_cnv_filtered <- TRUE

# Inférence de ploïdie
ScIGMA_data$ploidy.mtx <- process_cnv_to_clonal_profile(
    ScIGMA_data$cnv_dp_filtered,
    active_clones,
    diploid_ref = "clone_01",
    exclude_clone = "small"
)

# 8. Run UMAP
message("Running UMAP...")
ScIGMA_data$seurat_object <- Seurat::RunUMAP(
    ScIGMA_data$seurat_object,
    dims = 1:(nrow(ScIGMA_data$seurat_object) - 2),
    min.dist = 0.15,
    n.neighbors = 30,
    seed.use = 42
)

# 9. Unsupervised clustering
message("Unsupervised clustering...")
ScIGMA_data <- run_protein_clustering(
    obj = ScIGMA_data,
    resolution = 0.15,
    seed = 42
)

message("==== Deep Dive Pipeline Complete ! ====")
