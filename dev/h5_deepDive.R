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

# 1. Chargement Out-of-Core HDF5
ScIGMA_data <- loadH5_HDF5_biocond(
    filepath = directory, # Assure-toi que la variable directory est bien définie
    sample_name = "aml_4_lines",
    omic_type = "DNA+protein"
)



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

if (length(target_variants) == 0) stop("Fatal: No pathogenic variants found for COMPASS.")

# 4. Extraction des Matrices depuis le MAE
compass_inputs <- build_compass_matrices(
    obj = ScIGMA_data,
    selected_variants = target_variants
)

# EXTRACTION GT: Le nouveau wrapper C++ requiert la matrice de Génotype
# On extrait, transpose (Cellules x Variants) et on annule les 3 (Missing)
mat_ref <- as.matrix(t(compass_inputs$M_ref))
storage.mode(mat_ref) <- "integer"

mat_alt <- as.matrix(t(compass_inputs$M_alt))
storage.mode(mat_alt) <- "integer"

mat_cna <- as.matrix(t(compass_inputs$C))
storage.mode(mat_cna) <- "integer"

# EXTRACTION GT : Suppression du t() pour garder nativement [Variants x Cells]
mat_gt <- as.matrix(SummarizedExperiment::assay(ScIGMA_data$mae[["dna_variants"]], "gt")[target_variants, , drop = FALSE])
mat_gt[mat_gt == 3L] <- NA
storage.mode(mat_gt) <- "integer"

variant_matrices <- list(
    REF = mat_ref,
    ALT = mat_alt,
    GT  = mat_gt
)


# 1. Extraction des métadonnées SNV
snv_rowData <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]))
snv_sub <- snv_rowData[target_variants, ]

vec_locus_names <- snv_sub$gene        # Ex: "NPM1"
vec_locus_chrom <- snv_sub$chrom         # Ex: "5" (Le C++ gère le préfixe "chr" en interne)

# 2. Extraction des métadonnées CNA
cna_rowData <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))

# Agrégation Spatiale (Garantir le format CHR_GENE pour le C++)
vec_region_names <- paste0(cna_rowData$chrom, "_", sapply(cna_rowData$dna_id, \(x) strsplit(x, "_")[[1]][3]))
vec_region_names <- unique(vec_region_names) # Alignés avec les colonnes de ta matrice C_mat

# NEW : Extraction robuste des chromosomes, strictement alignée sur vec_region_names
vec_region_chrom <- sapply(vec_region_names, \(x) strsplit(x, "_")[[1]][1], USE.NAMES = FALSE)
vec_region_chrom <- sub("^chr", "", vec_region_chrom, ignore.case = TRUE) # Nettoyage de sécurité

prefix_out <- "results/compass_output/dut_env"

# Handle if dim(CNV) ≠ dim(SNV)
if (ncol(variant_matrices$REF) != ncol(mat_cna)){
    use_cna <- FALSE
} else {
    use_cna <- TRUE
}

# 3. L'appel final blindé
run_compass_mcmc(
    variant_matrices   = variant_matrices,
    locus_regions      = compass_inputs$locus_regions,
    region_matrix      = mat_cna,
    output_prefix      = prefix_out,
    locus_names        = vec_locus_names,
    locus_chromosomes  = vec_locus_chrom,
    region_names       = vec_region_names,
    region_chromosomes = vec_region_chrom, # <-- NEW INJECTION
    chains             = 4L,
    chain_length       = 300L,
    patient_sex        = "female",
    use_cna            = use_cna
)









ScIGMA_data$seurat_object <- RunUMAP(ScIGMA_data$seurat_object,
                                     dims = 1:(nrow(ScIGMA_data$seurat_object)-2),
                                     min.dist = 0.15,
                                     n.neighbors = 30,
                                     future.seed=TRUE)
