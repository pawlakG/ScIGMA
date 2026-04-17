devtools::load_all()
library(pbapply)
library(data.table)
# [ NODE_ACCESS : Process data from COMPASS article ]
# ----------------------------------------------------- _
data_folder <- "../COMPASS/data/preprocessed_data_AML_Morita2020/"
data_ids <- list.files(data_folder, pattern = "^AML")
data_ids <- data_ids[!grepl("PDX", data_ids)] # remove PDX
data_ids <- data_ids[!grepl("AML-88", data_ids)] # remove PDX
data_ids_samples <- setNames(data_ids, nm = sapply(data_ids, \(x){
    strsplit(x = x, "_")[[1]][1]
}))

data_id_processedFiles <- pbsapply(unique(names(data_ids_samples)) , \(x){
    tmp_files <- data_ids_samples[grepl(x, names(data_ids_samples))]
    print(tmp_files)
    tmp_variant <- paste0(data_folder, tmp_files[grepl("variant", tmp_files)])
    tmp_regions <- paste0(data_folder, tmp_files[grepl("regions", tmp_files)])
    prepare_compass_from_csv(variants_file = tmp_variant, regions_file = tmp_regions)
}, USE.NAMES = TRUE, simplify = FALSE)

data_id_processedFiles_compass <- data_id_processedFiles[c("AML-59-001", "AML-99-001", "AML-83-002", "AML-101-001")]


# UPDATED
# File: scripts/run_compass_multisample_grid.R

library(parallel)

# 1. Définition de la grille de paramètres MCMC
param_grid <- list(
    list(chains = 2L, length = 5L),
    list(chains = 4L, length = 10L),
    list(chains = 8L, length = 50L),
    list(chains = 10L, length = 100L)
    # list(chains = 14L, length = 200L),
    # list(chains = 16L, length = 400L),
    # list(chains = 18L, length = 600L),
    # list(chains = 4L, length = 5000L)
)

# 2. Construction du plan d'exécution (4 datasets x 5 configs = 20 jobs)
# On suppose que data_id_processedFiles_compass est déjà chargé en mémoire
dataset_names <- names(data_id_processedFiles_compass)

job_plan <- list()
for (ds_name in dataset_names) {
    for (params in param_grid) {
        job_plan[[length(job_plan) + 1]] <- list(
            dataset = ds_name,
            chains  = params$chains,
            length  = params$length
        )
    }
}

# 3. Sécurité du répertoire racine
base_output_dir <- "results/compass_output"
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

message(sprintf("Déploiement de %d jobs MCMC au total (Forking UNIX)...", length(job_plan)))
t0 <- Sys.time()

# 4. Exécution asynchrone sécurisée
# On bride volontairement mc.cores à 4 pour ne pas crasher le CPU M4 Pro
# si le backend C++ utilise de l'OpenMP sous le capot.
results <- parallel::mclapply(
    X = seq_along(job_plan),
    FUN = function(i) {
        job <- job_plan[[i]]

        # Isolation des I/O : Un sous-dossier par échantillon
        sample_dir <- file.path(base_output_dir, job$dataset)
        dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

        # Identification stricte de la sortie
        run_name <- sprintf("%s_c%d_l%d", job$dataset, job$chains, job$length)
        prefix_out <- file.path(sample_dir, run_name)

        # Extraction des pointeurs mémoires depuis la liste consolidée
        compass_inputs <- data_id_processedFiles_compass[[job$dataset]]

        message(sprintf("[Worker %02d/%02d] Démarrage -> %s", i, length(job_plan), run_name))

        # Barrière tryCatch pour isoler les Segfaults potentiels d'un dataset corrompu
        run_status <- tryCatch({
            # run_compass_mcmc(
            #     variant_matrices = compass_inputs$variant_matrices,
            #     locus_regions    = compass_inputs$locus_regions,
            #     region_matrix    = compass_inputs$region_matrix,
            #     output_prefix    = prefix_out,
            #     chains           = job$chains,
            #     chain_length     = job$length,
            #     # Attention : Assure-toi que "female" s'applique à tous les patients de la liste.
            #     # Sinon, il faudra l'injecter dynamiquement depuis les métadonnées.
            #     patient_sex      = "female"
            # )

            vec_locus_names <- rownames(compass_inputs$variant_matrices$REF)
            vec_locus_chrom <- as.character(compass_inputs$locus_regions)

            # 2. Extraction des métadonnées CNA
            # cna_rowData <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
            # Agrégation Spatiale (Garantir le format CHR_GENE pour le C++)
            # vec_region_names <- paste0(cna_rowData$chrom, "_", cna_rowData$gene) # Ex: "17_TP53"
            vec_region_names <- rownames(compass_inputs$region_matrix)
            vec_region_names <- unique(vec_region_names) # Alignés avec les colonnes de ta matrice C_mat




            run_compass_mcmc(
                variant_matrices  = compass_inputs$variant_matrices,
                locus_regions     = compass_inputs$locus_regions,
                region_matrix     = compass_inputs$region_matrix,
                output_prefix     = prefix_out,
                locus_names       = vec_locus_names,
                locus_chromosomes = vec_locus_chrom,
                region_names      = vec_region_names,
                chains            = job$chains,
                chain_length      = job$length,
                patient_sex       = "female"
            )
        }, error = function(e) {
            warning(sprintf("\n[Worker %02d] Échec critique sur %s: %s", i, run_name, e$message))
            return(FALSE)
        })

        return(data.frame(
            sample = job$dataset,
            run_id = run_name,
            chains = job$chains,
            length = job$length,
            success = run_status,
            stringsAsFactors = FALSE
        ))
    },
    mc.cores = min(4L, parallel::detectCores() - 1L),
    mc.preschedule = FALSE # Répartition dynamique pour lisser la charge des jobs longs
)

t1 <- Sys.time()
message("Inférence de masse terminée. Temps total : ", round(difftime(t1, t0, units = "mins"), 2), " minutes.")


# [ NODE_ACCESS : One sample ]
# ----------------------------------------------------- _


# Isolation des I/O : Un sous-dossier par échantillon
job <- list(dataset = "AML-59-001", chains = 4, length = 200)
sample_dir <- file.path(base_output_dir, job$dataset)
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

# Identification stricte de la sortie
run_name <- sprintf("%s_c%d_l%d", job$dataset, job$chains, job$length)
prefix_out <- file.path(sample_dir, run_name)

# Extraction des pointeurs mémoires depuis la liste consolidée
compass_inputs <- data_id_processedFiles_compass[[job$dataset]]


# Barrière tryCatch pour isoler les Segfaults potentiels d'un dataset corrompu
run_status <- tryCatch({

    vec_locus_names <- rownames(compass_inputs$variant_matrices$REF)

    vec_locus_chrom <- as.character(compass_inputs$locus_regions)


    # 2. Extraction des métadonnées CNA
    # cna_rowData <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
    # Agrégation Spatiale (Garantir le format CHR_GENE pour le C++)
    # vec_region_names <- paste0(cna_rowData$chrom, "_", cna_rowData$gene) # Ex: "17_TP53"
    vec_region_names <- rownames(compass_inputs$region_matrix)
    vec_region_names <- unique(vec_region_names) # Alignés avec les colonnes de ta matrice C_mat

    print("vec_locus_names")
    print(vec_locus_names)
    print("vec_locus_chrom")
    print(vec_locus_chrom)
    print("vec_region_names")
    print(vec_region_names)
    stop()

    run_compass_mcmc(
        variant_matrices  = compass_inputs$variant_matrices,
        locus_regions     = compass_inputs$locus_regions,
        region_matrix     = compass_inputs$region_matrix,
        output_prefix     = prefix_out,
        locus_names       = vec_locus_names,
        locus_chromosomes = vec_locus_chrom,
        region_names      = vec_region_names,
        chains            = job$chains,
        chain_length      = job$length,
        patient_sex       = "female"
    )
}, error = function(e) {
    warning(sprintf("\n[Worker %02d] Échec critique sur %s: %s", i, run_name, e$message))
    return(FALSE)
})

