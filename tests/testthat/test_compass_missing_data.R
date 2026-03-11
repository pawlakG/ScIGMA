# File: tests/test_compass_missing_data.R
# Exécution stricte en console APRÈS devtools::load_all()

t0 <- Sys.time()

# 1. Génération de matrices avec données massives manquantes (Dropout)
n_loci <- 5L # Diminué pour une validation claire de la topologie
n_cells <- 6000L
n_regions <- 4L

set.seed(1337)

# Noms de variants "Publication Grade" (Oncologie classique)
locus_names <- c(
    "TP53_p.R175H",
    "KRAS_p.G12D",
    "PIK3CA_p.H1047R",
    "EGFR_p.L858R",
    "PTEN_p.R130G"
)

# Génération des données de base avec rownames
mat_ref <- matrix(sample(0L:50L, n_loci * n_cells, replace = TRUE), nrow = n_loci, dimnames = list(locus_names, NULL))
mat_alt <- matrix(sample(0L:10L, n_loci * n_cells, replace = TRUE), nrow = n_loci, dimnames = list(locus_names, NULL))
mat_gt  <- matrix(sample(0L:3L, n_loci * n_cells, replace = TRUE), nrow = n_loci, dimnames = list(locus_names, NULL))

# Injection brutale de 30% de NAs pour simuler l'absence d'information
idx_na <- sample(seq_along(mat_ref), size = floor(0.3 * length(mat_ref)))
mat_ref[idx_na] <- NA_integer_
mat_alt[idx_na] <- NA_integer_
mat_gt[idx_na]  <- NA_integer_

list_variants_na <- list(
    REF = mat_ref,
    ALT = mat_alt,
    GT  = mat_gt
)

# Matrices de régions
mat_regions <- matrix(
    sample(50L:200L, n_regions * n_cells, replace = TRUE),
    nrow = n_regions
)

# Mapping Base-0 (Les 3 premiers variants dans la région 0, les 2 derniers dans la région 1)
vec_mapping <- c(0L, 0L, 0L, 1L, 1L)

dir_output <- tempdir()
# 2. Exécution du test (Le MCMC doit combler les trous)
message("--- Démarrage du test MCMC avec matrices trouées (NAs) ---")
prefix_out <- file.path(dir_output, "test_mcmc_missing")

status_execution <- tryCatch({
    run_compass_mcmc(
        variant_matrices = list_variants_na,
        locus_regions = vec_mapping,
        region_matrix = mat_regions,
        output_prefix = prefix_out,
        chains = 2L,
        chain_length = 50L,
        patient_sex = "female"
    )
}, error = function(e) {
    message("ÉCHEC INTERCEPTÉ : ", e$message)
    FALSE
})

# 3. Validation de l'imputation en aval
if ( status_execution ) {
    message("SUCCÈS : Le backend C++ a digéré les zéros de couverture.")

    message("--- Démarrage de l'imputation probabiliste ---")

    # On extrait une nouvelle matrice sans altérer list_variants_na
    gt_imputed <- impute_compass_genotypes(
        gt_matrix = list_variants_na$GT,
        output_prefix = prefix_out,
        min_probability = 0.80
    )

    total_nas_initial <- sum(is.na(list_variants_na$GT))
    total_nas_final <- sum(is.na(gt_imputed))

    message(sprintf("NAs dans la matrice d'origine (intacte) : %d", total_nas_initial))
    message(sprintf("NAs dans la NOUVELLE matrice (gt_imputed) : %d", total_nas_final))

    if ( total_nas_final < total_nas_initial ) {
        message("VALIDATION : La nouvelle matrice complétée a été générée avec succès.")
    } else {
        message("ATTENTION : Aucune imputation réalisée (probabilités d'assignation trop faibles).")
    }
}

t1 <- Sys.time()

message(paste0("Time taken : ", format(difftime(t1, t0), units = "auto")))
