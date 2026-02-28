# UPDATED
# File: tests/test_compass_in_memory.R
# Exécution stricte en console APRÈS devtools::load_all()

# 1. Génération de matrices factices (Mock Data)
n_loci <- 5L
n_cells <- 10L
n_regions <- 2L

set.seed(42)

# Matrices de variants (Lignes = Loci, Colonnes = Cellules)
mat_ref <- matrix(sample(0L:50L, n_loci * n_cells, replace = TRUE), nrow = n_loci)
mat_alt <- matrix(sample(0L:10L, n_loci * n_cells, replace = TRUE), nrow = n_loci)
mat_gt <- matrix(sample(0L:3L, n_loci * n_cells, replace = TRUE), nrow = n_loci)

list_variants <- list(
    REF = mat_ref,
    ALT = mat_alt,
    GT = mat_gt
)

# Matrices de régions (Lignes = Régions, Colonnes = Cellules)
mat_regions <- matrix(
    sample(50L:200L, n_regions * n_cells, replace = TRUE),
    nrow = n_regions
)

# 2. Le point critique : Mapping Locus -> Région (BASE-0)
vec_mapping <- c(0L, 0L, 0L, 1L, 1L)

dir_output <- tempdir()

# 3. Exécution du test (Scénario Réaliste avec CNA)
message("--- Démarrage du test MCMC In-Memory ---")
prefix_out <- file.path(dir_output, "test_mcmc")

status_execution <- tryCatch({
    run_compass_mcmc(
        variant_matrices = list_variants,
        locus_regions = vec_mapping,
        region_matrix = mat_regions,
        output_prefix = prefix_out,
        chains = 2L,
        chain_length = 10L,
        patient_sex = "female"
    )
}, error = function(e) {
    message("ÉCHEC INTERCEPTÉ : ", e$message)
    FALSE
})

if ( status_execution ) {
    message("SUCCÈS : Le backend C++ a digéré les matrices et rendu la main.")
    print(list.files(dir_output, pattern = "test_mcmc"))
}

# 4. Exécution du test (Validation du verrouillage Fail-Fast sans CNA)
message("\n--- Démarrage du test de verrouillage API (Matrice Vide) ---")
prefix_out_no_cna <- file.path(dir_output, "test_mcmc_no_cna")

status_execution_no_cna <- tryCatch({
    run_compass_mcmc(
        variant_matrices = list_variants,
        locus_regions = vec_mapping,
        region_matrix = matrix(integer(0), nrow = 0, ncol = 0),
        output_prefix = prefix_out_no_cna,
        chains = 2L,
        chain_length = 10L,
        patient_sex = "female"
    )
}, error = function(e) {
    # L'erreur native R que nous avons injectée doit apparaître ici
    message("SUCCÈS DU VERROUILLAGE (Rejet propre par R) : ", e$message)
    FALSE
})
