# File: R/run_compass.R

#' @useDynLib ScIGMA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Execute COMPASS MCMC Phylogeny (In-Memory)
#'
#' @description
#' Direct interface to the ETHZ C++ COMPASS backend. Runs natively in R memory.
#' This function strictly requires Copy Number Alteration (CNA) data.
#'
#' @param variant_matrices List. Contains 'REF', 'ALT', and 'GT' matrices.
#' @param locus_regions Integer vector. Maps each locus to a region index (0-based).
#' @param region_matrix Integer matrix. Region counts (regions x cells). Required.
#' @param output_prefix Character. Prefix for output tree and data.
#' @param chains Integer. Number of MCMC chains.
#' @param chain_length Integer. Iterations per chain.
#' @param patient_sex Character. "female" or "male".
#'
#' @return Logical TRUE on success.
#' @export
run_compass_mcmc <- function(
        variant_matrices,
        locus_regions,
        region_matrix,
        output_prefix,
        chains = 4L,
        chain_length = 5000L,
        patient_sex = "female"
) {

    # Barrière de sécurité stricte
    if ( length(region_matrix) == 0 || nrow(region_matrix) == 0 ) {
        stop("ScIGMA requires a valid CNA region matrix to run COMPASS.")
    }

    # Neutralisation des NAs
    variant_matrices$REF[is.na(variant_matrices$REF)] <- 0L
    variant_matrices$ALT[is.na(variant_matrices$ALT)] <- 0L
    variant_matrices$GT[is.na(variant_matrices$GT)] <- 0L

    # NEW : Extraction des noms biologiques canoniques
    locus_names <- rownames(variant_matrices$REF)
    if ( is.null(locus_names) ) {
        locus_names <- paste0("Locus_", seq_len(nrow(variant_matrices$REF)) - 1L)
    }

    message("Initializing COMPASS C++ backend (In-Memory)...")

    execution_status <- tryCatch({
        run_compass_inference_cpp(
            ref_counts = variant_matrices$REF,
            alt_counts = variant_matrices$ALT,
            genotypes = variant_matrices$GT,
            locus_region_mapping = as.integer(locus_regions),
            region_counts = region_matrix,
            locus_names = as.character(locus_names), # NEW : Injection vers Rcpp
            output_prefix = output_prefix,
            n_chains = as.integer(chains),
            chain_length = as.integer(chain_length),
            use_cna = TRUE,
            sex = patient_sex
        )
        TRUE
    }, error = function(e) {
        stop("COMPASS execution failed: ", e$message)
    })

    return(execution_status)
}
