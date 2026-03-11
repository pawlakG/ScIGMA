# UPDATED
# File: R/impute_genotypes.R
# Modification de la signature et de l'isolation de la matrice

#' Impute Missing Genotypes from COMPASS Phylogeny
#'
#' @description
#' Extracts inferred genotypes from the C++ MCMC outputs and projects them
#' onto the missing data (NAs). Returns a new standalone matrix.
#'
#' @param gt_matrix Integer matrix. The original GT matrix containing NAs.
#' @param output_prefix Character. The prefix used during MCMC execution.
#' @param min_probability Numeric. Minimum MCMC assignment probability (0 to 1).
#'
#' @return Integer matrix. A new imputed GT matrix.
#' @export
impute_compass_genotypes <- function(
        gt_matrix,
        output_prefix,
        min_probability = 0.90
) {
    file_assign <- paste0(output_prefix, "_cellAssignments.tsv")
    file_probs  <- paste0(output_prefix, "_cellAssignmentProbs.tsv")
    file_nodes  <- paste0(output_prefix, "_nodes_genotypes.tsv")

    if ( !file.exists(file_assign) || !file.exists(file_nodes) ) {
        stop("Fatal: COMPASS MCMC output files not found. Run inference first.")
    }

    df_assign <- read.table(file_assign, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df_probs  <- read.table(file_probs, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df_nodes  <- read.table(file_nodes, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

    # Création d'une copie stricte et isolée
    mat_imputed <- gt_matrix

    cell_names_cpp <- paste0("Cell_", seq_len(ncol(mat_imputed)) - 1L)

    cells_imputed <- 0L
    na_filled <- 0L

    for ( j in seq_len(ncol(mat_imputed)) ) {
        target_cell <- cell_names_cpp[j]

        is_doublet <- df_assign$doublet[df_assign$cell == target_cell]
        if ( is_doublet == "yes" ) next

        cell_probs <- as.numeric(df_probs[df_probs$cell == target_cell, -1])
        if ( max(cell_probs) < min_probability ) next

        node_idx <- df_assign$node[df_assign$cell == target_cell]
        node_name <- paste0("Node ", node_idx)

        idx_na <- which(is.na(mat_imputed[, j]))

        if ( length(idx_na) > 0 ) {
            mat_imputed[idx_na, j] <- as.integer(df_nodes[node_name, idx_na])
            cells_imputed <- cells_imputed + 1L
            na_filled <- na_filled + length(idx_na)
        }
    }

    message(
        sprintf(
            "Imputation complete: %d NAs filled across %d high-confidence cells.",
            na_filled,
            cells_imputed
        )
    )

    return(mat_imputed)
}
