# NEW
# File: scripts/prepare_compass_inputs.R

#' Parse COMPASS CSV outputs into memory-safe matrices for Rcpp
#'
#' @param variants_file Path to the [sample]_variants.csv file
#' @param regions_file Path to the [sample]_regions.csv file
#'
#' @return A list formatted strictly for run_compass_mcmc() API
#' @export
prepare_compass_from_csv <- function(variants_file, regions_file) {
    # Fast I/O and vectorized string ops
    if (!requireNamespace("data.table", quietly = TRUE) ||
        !requireNamespace("stringr", quietly = TRUE)) {
        stop("Packages 'data.table' and 'stringr' are required for fast parsing.")
    }

    message("Loading CSV files into RAM...")
    var_df <- data.table::fread(variants_file, data.table = FALSE)
    reg_df <- data.table::fread(regions_file, header = FALSE, data.table = FALSE)

    # ---- 1. Variants matrix extraction ----
    message("Parsing REF:ALT:GT multiplexed strings...")
    cell_cols <- 8:ncol(var_df)
    var_mat_raw <- as.matrix(var_df[, cell_cols])
    rownames(var_mat_raw) <- var_df$NAME

    # Fast vectorized split (C++ backed via stringr)
    parsed <- stringr::str_split_fixed(as.vector(var_mat_raw), ":", n = 3)

    mat_ref <- matrix(as.integer(parsed[, 1]),
                      nrow = nrow(var_mat_raw),
                      dimnames = dimnames(var_mat_raw))

    mat_alt <- matrix(as.integer(parsed[, 2]),
                      nrow = nrow(var_mat_raw),
                      dimnames = dimnames(var_mat_raw))

    mat_gt  <- matrix(as.integer(parsed[, 3]),
                      nrow = nrow(var_mat_raw),
                      dimnames = dimnames(var_mat_raw))

    # Neutralize Missing genotypes (3 = Missing in MissionBio topology)
    mat_gt[mat_gt == 3L | is.na(mat_gt)] <- NA_integer_

    # Rcpp Segfault prevention: Explicit coercion to integer
    storage.mode(mat_ref) <- "integer"
    storage.mode(mat_alt) <- "integer"
    storage.mode(mat_gt)  <- "integer"

    # ---- 2. Region matrix extraction ----
    message("Parsing CNA region read counts...")
    cna_mat <- as.matrix(reg_df[, -1, drop = FALSE])
    rownames(cna_mat) <- reg_df[[1]]
    colnames(cna_mat) <- colnames(var_mat_raw) # Strict cell alignment
    storage.mode(cna_mat) <- "integer"

    # ---- 3. Topological Mapping (Base-0) ----
    message("Computing Locus to Region mapping...")
    # MissionBio CSV prepends chroms (e.g. "5_NPM1"). We extract the true gene name.
    cna_pure_names <- sub("^[^_]+_", "", rownames(cna_mat))

    # Base-0 offset for C++ indexing
    locus_regions <- match(var_df$REGION, cna_pure_names) - 1L

    # Fail-Fast mechanism : Bloque l'exécution si les fichiers sont désynchronisés
    if (any(is.na(locus_regions))) {
        unmapped <- unique(var_df$REGION[is.na(locus_regions)])
        stop("Fatal: Variants mapped to unreferenced CNA regions: ",
             paste(unmapped, collapse = ", "))
    }

    message("COMPASS In-Memory structures successfully generated.")
    return(list(
        variant_matrices = list(REF = mat_ref, ALT = mat_alt, GT = mat_gt),
        locus_regions    = locus_regions,
        region_matrix    = cna_mat
    ))
}

#' Reconstruct imputed single-cell genotypes from COMPASS outputs
#'
#' @param prefix_out Character. Le chemin et préfixe utilisé dans run_compass_mcmc
#' @return Une matrice [Cellules x Locus] de type Integer (0 = REF, 1 = HET, 2 = HOM_ALT)
#' @export
get_imputed_genotypes <- function(prefix_out) {

    # 1. Localisation des artefacts générés par le C++
    nodes_gt_file <- paste0(prefix_out, "_nodes_genotypes.tsv")
    cell_assign_file <- paste0(prefix_out, "_cellAssignments.tsv")

    if (!file.exists(nodes_gt_file) || !file.exists(cell_assign_file)) {
        stop(sprintf("Fichiers introuvables pour le préfixe : %s. L'inférence a-t-elle convergé ?", prefix_out))
    }

    # 2. Chargement en RAM
    nodes_gt <- read.delim(nodes_gt_file, stringsAsFactors = FALSE)
    cell_assign <- read.delim(cell_assign_file, stringsAsFactors = FALSE)

    # 3. Filtrage strict des doublets
    # On exclut les cellules flaggées comme doublets par le modèle pour garantir la pureté clonale
    singlets <- cell_assign[cell_assign$doublet == "no", ]

    if (nrow(singlets) == 0) {
        stop("Aucun singulet trouvé. Matrice inexploitable.")
    }

    # 4. Préparation de la clé de jointure (Format: "Node X")
    rownames(nodes_gt) <- nodes_gt$node
    target_nodes <- paste0("Node ", singlets$node)

    # 5. Projection de la matrice des nœuds vers les cellules (Broadcasting)
    imputed_mat <- as.matrix(nodes_gt[target_nodes, -1, drop = FALSE])

    # 6. Restauration de l'identité et du typage
    rownames(imputed_mat) <- singlets$cell
    storage.mode(imputed_mat) <- "integer"

    message(sprintf("Matrice imputée reconstruite : %d cellules x %d variants.",
                    nrow(imputed_mat), ncol(imputed_mat)))

    return(imputed_mat)
}
