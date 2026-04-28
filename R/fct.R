
sanitize_filename <- function(filename, replacement = "_") {

    # 1. (Optional) Handle accents/special characters
    # Converts "Hélène" to "Helene".
    # Note: `iconv` results can vary slightly depending on OS locale.
    filename <- iconv(filename, to = "ASCII//TRANSLIT")

    # 2. Replace spaces with the replacement character
    filename <- gsub(" ", replacement, filename)

    # 3. Replace any character that is NOT alphanumeric, a dot, or a dash
    # Regex explanation:
    # [^...]      : Match any character NOT in this set
    # [:alnum:]   : Alphanumeric characters (letters and numbers)
    # \\.         : A literal dot (needs double escape)
    # -           : A literal dash
    filename <- gsub("[^[:alnum:]\\.-]", replacement, filename)

    # 4. Remove repeated instances of the replacement character
    # Example: "My___File" becomes "My_File"
    pattern_repeat <- paste0(replacement, "+")
    filename <- gsub(pattern_repeat, replacement, filename)

    # 5. Remove the replacement character from the start or end of the string
    # Example: "_file_" becomes "file"
    pattern_edge <- paste0("^", replacement, "|", replacement, "$")
    filename <- gsub(pattern_edge, "", filename)

    return(filename)
}


# Dans la définition de R6Class "ScIGMA_object" -> public list:

#' @description
#' Update cluster labels safely.
#' @param new_labels A named vector where names are old levels and values are new labels.
update_cluster_labels = function(new_labels) {
    # Vérification défensive
    if (is.null(self$meta_data) || is.null(self$dna.clones)) {
        stop("No clusters found in meta_data")
    }

    # On s'assure que c'est un facteur
    current_clusters <- as.factor(self$dna.clones)

    # Mise à jour rapide des niveaux (levels)
    # match() est extrêmement rapide pour faire la correspondance
    # On ne change que les niveaux présents dans new_labels
    current_levels <- levels(current_clusters)

    # On identifie les indices des niveaux à changer
    indices <- match(names(new_labels), current_levels)

    # On applique les changements uniquement sur les niveaux valides
    valid_updates <- !is.na(indices)
    if (any(valid_updates)) {
        levels(current_clusters)[indices[valid_updates]] <- new_labels[valid_updates]

        # Sauvegarde dans l'objet
        self$dna.clones <- current_clusters

        # Message log
        message(sprintf("Updated %d cluster labels.", sum(valid_updates)))
    }
}


protein_run_pca <- function(ScIGMA_data) {
    message("Running Protein PCA Pipeline...")

    # 1. Extraction (Protéines x Cellules)
    prot_counts <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], "counts")

    # 2. Normalisation In-Memory (Sécurisé pour les protéines)
    norm_mat <- normalize_linear_regression(as.matrix(prot_counts), jitter = 0.5)

    # 2.1 Mise à jour des métadonnées
    S4Vectors::metadata(ScIGMA_data$mae)$protein_normalize_method <- "Normalized"

    # 2.2 Mise à jour de l'indicateut de filtrage des proteines
    ScIGMA_data$protein.filtered <- TRUE



    # 3. Injection canonique dans le MAE
    SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], "normalized") <- norm_mat

    # 4. Construction de l'objet Seurat (AUCUNE TRANSPOSITION REQUISE)
    min_features_threshold <- floor(sqrt(nrow(norm_mat))) # nrow = nb protéines

    seurat_obj <- Seurat::CreateSeuratObject(
        counts = norm_mat,
        project = "ScIGMA_proteins",
        min.cells = 3,
        min.features = min_features_threshold
    )

    # Assigne les données normalisées au slot "data" (Standard Seurat v4/v5)
    seurat_obj <- Seurat::SetAssayData(seurat_obj, layer = "data", new.data = norm_mat)

    # 5. Pipeline Dimensionnel
    seurat_obj <- Seurat::FindVariableFeatures(
        seurat_obj,
        selection.method = "vst",
        nfeatures = nrow(seurat_obj)
    )

    seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj))

    # Sécurité algorithmique : Les PC max ne peuvent pas excéder le nombre de protéines
    max_npcs <- min(nrow(seurat_obj) - 1, ncol(seurat_obj) - 1, 50)

    seurat_obj <- Seurat::RunPCA(
        seurat_obj,
        features = Seurat::VariableFeatures(object = seurat_obj),
        npcs = max_npcs,
        verbose = FALSE
    )

    seurat_obj <- Seurat::FindNeighbors(
        seurat_obj,
        dims = 1:max_npcs,
        verbose = FALSE
    )

    invisible(seurat_obj)
}


#' Sanitize manual gating names to ensure strictly valid R strings
#' @param gate_name Character vector of raw gate names
#' @return Character vector of sanitized names
sanitize_gate_name <- function(gate_name) {
    clean_name <- trimws(gate_name)
    clean_name <- gsub("[^A-Za-z0-9]", "_", clean_name)
    clean_name <- gsub("_+", "_", clean_name)
    clean_name <- gsub("^_|_$", "", clean_name)
    return(clean_name)
}
