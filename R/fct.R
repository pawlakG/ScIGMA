
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


protein_run_pca <- function(ScIGMA_data){

    ScIGMA_data$protein.mtx.filtered.normalized <- normalize_linear_regression(as.matrix(ScIGMA_data$protein.mtx), jitter = 0.5)

    ScIGMA_data$seurat_object <- CreateSeuratObject(counts = t(ScIGMA_data$protein.mtx.filtered.normalized) ,
                                                    project = "ScIGMA_data",
                                                    min.cells = 3,
                                                    min.features = floor(sqrt(ncol(ScIGMA_data$protein.mtx.filtered.normalized))))

    ScIGMA_data$seurat_object@assays$RNA$data <- t(ScIGMA_data$protein.mtx.filtered.normalized)

    ScIGMA_data$seurat_object <- FindVariableFeatures(ScIGMA_data$seurat_object,
                                                      selection.method = "vst",
                                                      nfeatures = nrow(ScIGMA_data$seurat_object))

    ScIGMA_data$seurat_object <- ScaleData(ScIGMA_data$seurat_object, features = rownames(ScIGMA_data$seurat_object))
    ScIGMA_data$seurat_object <- RunPCA(ScIGMA_data$seurat_object,
                                        features = VariableFeatures(object = ScIGMA_data$seurat_object),
                                        npcs = nrow(ScIGMA_data$seurat_object)-2)

    ScIGMA_data$seurat_object <- FindNeighbors(ScIGMA_data$seurat_object, dims = 1:(nrow(ScIGMA_data$seurat_object)-2))
    return(ScIGMA_data)
}


