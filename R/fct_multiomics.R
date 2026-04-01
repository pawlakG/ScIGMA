#' Auto-annotate Seurat clusters based on surface markers
#' @param seurat_obj A Seurat object with protein assay
#' @param assay Name of the protein assay
auto_annotate_clusters <- function(seurat_obj, assay = "RNA") {
    # Dictionnaire standard (À adapter librement selon ton panel de protéines)
    marker_dict <- list(
        "T cells" = c("CD3", "CD4", "CD8", "CD3E"),
        "B cells" = c("CD19", "CD20", "MS4A1"),
        "Monocytes" = c("CD14", "CD11b", "ITGAM", "CD16", "FCGR3A"),
        "NK cells" = c("CD56", "NCAM1"),
        "HSC/Blasts" = c("CD34", "CD117", "KIT"),
        "Erythroid" = c("CD71", "TFRC", "CD235a", "GYPA"),
        "Platelets" = c("CD41", "ITGA2B", "CD61", "ITGB3")
    )

    # Extraction de l'expression moyenne par cluster
    # [[1]] permet de gérer dynamiquement le nom de l'assay renvoyé par Seurat
    cluster_avgs <- Seurat::AverageExpression(seurat_obj, assays = assay)[[1]]
    available_markers <- rownames(cluster_avgs)

    predicted_labels <- character(ncol(cluster_avgs))
    names(predicted_labels) <- colnames(cluster_avgs)

    for (clust in colnames(cluster_avgs)) {
        best_score <- -Inf
        best_label <- "Unknown"

        for (cell_type in names(marker_dict)) {
            intersect_markers <- intersect(marker_dict[[cell_type]], available_markers)

            if (length(intersect_markers) > 0) {
                # Score = Moyenne de l'expression des marqueurs détectés
                score <- mean(cluster_avgs[intersect_markers, clust])
                if (score > best_score && score > 0.5) { # Seuil minimum empirique
                    best_score <- score
                    best_label <- cell_type
                }
            }
        }
        predicted_labels[clust] <- best_label
    }

    # Concaténation (ex: "Cluster 0: T cells")
    final_labels <- paste0("Cluster ", names(predicted_labels), ": ", predicted_labels)
    names(final_labels) <- names(predicted_labels)

    return(final_labels)
}
