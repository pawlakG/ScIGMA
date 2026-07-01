#' Extract variant genotypes from MultiAssayExperiment or SummarizedExperiment
#' @param mae_data MultiAssayExperiment or SummarizedExperiment object
#' @param variant_id Character string of the targeted variant
#' @param use_compass Logical to select imputed vs raw assay
#' @return data.frame with Barcode and Variant_Genotype columns
extract_variant_genotypes <- function(mae_data, variant_id, use_compass) {
    if (use_compass) {
        # imputed_mtx_tmp <- S4Vectors::metadata(mae_data)
        # imputed_mtx <- imputed_mtx_tmp$compass$imputed_gt
        imputed_mtx <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "compass_imputed")

        if (!is.null(imputed_mtx) && all(variant_id %in% rownames(imputed_mtx))) {
            variant_vector <- imputed_mtx[variant_id, ]
        } else {
            variant_vector <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "gt")[variant_id, ]
        }
    } else {
        # Extraction standard (Raw) depuis le SummarizedExperiment DNA
        variant_vector <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "gt")[variant_id, ]
    }

    extracted_df <- variant_vector |>
        as.data.frame() |>
        rownames_to_column("cell_barcode")
    message("inside extract_variant before mutate variant_vector")
    # print(head(extracted_df))
    # print(str(extracted_df))
    # print(table(extracted_df$variant_vector))
    message("\n")

    extracted_df <- extracted_df %>% dplyr::mutate(variant_vector = dplyr::recode(variant_vector,
        "0" = "WT",
        "1" = "HET",
        "2" = "HOM",
        "3" = "Missing/ADO"
    ))
    colnames(extracted_df)[colnames(extracted_df) == "variant_vector"] <- "Variant_Genotype"
    colnames(extracted_df)[colnames(extracted_df) == "cell_barcode"] <- "Barcode"
    message("inside extract_variant after mutate variant_vector")
    # print(head(extracted_df))
    # print(str(extracted_df))
    # print(table(extracted_df$Variant_Genotype))
    message("\n")

    return(extracted_df)
}

#' Compute genotype distribution for a specific cell population
#' @param mae_data MultiAssayExperiment object
#' @param variant_ids Character vector of targeted variants
#' @param cell_barcodes Character vector of barcodes in the population
#' @param use_compass Logical to use imputed data
#' @noRd
compute_population_genotype_distribution <- function(mae_data,
                                                    variant_ids,
                                                    cell_barcodes,
                                                    use_compass,
                                                    seurat_cluster = NULL) {
    if (use_compass) {
        mtx_source <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "compass_imputed")
    } else {
        mtx_source <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "gt")
    }

    # 2. Filtrage et Extraction (Garder uniquement les cellules de la gate)
    common_cells <- intersect(cell_barcodes, colnames(mtx_source))
    if (length(common_cells) == 0) {
        return(data.frame())
    }

    sub_mtx <- mtx_source[variant_ids, common_cells, drop = FALSE]

    dist_df <- as.data.frame(as.matrix(sub_mtx)) |>
        tibble::rownames_to_column("Variant_ID") |>
        tidyr::pivot_longer(-Variant_ID, names_to = "Barcode", values_to = "Code") |>
        dplyr::mutate(
            Variant_Genotype = dplyr::recode(as.character(Code),
                "0" = "WT", "1" = "HET",
                "2" = "HOM", "3" = "Missing/ADO"
            )
        )

    if (is.null(seurat_cluster)) {
        final_stats <- dist_df |>
            dplyr::group_by(Variant_ID, Variant_Genotype) |>
            dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
            dplyr::group_by(Variant_ID)
    } else {
        dist_df <- dist_df |>
            dplyr::mutate(
                Cluster = seurat_cluster[Barcode]
            )
        final_stats <- dist_df |>
            dplyr::group_by(Variant_ID, Cluster, Variant_Genotype) |>
            dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
            dplyr::group_by(Variant_ID, Cluster)
    }

    final_stats <- final_stats |>
        dplyr::mutate(
            Total_In_Variant = sum(Count),
            Percentage = (Count / Total_In_Variant) * 100
        ) |>
        dplyr::ungroup()


    levels_genotypes <- c("WT", "HET", "HOM", "Missing/ADO")
    final_stats$Variant_Genotype <- factor(final_stats$Variant_Genotype, levels = levels_genotypes)

    return(final_stats)
}

#' Generate Multiomics UMAP with DNA Clones
#' @export
generate_multi_umap_dna_clones <- function(obj, use_compass = TRUE) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) stop("UMAP not computed.")
    if (is.null(obj$dna.clones)) stop("DNA clones not assigned.")
    
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$barcode <- rownames(umap_df)
    
    umap_df$dna_clones <- as.character(obj$dna.clones[umap_df$barcode])
    umap_df$dna_clones[is.na(umap_df$dna_clones)] <- "Missing"
    
    p <- ggplot2::ggplot(umap_df, ggplot2::aes(x = umap_1, y = umap_2, color = dna_clones)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")
    
    if (!is.null(obj$dna_clone_colors)) {
        p <- p + ggplot2::scale_color_manual(values = c(obj$dna_clone_colors, "Missing" = "#e0e0e0"))
    }
    return(p)
}

#' Generate Multiomics UMAP with DNA Genotype
#' @export
generate_multi_umap_dna_genotype <- function(obj, variant_id, use_compass = TRUE) {
    if (is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)) stop("UMAP not computed.")
    
    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$Barcode <- rownames(umap_df)
    
    geno_df <- extract_variant_genotypes(obj$mae, variant_id, use_compass)
    plot_df <- merge(umap_df, geno_df, by = "Barcode", all.x = TRUE)
    plot_df$Variant_Genotype[is.na(plot_df$Variant_Genotype)] <- "Missing/ADO"
    plot_df$Variant_Genotype <- factor(plot_df$Variant_Genotype, levels = c("WT", "HET", "HOM", "Missing/ADO"))
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = umap_1, y = umap_2, color = Variant_Genotype)) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::scale_color_manual(values = c("WT" = "#d3d3d3", "HET" = "#f39c12", "HOM" = "#c0392b", "Missing/ADO" = "#e0e0e0"))
    return(p)
}

#' Generate Multiomics Barplot of DNA Clones vs Unsupervised Clusters
#' @export
generate_multi_barplot_clones_per_cluster <- function(obj) {
    if (is.null(obj$seurat_object) || !"seurat_clusters" %in% colnames(obj$seurat_object@meta.data)) stop("Clusters not computed.")
    if (is.null(obj$dna.clones)) stop("DNA clones not assigned.")
    
    clusters <- obj$seurat_object@meta.data$seurat_clusters
    clones <- obj$dna.clones[names(clusters)]
    clones[is.na(clones)] <- "Missing"
    
    plot_df <- data.frame(Cluster = clusters, Clone = clones)
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Cluster, fill = Clone)) +
        ggplot2::geom_bar(position = "fill") +
        ggplot2::ylab("Proportion") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")
        
    if (!is.null(obj$dna_clone_colors)) {
        p <- p + ggplot2::scale_fill_manual(values = c(obj$dna_clone_colors, "Missing" = "#e0e0e0"))
    }
    return(p)
}

#' Generate Multiomics Barplot of DNA Variants vs Protein Clusters
#' @export
generate_multi_barplot_variant_genotype <- function(obj, variant_ids, use_compass = TRUE) {
    if (is.null(obj$seurat_object) || !"seurat_clusters" %in% colnames(obj$seurat_object@meta.data)) stop("Clusters not computed.")
    
    clusters <- obj$seurat_object@meta.data$seurat_clusters
    stats_df <- compute_population_genotype_distribution(
        mae_data = obj$mae, 
        variant_ids = variant_ids, 
        cell_barcodes = names(clusters), 
        use_compass = use_compass, 
        seurat_cluster = clusters
    )
    
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = Cluster, y = Percentage, fill = Variant_Genotype)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::facet_wrap(~Variant_ID) +
        ggplot2::ylab("Percentage (%)") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::scale_fill_manual(values = c("WT" = "#d3d3d3", "HET" = "#f39c12", "HOM" = "#c0392b", "Missing/ADO" = "#e0e0e0"))
    return(p)
}

#' Generate Multiomics Barplot of Gating Subpopulation Genotype vs Variants
#' @export
generate_multi_barplot_gating_genotype <- function(obj, variant_ids, use_compass = TRUE) {
    if (is.null(obj$protein_gates)) stop("No protein gates defined.")
    
    all_cells <- unique(unlist(obj$protein_gates))
    # We assign cells to their last defined gate if overlapping
    gate_assignment <- character(length(all_cells))
    names(gate_assignment) <- all_cells
    
    for (gate_name in names(obj$protein_gates)) {
        gate_cells <- obj$protein_gates[[gate_name]]
        gate_assignment[gate_cells] <- gate_name
    }
    
    stats_df <- compute_population_genotype_distribution(
        mae_data = obj$mae, 
        variant_ids = variant_ids, 
        cell_barcodes = names(gate_assignment), 
        use_compass = use_compass, 
        seurat_cluster = gate_assignment
    )
    
    # Rename Cluster to Gate for plot
    stats_df$Gate <- stats_df$Cluster
    
    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = Gate, y = Percentage, fill = Variant_Genotype)) +
        ggplot2::geom_bar(stat = "identity", position = "stack") +
        ggplot2::facet_wrap(~Variant_ID) +
        ggplot2::ylab("Percentage (%)") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::scale_fill_manual(values = c("WT" = "#d3d3d3", "HET" = "#f39c12", "HOM" = "#c0392b", "Missing/ADO" = "#e0e0e0"))
    return(p)
}
