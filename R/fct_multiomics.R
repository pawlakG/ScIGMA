#' Extract variant genotypes from MultiAssayExperiment or SummarizedExperiment
#' @param mae_data MultiAssayExperiment or SummarizedExperiment object
#' @param variant_id Character string of the targeted variant
#' @param use_compass Logical to select imputed vs raw assay
#' @return data.frame with Barcode and Variant_Genotype columns
extract_variant_genotypes <- function(mae_data, variant_id, use_compass) {
    if (use_compass) {
        imputed_mtx <- SummarizedExperiment::assay(
            mae_data[["dna_variants"]],
            "compass_imputed"
        )

        if (
            !is.null(imputed_mtx) && all(variant_id %in% rownames(imputed_mtx))
        ) {
            variant_vector <- imputed_mtx[variant_id, ]
        } else {
            variant_vector <- SummarizedExperiment::assay(
                mae_data[["dna_variants"]],
                "gt"
            )[variant_id, ]
        }
    } else {
        # Standard extraction (Raw) from the DNA SummarizedExperiment
        variant_vector <- SummarizedExperiment::assay(
            mae_data[["dna_variants"]],
            "gt"
        )[variant_id, ]
    }

    extracted_df <- variant_vector |>
        as.data.frame() |>
        rownames_to_column("cell_barcode")
    colnames(extracted_df) <- c("cell_barcode", "variant_vector")
    message("inside extract_variant before mutate variant_vector")
    message("\n")

    extracted_df <- extracted_df %>%
        dplyr::mutate(
            variant_vector = dplyr::recode(
                variant_vector,
                "0" = "WT",
                "1" = "HET",
                "2" = "HOM",
                "3" = "Missing/ADO"
            )
        )
    colnames(extracted_df)[
        colnames(extracted_df) == "variant_vector"
    ] <- "Variant_Genotype"
    colnames(extracted_df)[
        colnames(extracted_df) == "cell_barcode"
    ] <- "Barcode"
    message("inside extract_variant after mutate variant_vector")
    message("\n")

    return(extracted_df)
}

#' Compute genotype distribution for a specific cell population
#' @param mae_data MultiAssayExperiment object
#' @param variant_ids Character vector of targeted variants
#' @param cell_barcodes Character vector of barcodes in the population
#' @param use_compass Logical to use imputed data
#' @noRd
compute_population_genotype_distribution <- function(
    mae_data,
    variant_ids,
    cell_barcodes,
    use_compass,
    seurat_cluster = NULL
) {
    if (use_compass) {
        mtx_source <- SummarizedExperiment::assay(
            mae_data[["dna_variants"]],
            "compass_imputed"
        )
    } else {
        mtx_source <- SummarizedExperiment::assay(
            mae_data[["dna_variants"]],
            "gt"
        )
    }

    # 2. Filtering and Extraction (Keep only gate cells)
    common_cells <- intersect(cell_barcodes, colnames(mtx_source))
    if (length(common_cells) == 0) {
        return(data.frame())
    }

    sub_mtx <- mtx_source[variant_ids, common_cells, drop = FALSE]

    dist_df <- as.data.frame(as.matrix(sub_mtx)) |>
        tibble::rownames_to_column("Variant_ID") |>
        tidyr::pivot_longer(
            -Variant_ID,
            names_to = "Barcode",
            values_to = "Code"
        ) |>
        dplyr::mutate(
            Variant_Genotype = dplyr::recode(
                as.character(Code),
                "0" = "WT",
                "1" = "HET",
                "2" = "HOM",
                "3" = "Missing/ADO"
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
    final_stats$Variant_Genotype <- factor(
        final_stats$Variant_Genotype,
        levels = levels_genotypes
    )

    return(final_stats)
}

#' Generate Multiomics UMAP with DNA Clones
#' @export
generate_multi_umap_dna_clones <- function(obj, use_compass = TRUE) {
    if (
        is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)
    ) {
        stop("UMAP not computed.")
    }
    if (is.null(obj$dna.clones)) {
        stop("DNA clones not assigned.")
    }

    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$barcode <- rownames(umap_df)

    umap_df$dna_clones <- as.character(obj$dna.clones[umap_df$barcode])
    umap_df$dna_clones[is.na(umap_df$dna_clones)] <- "Missing"

    p <- ggplot2::ggplot(
        umap_df,
        ggplot2::aes(x = umap_1, y = umap_2, color = dna_clones)
    ) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")

    if (!is.null(obj$dna_clone_colors)) {
        p <- p +
            ggplot2::scale_color_manual(
                values = c(obj$dna_clone_colors, "Missing" = "#e0e0e0")
            )
    }
    return(p)
}

#' Generate Multiomics UMAP with DNA Genotype
#' @export
generate_multi_umap_dna_genotype <- function(
    obj,
    variant_id,
    use_compass = TRUE
) {
    if (
        is.null(obj$seurat_object) || is.null(obj$seurat_object@reductions$umap)
    ) {
        stop("UMAP not computed.")
    }

    umap_df <- as.data.frame(obj$seurat_object@reductions$umap@cell.embeddings)
    umap_df$Barcode <- rownames(umap_df)

    geno_df <- extract_variant_genotypes(obj$mae, variant_id, use_compass)
    plot_df <- merge(umap_df, geno_df, by = "Barcode", all.x = TRUE)
    plot_df$Variant_Genotype[is.na(plot_df$Variant_Genotype)] <- "Missing/ADO"
    plot_df$Variant_Genotype <- factor(
        plot_df$Variant_Genotype,
        levels = c("WT", "HET", "HOM", "Missing/ADO")
    )

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = umap_1, y = umap_2, color = Variant_Genotype)
    ) +
        ggplot2::geom_point(size = 1.5, alpha = 0.8) +
        ggplot2::xlab("UMAP 1") +
        ggplot2::ylab("UMAP 2") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::scale_color_manual(
            values = c(
                "WT" = "#d3d3d3",
                "HET" = "#f39c12",
                "HOM" = "#c0392b",
                "Missing/ADO" = "#e0e0e0"
            )
        )
    return(p)
}

#' Generate Multiomics Barplot of DNA Clones vs Unsupervised Clusters
#' @export
generate_multi_barplot_clones_per_cluster <- function(obj) {
    if (
        is.null(obj$seurat_object) ||
            !"seurat_clusters" %in% colnames(obj$seurat_object@meta.data)
    ) {
        stop("Clusters not computed.")
    }
    if (is.null(obj$dna.clones)) {
        stop("DNA clones not assigned.")
    }

    clusters <- obj$seurat_object@meta.data$seurat_clusters
    names(clusters) <- rownames(obj$seurat_object@meta.data)
    clones <- obj$dna.clones[names(clusters)]
    clones[is.na(clones)] <- "Missing"

    plot_df <- data.frame(Cluster = clusters, Clone = clones) |>
        dplyr::group_by(Cluster, Clone) |>
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
        dplyr::group_by(Cluster) |>
        dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
        dplyr::ungroup()

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = Cluster, y = Percentage, fill = Clone)
    ) +
        # ggplot2::geom_bar(position = "fill") +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::ylab("Percentage (%)") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")

    if (!is.null(obj$dna_clone_colors)) {
        p <- p +
            ggplot2::scale_fill_manual(
                values = c(obj$dna_clone_colors, "Missing" = "#e0e0e0")
            )
    }
    return(p)
}

#' Generate Multiomics Barplot of DNA Variants vs Protein Clusters
#' @export
generate_multi_barplot_variant_genotype <- function(
    obj,
    variant_ids,
    use_compass = TRUE
) {
    if (
        is.null(obj$seurat_object) ||
            !"seurat_clusters" %in% colnames(obj$seurat_object@meta.data)
    ) {
        stop("Clusters not computed.")
    }

    clusters <- obj$seurat_object@meta.data$seurat_clusters
    names(clusters) <- rownames(obj$seurat_object@meta.data)

    stats_df <- compute_population_genotype_distribution(
        mae_data = obj$mae,
        variant_ids = variant_ids,
        cell_barcodes = names(clusters),
        use_compass = use_compass,
        seurat_cluster = clusters
    ) |>
        as.data.frame()

    variant_colors <- c(
        "WT" = "#440154",
        "HET" = "#21918c",
        "HOM" = "#fde725",
        "Missing/ADO" = "#e0e0e0"
    )

    variants_filtered_tmp <- obj$variants.filtered |>
        dplyr::select(-dplyr::any_of("variant_id")) |>
        tibble::rownames_to_column("Variant_ID")

    plot_df_joined <- dplyr::left_join(
        stats_df,
        variants_filtered_tmp,
        by = "Variant_ID"
    )

    id_mapping <- plot_df_joined |>
        dplyr::select(Variant_ID, gene, cdna) |>
        dplyr::distinct() |>
        dplyr::mutate(
            cdna_clean = ifelse(is.na(cdna) | cdna == "", "Unknown", cdna),
            Variant_Base = paste(gene, cdna_clean, sep = " - ")
        ) |>
        dplyr::group_by(Variant_Base) |>
        dplyr::mutate(n_variants = dplyr::n()) |>
        dplyr::ungroup() |>
        dplyr::mutate(
            Variant = ifelse(
                n_variants > 1,
                paste0(Variant_Base, " [", Variant_ID, "]"),
                Variant_Base
            )
        )

    plot_df_joined <- plot_df_joined |>
        dplyr::left_join(
            id_mapping |> dplyr::select(Variant_ID, Variant),
            by = "Variant_ID"
        )

    plot_df_joined <- plot_df_joined |>
        dplyr::arrange(as.numeric(Cluster), Variant) |>
        dplyr::mutate(
            X_Axis_Label = paste0(
                "<b>C",
                Cluster,
                "</b> | ",
                Variant
            ),
            X_Axis_Label = factor(
                X_Axis_Label,
                levels = unique(X_Axis_Label)
            )
        ) |>
        tidyr::complete(
            tidyr::nesting(X_Axis_Label, Cluster, Variant),
            Variant_Genotype,
            fill = list(Percentage = 0, Count = 0)
        ) |>
        dplyr::arrange(X_Axis_Label, Variant_Genotype)

    prism_axis_style <- list(
        showline = TRUE,
        linewidth = 2,
        linecolor = "black",
        ticks = "outside",
        tickwidth = 2,
        tickcolor = "black",
        ticklen = 5,
        mirror = FALSE,
        titlefont = list(size = 14, color = "black", family = "Arial"),
        tickfont = list(size = 12, color = "black", family = "Arial")
    )

    p <- plotly::plot_ly(
        data = plot_df_joined,
        x = ~X_Axis_Label,
        y = ~Percentage,
        color = ~Variant_Genotype,
        colors = variant_colors,
        type = "bar",
        height = 700,
        text = ~ ifelse(
            Percentage >= 5,
            paste0(round(Percentage, 1), "%"),
            ""
        ),
        textposition = "inside",
        insidetextanchor = "middle",
        textfont = list(size = 11, color = "black", family = "Arial"),
        constraintext = "none",
        hovertext = ~ paste0(
            "<b>Cluster:</b> ",
            Cluster,
            "<br>",
            "<b>Variant:</b> ",
            Variant,
            "<br>",
            "<b>Genotype:</b> ",
            Variant_Genotype,
            "<br>",
            "<b>Frequency:</b> ",
            round(Percentage, 1),
            "% (n=",
            Count,
            ")"
        ),
        hoverinfo = "text"
    ) |>
        plotly::layout(
            barmode = "stack",
            xaxis = c(
                list(title = "", tickangle = -45, automargin = TRUE),
                prism_axis_style
            ),
            yaxis = c(
                list(title = "<b>Frequency (%)</b>", range = c(0, 105)),
                prism_axis_style
            ),
            legend = list(title = list(text = "<b>Genotype</b>")),
            margin = list(b = 140, t = 50)
        ) |>
        plotly::config(displaylogo = FALSE)

    return(p)
}

#' Generate Multiomics Barplot of Gating Subpopulation Genotype vs Variants
#' @export
generate_multi_barplot_gating_genotype <- function(
    obj,
    variant_ids,
    use_compass = TRUE,
    target_gate = "All"
) {
    if (is.null(obj$protein_gates)) {
        stop("No protein gates defined.")
    }
    if (!target_gate == "All" & !target_gate %in% names(obj$protein_gates)) {
        stop(paste("Gate", target_gate, "not found in protein gates."))
    }

    if (target_gate == "All") {
        gate_cells <- unique(unlist(obj$protein_gates))
    } else {
        gate_cells <- obj$protein_gates[[target_gate]]
    }
    gate_assignment <- rep(target_gate, length(gate_cells))
    names(gate_assignment) <- gate_cells

    stats_df <- compute_population_genotype_distribution(
        mae_data = obj$mae,
        variant_ids = variant_ids,
        cell_barcodes = names(gate_assignment),
        use_compass = use_compass,
        seurat_cluster = gate_assignment
    ) |>
        as.data.frame()

    # Join with variant information to get clean names
    variants_filtered_tmp <- obj$variants.filtered |>
        dplyr::select(-dplyr::any_of("variant_id")) |>
        tibble::rownames_to_column("Variant_ID")

    plot_df <- dplyr::left_join(
        stats_df,
        variants_filtered_tmp,
        by = "Variant_ID"
    )

    plot_df <- plot_df |>
        dplyr::mutate(
            cdna_clean = ifelse(is.na(cdna) | cdna == "", "Unknown", cdna),
            Variant = paste(gene, cdna_clean, sep = " - ")
        ) |>
        dplyr::group_by(Variant) |>
        dplyr::mutate(n_variants = dplyr::n_distinct(Variant_ID)) |>
        dplyr::ungroup() |>
        dplyr::mutate(
            Variant = ifelse(
                n_variants > 1,
                paste0(Variant, " [", Variant_ID, "]"),
                Variant
            )
        ) |>
        tidyr::complete(
            Variant,
            Variant_Genotype,
            fill = list(Percentage = 0, Count = 0)
        )

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = Variant, y = Percentage, fill = Variant_Genotype)
    ) +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::ylab("Percentage (%)") +
        ggplot2::xlab(NULL) +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold") +
        ggplot2::scale_fill_manual(
            values = c(
                "WT" = "#d3d3d3",
                "HET" = "#f39c12",
                "HOM" = "#c0392b",
                "Missing/ADO" = "#e0e0e0"
            ),
            drop = FALSE
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

    return(p)
}

#' Generate Multiomics Barplot of DNA Clones vs Gating Subpopulations
#' @export
generate_multi_barplot_clones_per_gate <- function(obj) {
    if (is.null(obj$protein_gates)) {
        stop("No protein gates defined.")
    }
    if (is.null(obj$dna.clones)) {
        stop("DNA clones not assigned.")
    }

    all_cells <- unique(unlist(obj$protein_gates))
    # We assign cells to their last defined gate if overlapping
    gate_assignment <- character(length(all_cells))
    names(gate_assignment) <- all_cells

    for (gate_name in names(obj$protein_gates)) {
        gate_cells <- obj$protein_gates[[gate_name]]
        gate_assignment[gate_cells] <- gate_name
    }

    clones <- obj$dna.clones[names(gate_assignment)]
    clones[is.na(clones)] <- "Missing"

    plot_df <- data.frame(Gate = gate_assignment, Clone = clones) |>
        dplyr::group_by(Gate, Clone) |>
        dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
        dplyr::group_by(Gate) |>
        dplyr::mutate(Percentage = Count / sum(Count) * 100) |>
        dplyr::ungroup()

    p <- ggplot2::ggplot(
        plot_df,
        ggplot2::aes(x = Gate, y = Percentage, fill = Clone)
    ) +
        # ggplot2::geom_bar(position = "fill") +
        ggplot2::geom_bar(stat = "identity", position = "dodge") +
        ggplot2::ylab("Percentage (%)") +
        ggprism::theme_prism(base_size = 14, base_fontface = "bold")

    if (!is.null(obj$dna_clone_colors)) {
        p <- p +
            ggplot2::scale_fill_manual(
                values = c(obj$dna_clone_colors, "Missing" = "#e0e0e0")
            )
    }
    return(p)
}
