#' Extract variant genotypes from MultiAssayExperiment or SummarizedExperiment
#' @param mae_data MultiAssayExperiment or SummarizedExperiment object
#' @param variant_id Character string of the targeted variant
#' @param use_compass Logical to select imputed vs raw assay
#' @return data.frame with Barcode and Variant_Genotype columns
extract_variant_genotypes <- function(mae_data, variant_id, use_compass) {
    if (use_compass) {
        # imputed_mtx_tmp <- S4Vectors::metadata(mae_data)
        # imputed_mtx <- imputed_mtx_tmp$compass$imputed_gt
        imputed_mtx <-  SummarizedExperiment::assay(mae_data[["dna_variants"]], "compass_imputed")

        if (!is.null(imputed_mtx) && all(variant_id %in% rownames(imputed_mtx))) {
            variant_vector <- imputed_mtx[variant_id, ]
        } else {
            variant_vector <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "gt")[variant_id, ]
        }
    } else {
        # Extraction standard (Raw) depuis le SummarizedExperiment DNA
        variant_vector <- SummarizedExperiment::assay(mae_data[["dna_variants"]], "gt")[variant_id, ]
    }

    extracted_df <- variant_vector |> as.data.frame() |> rownames_to_column("cell_barcode")
    message("inside extract_variant before mutate variant_vector")
    print(head(extracted_df))
    print(str(extracted_df))
    print(table(extracted_df$variant_vector))
    message("\n")

    extracted_df <- extracted_df %>% dplyr::mutate(variant_vector = dplyr::recode(variant_vector,
                                                                                  "0" = "WT",
                                                                                  "1" = "HET",
                                                                                  "2" = "HOM",
                                                                                  "3" = "Missing/ADO"))
    colnames(extracted_df)[colnames(extracted_df) == "variant_vector"] <- "Variant_Genotype"
    colnames(extracted_df)[colnames(extracted_df) == "cell_barcode"] <- "Barcode"
    message("inside extract_variant after mutate variant_vector")
    print(head(extracted_df))
    print(str(extracted_df))
    print(table(extracted_df$Variant_Genotype))
    message("\n")

    return(extracted_df)
}

#' Compute genotype distribution for a specific cell population
#' @param mae_data MultiAssayExperiment object
#' @param variant_ids Character vector of targeted variants
#' @param cell_barcodes Character vector of barcodes in the population
#' @param use_compass Logical to use imputed data
#' @export
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
    if (length(common_cells) == 0) return(data.frame())

    sub_mtx <- mtx_source[variant_ids, common_cells, drop = FALSE]

    dist_df <- as.data.frame(as.matrix(sub_mtx)) |>
        tibble::rownames_to_column("Variant_ID") |>
        tidyr::pivot_longer(-Variant_ID, names_to = "Barcode", values_to = "Code")|>
        dplyr::mutate(
            Variant_Genotype = dplyr::recode(as.character(Code),
                                             "0" = "WT", "1" = "HET",
                                             "2" = "HOM", "3" = "Missing/ADO")
        )

    if (is.null(seurat_cluster)){
        final_stats <- dist_df |>
            dplyr::group_by(Variant_ID, Variant_Genotype)|>
            dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
            dplyr::group_by(Variant_ID)
    } else {
        dist_df <- dist_df |>
            dplyr::mutate(
                Cluster = seurat_cluster[Barcode]
            )
        final_stats <- dist_df |>
            dplyr::group_by(Variant_ID, Cluster, Variant_Genotype)|>
            dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
            dplyr::group_by(Variant_ID, Cluster)
    }

    final_stats <- final_stats  |>
        dplyr::mutate(
            Total_In_Variant = sum(Count),
            Percentage = (Count / Total_In_Variant) * 100
        ) |>
        dplyr::ungroup()


    levels_genotypes <- c("WT", "HET", "HOM", "Missing/ADO")
    final_stats$Variant_Genotype <- factor(final_stats$Variant_Genotype, levels = levels_genotypes)

    return(final_stats)
}
