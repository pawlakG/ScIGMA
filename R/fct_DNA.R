
#' --------------------------------------------------------------- #
#' Filter variants and cells in a ScIGMA_object (R6)
#'
#' @description
#' This function applies quality-based filtering on variants and cells stored
#' in a `ScIGMA_object` (R6). Filtering is performed based on minimum depth (DP),
#' genotype quality (GQ), and VAF/GT consistency. Variants are retained only if
#' they are sufficiently covered across cells and observed in a minimum fraction
#' of mutated cells. Cells are retained only if they contain a sufficient number
#' of evaluable variants. The function returns a **filtered clone** of the input
#' object, leaving the original object unchanged.
#'
#' @param obj A `ScIGMA_object` (R6) with the following fields:
#'   - `meta.data`
#'   - `cell.ids`, `cell.labels`
#'   - `variants`, `variant.filter`
#'   - `vaf.mtx`, `gt.mtx`, `dp.mtx`, `gq.mtx`
#'   - `amps`, `amp.normalize.method`, `amp.mtx`
#'   - `proteins`, `protein.normalize.method`, `protein.mtx`
#'
#' @param min.dp Integer. Minimum sequencing depth (DP) per cell/variant to be considered.
#'   Default: `10`.
#' @param min.gq Integer. Minimum genotype quality (GQ) threshold.
#'   Default: `30`.
#' @param vaf.ref Numeric. Maximum tolerated VAF for reference genotypes (GT = 0).
#'   Default: `5`.
#' @param vaf.hom Numeric. Minimum expected VAF for homozygous alternate genotypes (GT = 2).
#'   Default: `95`.
#' @param vaf.het Numeric. Minimum expected VAF for heterozygous genotypes (GT = 1).
#'   Default: `35`.
#' @param min.cell.pt Numeric. Minimum percentage of cells required to cover a variant.
#'   Default: `50`.
#' @param min.mut.cell.pt Numeric. Minimum percentage of cells required to be mutated
#'   (GT = 1 or 2) for a variant to be retained.
#'   Default: `1`.
#'
#' @return A new `ScIGMA_object` (clone) containing only the filtered variants and cells.
#'   The original object is not modified.
#'
#' @examples
#' \dontrun{
#' filtered_obj <- filter_variant_ScIGMA(my_scigma_obj,
#'                                       min.dp = 15,
#'                                       min.gq = 40,
#'                                       vaf.ref = 3,
#'                                       vaf.hom = 90,
#'                                       vaf.het = 30,
#'                                       min.cell.pt = 60,
#'                                       min.mut.cell.pt = 5)
#' }
#'
#' @seealso
#' - `DelayedMatrixStats` for efficient matrix operations when using DelayedArray objects.
#' - The original `filter_variant` function from the optima package for conceptual inspiration.
#'
#' @export
filter_variant_ScIGMA <- function(
        obj,
        min.dp = 10, min.gq = 30,
        vaf.ref = 5, vaf.hom = 95, vaf.het = 35,
        min.cell.pt = 50, min.mut.cell.pt = 1
){
    # Setup multiprocessing
    bp <- if (.Platform$OS.type == "windows") SnowParam(workers = parallel::detectCores()-1)
    else MulticoreParam(workers = parallel::detectCores()-1)


    message("Filtering variant")
    stopifnot(methods::is(obj, "R6"))
    # Shortcuts
    vaf.mtx <- obj$vaf.mtx
    gt.mtx  <- obj$gt.mtx
    dp.mtx  <- obj$dp.mtx
    gq.mtx  <- obj$gq.mtx


    # 1) Quality masks / VAF vs GT consistency
    dp.tf      <- dp.mtx < min.dp
    gq.tf      <- gq.mtx < min.gq
    vaf.ref.tf <- (vaf.mtx > vaf.ref) & (gt.mtx == 0)
    vaf.hom.tf <- (vaf.mtx < vaf.hom) & (gt.mtx == 2)
    vaf.het.tf <- (vaf.mtx < vaf.het) & (gt.mtx == 1)

    keep <- !(dp.tf | gq.tf | vaf.ref.tf | vaf.hom.tf | vaf.het.tf)

    # 2) Mark unreliable genotypes as 3 (unknown), VAF set to -1
    gt.mtx[!keep] <- 3L
    vaf.mtx[gt.mtx == 3L] <- -1

    num.cells    <- nrow(gt.mtx)
    num.variants <- ncol(gt.mtx)

    # 3) Variant filtering criteria
    gt.mtx_realized <- DelayedArray::realize(gt.mtx)
    cell.cover.per.variant <- colSums2(gt.mtx_realized != 3L, BPPARAM = bp)
    mut.cells.per.variant  <- colSums2((gt.mtx_realized == 1L) | (gt.mtx_realized == 2L), BPPARAM = bp)

    cell.num.keep.tf <- cell.cover.per.variant > num.cells * (min.cell.pt/100)
    mut.cell.num.keep.tf <- mut.cells.per.variant > num.cells * (min.mut.cell.pt/100)
    variant.keep.tf <- cell.num.keep.tf & mut.cell.num.keep.tf

    # 4) Cell filtering criteria
    cell.variants.keep.tf <- rowSums2(gt.mtx_realized != 3L, BPPARAM = bp) > num.variants * (min.cell.pt/100)
    # 5) Subsetting
    keep_cells    <- which(cell.variants.keep.tf)
    keep_variants <- which(variant.keep.tf)

    # 6) Protein handling
    # proteins_are_none <- (length(obj$proteins) == 1L && identical(obj$proteins[1], "non-protein"))
    # if (proteins_are_none) {
    #     my.protein.mtx <- obj$protein.mtx
    # }
    if (length(obj$proteins) == 1) {
        my.protein.mtx <- obj$protein.mtx
    } else {
        my.protein.mtx <- obj$protein.mtx[keep_cells, , drop = FALSE]
    }
    # 7) Build a filtered clone
    filtered <- obj$clone(deep = TRUE)
    filtered$variant.filter <- "filtered"

    filtered$cell.ids.filtered     <- obj$cell.ids[keep_cells]
    filtered$cell.labels.filtered  <- obj$cell.labels[keep_cells]

    filtered$variants.filtered     <- obj$variants[variant.keep.tf]

    filtered$vaf.mtx.filtered      <- obj$vaf.mtx[keep_cells, keep_variants, drop = FALSE]
    filtered$gt.mtx.filtered       <- obj$gt.mtx [keep_cells, keep_variants, drop = FALSE]
    filtered$dp.mtx.filtered       <- obj$dp.mtx [keep_cells, keep_variants, drop = FALSE]
    filtered$gq.mtx.filtered       <- obj$gq.mtx [keep_cells, keep_variants, drop = FALSE]

    filtered$amp.mtx.filtered      <- obj$amp.mtx[keep_cells, , drop = FALSE]

    filtered$protein.mtx.filtered  <- my.protein.mtx

    # Modify variant.filter.mask according to filtered cells
    filtered$dna.variant.filter.mask.filtered <- filtered$dna.variant.filter.mask[filtered$cell.ids.filtered,]

    # Summary
    removed_cells    <- length(cell.variants.keep.tf) - sum(cell.variants.keep.tf)
    removed_variants <- length(variant.keep.tf) - sum(variant.keep.tf)
    message("Number of cells removed: ", removed_cells)
    message("Number of variants removed: ", removed_variants)

    invisible(filtered)
}

#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom DelayedMatrixStats rowSums2 colSums2
#' @importFrom DelayedArray realize
#' @importFrom BiocParallel SnowParam MulticoreParam
filter_variant_ScIGMA_mae <- function(
        obj,
        min.dp = 10, min.gq = 30,
        vaf.ref = 5, vaf.hom = 95, vaf.het = 35,
        min.cell.pt = 50, min.mut.cell.pt = 1
) {

    message("Filtering variants and cells (Matrix: Features x Cells)...")
    if (!inherits(obj, "ScIGMA_object")) stop("obj must be a ScIGMA_object.")

    bp <- if (.Platform$OS.type == "windows") {
        BiocParallel::SnowParam(workers = parallel::detectCores() - 1)
    } else {
        BiocParallel::MulticoreParam(workers = parallel::detectCores() - 1)
    }

    # ---- 1. Extraction (Variants x Cellules) ----
    vaf_mtx <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "vaf")
    gt_mtx  <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "gt")
    dp_mtx  <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "dp")
    gq_mtx  <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "gq")

    # ---- 2. Quality masks ----
    dp_tf      <- dp_mtx < min.dp
    gq_tf      <- gq_mtx < min.gq
    vaf_ref_tf <- (vaf_mtx > vaf.ref) & (gt_mtx == 0L)
    vaf_hom_tf <- (vaf_mtx < vaf.hom) & (gt_mtx == 2L)
    vaf_het_tf <- (vaf_mtx < vaf.het) & (gt_mtx == 1L)

    keep_mask <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)

    gt_mtx[!keep_mask] <- 3L
    vaf_mtx[gt_mtx == 3L] <- -1

    num_variants <- nrow(gt_mtx)
    num_cells    <- ncol(gt_mtx)

    gt_mtx_realized <- DelayedArray::realize(gt_mtx)

    cell_cover_per_variant <- DelayedMatrixStats::rowSums2(gt_mtx_realized != 3L, BPPARAM = bp)
    mut_cells_per_variant  <- DelayedMatrixStats::rowSums2((gt_mtx_realized == 1L) | (gt_mtx_realized == 2L), BPPARAM = bp)

    # sum(rowSums2((test == 1L) | (test == 2L), BPPARAM = bp) > ncol(test)*0.01)


    variant_keep_tf <- (cell_cover_per_variant > num_cells * (min.cell.pt / 100)) &
        (mut_cells_per_variant > num_cells * (min.mut.cell.pt / 100))

    cell_keep_tf <- DelayedMatrixStats::colSums2(gt_mtx_realized != 3L, BPPARAM = bp) > num_variants * (min.cell.pt / 100)

    # ---- 5. MAE Subsetting ----
    mae_filtered <- obj$mae[, cell_keep_tf, drop = FALSE]

    mae_filtered[["dna_variants"]] <- mae_filtered[["dna_variants"]][variant_keep_tf, , drop = FALSE]

    SummarizedExperiment::assay(mae_filtered[["dna_variants"]], "vaf") <- vaf_mtx[variant_keep_tf, cell_keep_tf, drop = FALSE]
    SummarizedExperiment::assay(mae_filtered[["dna_variants"]], "gt")  <- gt_mtx[variant_keep_tf, cell_keep_tf, drop = FALSE]

    # ---- 6. R6 Instantiation ----
    filtered_obj <- ScIGMA_object$new(
        mae = mae_filtered,
        mae_raw = obj$mae_raw,
        backing_files = obj$backing_files,
        filetype = obj$filetype,
        seurat_object = obj$seurat_object
    )

    removed_cells <- num_cells - sum(cell_keep_tf)
    removed_variants <- num_variants - sum(variant_keep_tf)
    message(sprintf("Filtered out: %d cells and %d variants.", removed_cells, removed_variants))

    invisible(filtered_obj)
}

# Normalize counts (amplicons x cells)
# Step 1: Each cell is normalized by its total read count (library size)
# Step 2: Each amplicon is normalized by its median value across cells
normalize_amplicon_counts <- function(count_matrix,
                                      scale_after_cell = 1,   # e.g., 1e6 for CPM, 1 for proportions
                                      epsilon = 1e-8,         # numerical stability to avoid division by zero
                                      use_nonzero_for_median = FALSE) {
    # --- Checks ---
    if (!class(count_matrix) == "DelayedMatrix") {
        stop("count_matrix must be a DelayedMatrix (amplicons x cells).")
    }
    count_matrix <- as.matrix(count_matrix)
    if (!is.numeric(count_matrix)) stop("count_matrix must be numeric.")
    if (is.null(rownames(count_matrix))) rownames(count_matrix) <- paste0("amplicon_", seq_len(nrow(count_matrix)))
    if (is.null(colnames(count_matrix))) colnames(count_matrix) <- paste0("cell_", seq_len(ncol(count_matrix)))


    message("Normalizing CNV")

    # --- Step 1: Cell-level normalization (library size) ---
    library_sizes <- colSums2(count_matrix, na.rm = TRUE)
    # avoid division by zero
    library_sizes_safe <- library_sizes + epsilon
    counts_cell_norm <- sweep(count_matrix, 2, library_sizes_safe, FUN = "/")
    # optional: scale to CPM (or other units)
    if (!is.null(scale_after_cell) && is.finite(scale_after_cell) && scale_after_cell != 1) {
        counts_cell_norm <- counts_cell_norm * scale_after_cell
    }

    # --- Step 2: Amplicon-level normalization (median across cells) ---
    if (isTRUE(use_nonzero_for_median)) {
        # compute median on non-zero values if available, otherwise use overall median
        amplicon_medians <- apply(counts_cell_norm, 1, function(x) {
            nz <- x[x > 0 & is.finite(x)]
            if (length(nz) > 0) stats::median(nz, na.rm = TRUE) else stats::median(x, na.rm = TRUE)
        })
    } else {
        amplicon_medians <- apply(counts_cell_norm, 1, stats::median, na.rm = TRUE)
    }
    amplicon_medians_safe <- amplicon_medians + epsilon

    normalized_matrix <- sweep(counts_cell_norm, 1, amplicon_medians_safe, FUN = "/")

    # --- Output ---
    attr(normalized_matrix, "library_sizes") <- library_sizes
    attr(normalized_matrix, "amplicon_medians_after_cell_norm") <- amplicon_medians
    normalized_matrix
}


#' Generate Clonal Labels from Single-Cell Genotype Matrix
#'
#' Infers clonal cluster assignments for single cells based on their mutational
#' signatures. Maps short variant IDs from the genotype matrix to full
#' nomenclature, decodes integer states (WT, HET, HOM, MISS), and aggregates
#' cells into unique clonal populations.
#'
#' @param ngt_matrix Integer matrix. Single-cell genotypes (rows = cells,
#'   columns = short variant IDs). Expected values: 0 (WT), 1 (HET), 2 (HOM),
#'   or 3/NA (MISS).
#' @param target_variants_df Data.frame. Must contain a `variant_id` column
#'   with full variant nomenclatures (e.g., used as ground truth for labels).
#'   (`MISS`) for any target variant are flagged as 'unassigned'. Default FALSE.
#'
#' @return A named list containing two elements:
#'   \item{cell_metadata}{Tibble mapping `cell_barcode` to `clonal_cluster_id`}
#'   \item{cluster_definitions}{Tibble mapping `clonal_cluster_id` to
#'   `genotype_signature`}
#'
#' @importFrom dplyr select mutate across all_of case_when filter distinct
#'   arrange left_join row_number
#' @importFrom tidyr unite
#' @importFrom stringr str_detect
#' @importFrom tibble as_tibble
#' @export
generate_clonal_labels <- function(ngt_matrix,
                                   target_variants_df,
                                   ignore_missing = FALSE) {

    # 1. Data preparation and variant mapping
    # Extract full nomenclature (ground truth for biological labels)
    full_variant_ids <- target_variants_df$variant_id
    # full_variant_ids <- rownames(target_variants_df)

    # Extract short IDs to match matrix column names
    short_variant_ids <- sub(x = full_variant_ids, pattern = "^([^:]+:)|^:", "")

    # Define mapping dictionary: short_id -> full_id
    variant_map <- setNames(full_variant_ids, short_variant_ids)

    # 2. Integrity check
    missing_variants <- setdiff(short_variant_ids, colnames(ngt_matrix))
    if (length(missing_variants) > 0) {
        stop(paste("Critical error: Missing variants in the genotype matrix:",
                   paste(missing_variants, collapse = ", ")))
    }

    # 3. State decoding function
    decode_genotype <- function(x) {
        case_when(
            x == 0 ~ "WT",
            x == 1 ~ "HET",
            x == 2 ~ "HOM",
            x == 3 | is.na(x) ~ "MISS",
            TRUE ~ "ERR"
        )
    }

    # 4. Vectorized processing
    processed_cells <- as_tibble(ngt_matrix, rownames = "cell_barcode") %>%
        select(cell_barcode, all_of(short_variant_ids)) %>%

        # Step A: Decode integers to character states
        mutate(across(all_of(short_variant_ids), decode_genotype)) %>%

        # Step B: Prefix states with full variant nomenclature via lookup map
        mutate(across(all_of(short_variant_ids),
                      ~ paste0(variant_map[cur_column()], ":", .))) %>%

        # Step C: Concatenate into a strict clonal signature string
        unite(col = "signature_string", all_of(short_variant_ids),
              sep = " & ", remove = FALSE)

    # 5. Missing data filtering
    if (ignore_missing) {
        processed_cells <- processed_cells %>%
            mutate(is_valid = !str_detect(signature_string, ":MISS"))
    } else {
        processed_cells <- processed_cells %>% mutate(is_valid = TRUE)
    }

    # 6. Build cluster lookup table (Definitions)
    cluster_definitions <- processed_cells %>%
        filter(is_valid) %>%
        group_by(signature_string) %>%
        mutate(n_cells = n()) %>%
        ungroup() %>%
        distinct(signature_string, n_cells) %>%
        arrange(desc(n_cells), signature_string) %>%
        mutate(clonal_cluster_id = sprintf("clone_%02d", row_number())) %>%
        select(clonal_cluster_id, genotype_signature = signature_string)

    # 7. Final assignment via relational join
    final_metadata <- processed_cells %>%
        left_join(cluster_definitions,
                  by = c("signature_string" = "genotype_signature")) %>%
        select(cell_barcode, clonal_cluster_id)

    # Handle unassigned cells (e.g., filtered due to missing data)
    final_metadata$clonal_cluster_id[is.na(final_metadata$clonal_cluster_id)] <-
        "unassigned"

    return(list(
        cell_metadata = final_metadata,
        cluster_definitions = cluster_definitions
    ))
}


#' Generate Heatmap of DNA Variant Genotypes
#'
#' This function processes genotype data for a set of selected DNA variants,
#' performs clonal clustering on samples without missing genotype information,
#' and generates a complex heatmap visualizing the genotypes across all samples,
#' organized by the inferred clonal clusters.
#'
#' @param obj A list or S4 object (e.g., a SingleCellExperiment or similar custom object)
#'   containing the necessary genotype information. It must contain:
#'   \itemize{
#'     \item `$gt.mtx`: A matrix of genotypes (e.g., 0, 1, 2 for WT, HET, HOM)
#'       where rows are samples (cells) and columns are variant IDs.
#'     \item `$dna.variant.filter.mask.filtered`: A boolean/integer matrix (same dimensions as `$gt.mtx`)
#'       indicating which genotypes should be considered as masked/filtered (e.g., unreliable).
#'   }
#' @param selected_variants_df A data frame containing information about the variants
#'   to be included in the heatmap. It must have a column named `variant_id`
#'   containing the full variant identifiers.
#' @param min_prop_cluster A numeric value (default: 0.01). The minimum proportion
#'   of total samples required for a cluster to be considered a major cluster.
#'   Clusters smaller than this proportion are grouped into a single "small" cluster.
#' @param heatmap_include_all_samples A logical value (default: TRUE). If TRUE,
#'   samples with missing/filtered variant data are included in the heatmap
#'   under a separate "missing" split. If FALSE, only samples with complete,
#'   unfiltered genotype data for the selected variants are included, and the
#'   "small" cluster is also excluded from the final heatmap.
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item `heatmap`: A `Heatmap` object generated by the `ComplexHeatmap` package.
#'     \item `clones`: A factor vector of clonal cluster assignments for the
#'       samples that had complete genotype information (used for clustering).
#'   }
#' @importFrom magrittr "%>%"
#' @importFrom forcats fct_c fct_infreq
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom colorBlindness paletteMartin
#' @importFrom stringr str_detect
#' @export
#'
#' @examples
#' # The example below is illustrative and will not run without the necessary
#' # input objects and the `generate_clonal_labels` function.
#' #
#' # obj <- list(
#' #   gt.mtx = matrix(sample(0:2, 100*5, replace = TRUE), nrow = 100,
#' #                   dimnames = list(paste0("Sample", 1:100), paste0("V", 1:5))),
#' #   dna.variant.filter.mask.filtered = matrix(sample(0:1, 100*5, replace = TRUE, prob=c(0.9, 0.1)), nrow = 100,
#' #                   dimnames = list(paste0("Sample", 1:100), paste0("V", 1:5)))
#' # )
#' # selected_variants_df <- data.frame(variant_id = paste0("ChrX:", paste0("V", 1:5)))
#' #
#' # # Assuming 'generate_clonal_labels' is defined elsewhere:
#' # # heatmap_result <- generate_dna_variant_heatmap(obj, selected_variants_df)
#' # # draw(heatmap_result$heatmap)
generate_dna_variant_heatmap <- function(obj,
                                         selected_variants_df,
                                         min_prop_cluster = 0.01,
                                         heatmap_include_all_samples = TRUE,
                                         use_imputed = FALSE) { # <-- NEW ARGUMENT

    target_variants <- selected_variants_df$variant_id

    short_variants <- sub(x = target_variants, pattern = "^([^:]+:)|^:", "")

    if (isTRUE(use_imputed)) {

        if (!"compass_imputed" %in% SummarizedExperiment::assayNames(obj$mae[["dna_variants"]])) {
            stop("Error : 'compass_imputed' Assay not found. Please rerun COMPASS inference.")
        }
        gt_full <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "compass_imputed")

        if (!all(short_variants %in% rownames(gt_full))) {
            stop("Error : Some selected variants do not exist in the matrix.")
        }

        msk_full <- matrix(0L, nrow = nrow(gt_full), ncol = ncol(gt_full),
                           dimnames = dimnames(gt_full))

    } else {
        # Extraction Out-of-Core native (Raw data)
        gt_full <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "gt")
        msk_full <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "variant_filter_mask")
    }

    gt <- t(as.matrix(gt_full[short_variants, , drop = FALSE]))
    msk <- t(as.matrix(msk_full[short_variants, , drop = FALSE])) != 0

    colnames(gt) <- short_variants
    colnames(msk) <- short_variants
    selected_variants <- short_variants

    common_rows <- intersect(rownames(gt), rownames(msk))
    common_cols <- intersect(colnames(gt), colnames(msk))

    if (length(common_rows) == 0L || length(common_cols) == 0L) {
        stop("Fatal: No overlap between gt matrix and mask.")
    }

    # 5. Application du NGT_MASK
    gt_filtered <- gt[common_rows, common_cols, drop = FALSE]

    # Set value 3 for filtered genotypes
    gt_filtered[cbind(row(msk)[msk], col(msk)[msk])] <- 3

    tmp_heamtap_matrix_filtered_noMissing <- gt_filtered[rowSums(gt_filtered == 3) == 0, , drop = FALSE]
    tmp_heamtap_matrix_filtered_withMissing <- gt_filtered[rowSums(gt_filtered == 3) > 0, , drop = FALSE]

    if (nrow(tmp_heamtap_matrix_filtered_withMissing) > 0) {
        tmp_heamtap_matrix_filtered_withMissing[tmp_heamtap_matrix_filtered_withMissing == 3] <- NA
    }

    results_clustering <- generate_clonal_labels(
        ngt_matrix = tmp_heamtap_matrix_filtered_noMissing,
        target_variants_df = selected_variants_df,
        ignore_missing = TRUE
    )

    clustered_samples <- setNames(results_clustering$cell_metadata$clonal_cluster_id,
                                  nm = results_clustering$cell_metadata$cell_barcode)

    ### Reassign too small samples clusters
    res_table_clusters <- table(clustered_samples)
    too_small_clusters <- names(res_table_clusters)[res_table_clusters / nrow(gt) < min_prop_cluster]
    clustered_samples[clustered_samples %in% too_small_clusters] <- "small"

    small_cluster <- clustered_samples[clustered_samples == "small"] |> as.factor()
    nonSmall_cluster <- clustered_samples[clustered_samples != "small"] |> sort() |> as.factor()
    nonSmall_cluster <- forcats::fct_infreq(nonSmall_cluster)
    levels(nonSmall_cluster) <- sprintf("clone_%02d", 1:length(levels(nonSmall_cluster)))

    clustered_samples <- forcats::fct_c(nonSmall_cluster, small_cluster)
    names(clustered_samples) <- c(names(nonSmall_cluster), names(small_cluster))

    clustered_samples <- forcats::fct_infreq(clustered_samples)

    if (!is.null(obj$dna_clones_renamed)) {
        clustered_samples <- obj$dna_clones_renamed
    }

    ### Reorder
    tmp_heamtap_matrix_filtered_noMissing_ordered <- tmp_heamtap_matrix_filtered_noMissing[names(sort(clustered_samples)), , drop = FALSE]

    sorted_clusters <- sort(clustered_samples)
    desired_levels <- levels(sorted_clusters)

    # rbind with samples with missing info
    if (heatmap_include_all_samples) {
        tmp_heamtap_matrix_filtered_complete_ordered <- rbind(tmp_heamtap_matrix_filtered_noMissing_ordered,
                                                              tmp_heamtap_matrix_filtered_withMissing)

        heatmap_split_vector <- c(as.character(sorted_clusters),
                                  rep("missing", nrow(tmp_heamtap_matrix_filtered_withMissing)))

        heatmap_split_vector <- factor(heatmap_split_vector,
                                       levels = c(desired_levels, "missing"))

    } else {
        tmp_heamtap_matrix_filtered_complete_ordered <- tmp_heamtap_matrix_filtered_noMissing_ordered

        heatmap_split_char <- as.character(sorted_clusters)
        keep_idx <- heatmap_split_char != "small"

        tmp_heamtap_matrix_filtered_complete_ordered <- tmp_heamtap_matrix_filtered_complete_ordered[keep_idx, , drop = FALSE]
        heatmap_split_char <- heatmap_split_char[keep_idx]

        final_levels <- desired_levels[desired_levels != "small"]
        heatmap_split_vector <- factor(heatmap_split_char, levels = final_levels)
    }

    # plot heatmap
    dna_variant_colorPalette <- setNames(c("#E9E8EC", "#BAB7D0", "#3C2692"), nm = c("0", "1", "2"))

    heatmap_true_levels <- levels(clustered_samples)[levels(clustered_samples) != "small"]
    annotationColor <- list(Cluster = setNames(c(colorBlindness::paletteMartin[-1][1:length(heatmap_true_levels)], "grey", "#333333"),
                                               nm = c(heatmap_true_levels, "missing", "small")))

    dna_variant_annotation <- ComplexHeatmap::rowAnnotation(
        Cluster = heatmap_split_vector,
        col = annotationColor,
        show_annotation_name = FALSE,
        show_legend = FALSE
    )

    # Verify colnames before plotting
    idx_colnames <- sapply(colnames(tmp_heamtap_matrix_filtered_complete_ordered), function(x) {
        which(grepl(x, selected_variants_df$variant_id))
    }, simplify = FALSE, USE.NAMES = FALSE) |> unlist()

    new_colnames <- paste0(selected_variants_df$gene[idx_colnames], "   \n", selected_variants_df$cdna[idx_colnames])

    ## Heatmap
    heatmap <- ComplexHeatmap::Heatmap(
        tmp_heamtap_matrix_filtered_complete_ordered,
        row_order = rownames(tmp_heamtap_matrix_filtered_complete_ordered),
        row_title_rot = 0,
        column_names_rot = 45,
        row_split = heatmap_split_vector,
        show_column_dend = FALSE,
        column_split = selected_variants,
        column_title = NULL,
        column_labels = new_colnames,
        cluster_rows = FALSE,
        show_row_names = FALSE,
        na_col = "black",
        col = dna_variant_colorPalette,
        show_heatmap_legend = TRUE,
        left_annotation = dna_variant_annotation,
        heatmap_legend_param = list(
            title = "Genotype",
            at = c("0", "1", "2"),
            labels = c("WT", "HET", "HOM"),
            grid_height = grid::unit(4, "cm")
        )
    )

    return(list("heatmap" = heatmap,
                "clones" = clustered_samples))
}

# NEW

#'
#' @description
#' Extracts and formats data from a ScIGMA_object MAE to feed the COMPASS MCMC backend.
#' Computes M_ref and M_alt from Depth and VAF. Aggregates amplicon depths by Gene to form
#' the region count matrix (C). Computes the 0-based topological mapping vector (locus_regions).
#'
#' @param obj A ScIGMA_object (R6) containing the MAE.
#' @param selected_variants Character vector. The exact variant IDs to include in the phylogeny.
#'
#' @return A list containing:
#'   - M_ref: Matrix (Cells x Variants) of Reference allele counts.
#'   - M_alt: Matrix (Cells x Variants) of Alternate allele counts.
#'   - C: Matrix (Cells x Genes) of total region depth.
#'   - locus_regions: Integer vector (0-based) mapping each variant to its column in C.
#'   - regions_names: Character vector of the region (Gene) names in order.
#'
#' @export
build_compass_matrices <- function(obj, selected_variants) {

    message("Building COMPASS topological and count matrices...")

    # Extraction depuis le MAE (Variants x Cellules)
    dp_full <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "dp")
    vaf_full <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "vaf")

    dp_sub <- as.matrix(dp_full[selected_variants, , drop = FALSE])
    vaf_sub <- as.matrix(vaf_full[selected_variants, , drop = FALSE])

    dp_mat <- t(dp_sub)
    vaf_mat <- t(vaf_sub)

    M_alt <- round(dp_mat * (vaf_mat / 100))
    M_ref <- dp_mat - M_alt

    missing_mask <- vaf_mat < 0
    M_alt[missing_mask] <- NA
    M_ref[missing_mask] <- NA

    amp_full <- SummarizedExperiment::assay(obj$mae[["amplicons"]], "counts")
    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(obj$mae[["amplicons"]]))

    genome_v <- S4Vectors::metadata(obj$mae)$genome_version
    if (is.null(genome_v)) genome_v <- "hg19"

    cnv_annotated <- annotate_genomic_regions(cnv_id_table, build = genome_v)
    cnv_annotated$symbol[is.na(cnv_annotated$symbol)] <- "Unknown"

    cnv_annotated$compass_region <- paste0(cnv_annotated$chrom, "_", cnv_annotated$symbol)

    amp_mat <- t(as.matrix(amp_full))

    C_mat <- aggregate_matrix_by_mappingTable(
        numeric_matrix = amp_mat,
        mapping_table = cnv_annotated,
        feature_column_name = dna_id,
        group_column_name = compass_region # <-- CHANGED
    )

    # ---- 3. Mapping Topologique (locus_regions) ----
    dna_id_table <- as.data.frame(SummarizedExperiment::rowData(obj$mae[["dna_variants"]]))
    snv_info <- dna_id_table[selected_variants, , drop = FALSE]

    tmp_x <- data.frame(variant_id = rownames(snv_info), amplicon = snv_info$amplicon)
    tmp_y <- cnv_annotated[, c("dna_id", "compass_region")]

    snv_to_gene <- merge(
        x = data.frame(variant_id = rownames(snv_info), amplicon = snv_info$amplicon),
        y = cnv_annotated[, c("dna_id", "compass_region")],
        by.x = "amplicon",
        by.y = "dna_id",
        all.x = TRUE,
        sort = FALSE
    )
    snv_to_gene <- tmp_x

    snv_to_gene$compass_region <- tmp_y$compass_region[match(tmp_x$amplicon, tmp_y$dna_id)]


    snv_to_gene <- snv_to_gene[match(selected_variants, snv_to_gene$variant_id), ]
    # snv_to_gene$compass_region <- paste(snv_to_gene$)
    snv_to_gene$compass_region[is.na(snv_to_gene$compass_region)] <- "0_Unknown"

    gene_cols <- colnames(C_mat)
    locus_regions <- match(snv_to_gene$compass_region, gene_cols) - 1L

    if (any(is.na(locus_regions))) {
        stop("Fatal: Some selected variants could not be mapped to a valid region in matrix C.")
    }

    message("Matrices successfully built for COMPASS.")

    return(list(
        M_ref = M_ref,
        M_alt = M_alt,
        C = C_mat,
        locus_regions = locus_regions,
        variants = selected_variants,
        regions = gene_cols
    ))
}

# NEW
# File: R/fct_compass_wrapper.R

#' Run COMPASS MCMC phylogeny inference
#'
#' @description
#' Triggers the COMPASS CLI via system call. Isolates MCMC memory overhead
#' from the main R process.
#'
#' @param output_dir Character. Path for COMPASS output.
#' @param compass_exec Character. Path or alias to COMPASS executable.
#' @param n_iters Integer. Number of MCMC iterations per chain.
#' @param restarts Integer. Number of random initializations.
#' @param threads Integer. Cores allocated for parallel chains.
#'
#' @return Logical. TRUE if successful.
#' @export
run_compass_mcmc <- function(input_dir, output_dir, compass_exec = "compass",
                             n_iters = 100000, restarts = 10, threads = 8) {
    if (!dir.exists(input_dir)) {
        stop("Fatal: input directory not found.")
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    args <- c(
        "--m-ref", file.path(input_dir, "M_ref.tsv"),
        "--m-alt", file.path(input_dir, "M_alt.tsv"),
        "--c-mat", file.path(input_dir, "C_mat.tsv"),
        "--locus-regions", file.path(input_dir, "locus_regions.txt"),
        "--out-dir", output_dir,
        "--iters", as.character(n_iters),
        "--restarts", as.character(restarts),
        "--threads", as.character(threads)
    )

    message("Spawning COMPASS backend process...")

    exit_status <- system2(
        command = compass_exec,
        args = args,
        stdout = file.path(output_dir, "compass.log"),
        stderr = file.path(output_dir, "compass.err")
    )

    if (exit_status != 0) {
        stop("Critical: MCMC process crashed. Inspect compass.err.")
    }

    invisible(TRUE)
}

#'
#' doublet cells globally from the MultiAssayExperiment object. Results are
#' stored non-destructively in the MAE metadata slot.
#'
#' @param scigma_data List containing a MultiAssayExperiment ('mae').
#' @param target_variants Character vector. Pathogenic variants to track.
#' @param chain_length Integer. Number of MCMC iterations.
#' @param output_dir Character. Directory for temporary files. Defaults to
#'   tempdir() for strict BioConductor compliance.
#' @return Updated scigma_data object containing only singlet cells and
#'   COMPASS results in MAE metadata.
#' @export
infer_clonal_architecture <- function(scigma_data, target_variants,
                                      chain_length = 500L,
                                      output_dir = tempdir()) {

    # 1. Matrix extraction
    compass_inputs <- build_compass_matrices(
        obj = scigma_data,
        selected_variants = target_variants
    )

    mat_ref <- t(as.matrix(compass_inputs$M_ref))
    storage.mode(mat_ref) <- "integer"

    mat_alt <- t(as.matrix(compass_inputs$M_alt))
    storage.mode(mat_alt) <- "integer"

    mat_cna <- t(as.matrix(compass_inputs$C))
    storage.mode(mat_cna) <- "integer"

    gt_assay <- SummarizedExperiment::assay(
        scigma_data$mae[["dna_variants"]], "gt"
    )
    # GT extraction already uses the correct order, but ensured here for consistency
    mat_gt <- as.matrix(gt_assay[target_variants, , drop = FALSE])
    mat_gt[mat_gt == 3L] <- NA
    storage.mode(mat_gt) <- "integer"

    variant_matrices <- list(REF = mat_ref, ALT = mat_alt, GT = mat_gt)

    # 2. Metadata extraction
    dna_se <- scigma_data$mae[["dna_variants"]]
    snv_row_data <- as.data.frame(SummarizedExperiment::rowData(dna_se))
    snv_sub <- snv_row_data[target_variants, ]

    vec_locus_names <- snv_sub$gene
    vec_locus_chrom <- snv_sub$chrom

    amp_se <- scigma_data$mae[["amplicons"]]
    cna_row_data <- as.data.frame(SummarizedExperiment::rowData(amp_se))

    get_gene <- function(x) strsplit(x, "_")[[1]][3]
    vec_region_names <- paste0(cna_row_data$chrom, "_",
                               sapply(cna_row_data$dna_id, get_gene))
    vec_region_names <- unique(vec_region_names)

    get_chrom <- function(x) strsplit(x, "_")[[1]][1]
    vec_region_chrom <- sapply(vec_region_names, get_chrom, USE.NAMES = FALSE)
    vec_region_chrom <- sub("^chr", "", vec_region_chrom, ignore.case = TRUE)

    # 3. BioConductor compliant I/O
    prefix_out <- file.path(
        output_dir, paste0("compass_", as.integer(Sys.time()))
    )

    use_cna <- if (ncol(variant_matrices$REF) != ncol(mat_cna)) FALSE else TRUE

    # 4. C++ Execution
    run_compass_mcmc(
        variant_matrices   = variant_matrices,
        locus_regions      = compass_inputs$locus_regions,
        region_matrix      = mat_cna,
        output_prefix      = prefix_out,
        locus_names        = vec_locus_names,
        locus_chromosomes  = vec_locus_chrom,
        region_names       = vec_region_names,
        region_chromosomes = vec_region_chrom,
        chains             = 4L,
        chain_length       = as.integer(chain_length),
        patient_sex        = "female",
        use_cna            = use_cna
    )

    # 5. Imputation and global filtering
    mat_imputed <- get_imputed_genotypes(prefix_out = prefix_out)
    cells_to_keep <- rownames(mat_imputed)

    if (length(cells_to_keep) == 0) {
        stop("Fatal: No singlet cells returned by COMPASS inference.")
    }

    scigma_data$mae <- scigma_data$mae[, cells_to_keep, ]

    # Matrix alignment
    mat_imputed_t <- t(mat_imputed)
    rownames(mat_imputed_t) <- target_variants
    colnames(mat_imputed_t) <- cells_to_keep

    # 6. Non-Destructive Storage in S4 metadata
    S4Vectors::metadata(scigma_data$mae)$compass <- list(
        imputed_gt = mat_imputed_t,
        singlet_barcodes = cells_to_keep,
        target_variants = target_variants
    )

    return(scigma_data)
}

#' Generate stable color palette for DNA clones
#' @param clone_levels Character vector or factor of clone names
#' @return Named character vector of hex colors
#' @importFrom viridisLite turbo
generate_clone_palette <- function(clone_levels) {
    unique_clones <- unique(as.character(clone_levels))
    unique_clones <- unique_clones[!is.na(unique_clones)]

    missing_tags <- c("Missing", "Missing/ADO", "NA", "Unknown", "missing")
    is_missing <- unique_clones %in% missing_tags

    base_clones <- sort(unique_clones[!is_missing])

    clone_colors <- c()
    if (length(base_clones) > 0) {
        clone_colors <- viridisLite::turbo(n = length(base_clones))
        names(clone_colors) <- base_clones
    }

    missing_colors <- c("Missing" = "#e0e0e0")
    clone_colors <- c(clone_colors, missing_colors)
    small_color <- c("small" = "#00000019")
    clone_colors <- c(clone_colors, small_color)

    return(clone_colors)
}
