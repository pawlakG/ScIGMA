
#' Aggregate Matrix Columns by Amplicon Mean
#'
#' This function aggregates columns of a numeric matrix (e.g., variant data)
#' into new columns based on mean values, using a provided mapping table (dna_id to amplicon).
#' It uses an optimized matrix multiplication approach, avoiding slow wide-to-long transformations
#' and row-wise loops, which is essential for large genomic datasets.
#'
#' @param numeric_matrix A numeric matrix where rows are samples and columns are
#'   the features to be aggregated (e.g., 'dna_id').
#' @param mapping_table A data frame or tibble with at least two columns:
#'   one for the original feature identifiers (matching `colnames(numeric_matrix)`)
#'   and one for the grouping identifier (e.g., 'amplicon').
#' @param feature_column_name The unquoted name of the column in `mapping_table`
#'   that contains the original feature identifiers (e.g., dna_id).
#' @param group_column_name The unquoted name of the column in `mapping_table`
#'   that contains the grouping identifiers (e.g., amplicon).
#'
#' @return A numeric matrix with the same number of rows (samples) but with
#'   columns aggregated by the mean of the specified groups.
#' @importFrom dplyr select filter pull rename_with
#' @importFrom Matrix sparse.model.matrix
#'
#' @examples
#' # Simulate data (5 samples, 10 DNA IDs, 3 Amplicons)
#' dna_ids_example <- paste0("chr1:115256", 568:577, ":C/T")
#' data_matrix <- matrix(
#'   runif(5 * 10, 0, 1),
#'   nrow = 5,
#'   ncol = 10,
#'   dimnames = list(paste0("Sample_", 1:5), dna_ids_example)
#' )
#'
#' mapping <- tibble(
#'   variant_id = dna_ids_example,
#'   amplicon_group = c(
#'     rep("NRAS_1", 3), rep("KIT_2", 4), rep("FLT3_3", 3)
#'   )
#' )
#'
#' aggregated_data <- aggregate_matrix_by_amplicon(
#'   numeric_matrix = data_matrix,
#'   mapping_table = mapping,
#'   feature_column_name = variant_id,
#'   group_column_name = amplicon_group
#' )
#' print(aggregated_data)
#'
aggregate_matrix_by_mappingTable<- function(numeric_matrix,
                                            mapping_table,
                                            feature_column_name,
                                            group_column_name) {

    # Ensure necessary packages are loaded
    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package 'Matrix' is required for performance. Please install it with install.packages('Matrix').")
    }

    # --- Input Validation and Alignment ---

    # Capture column names passed by the user as characters
    feature_col <- deparse(substitute(feature_column_name))
    group_col <- deparse(substitute(group_column_name))

    # Check if all matrix columns are present in the mapping table
    matrix_columns <- colnames(numeric_matrix)

    # Filter mapping table to include only the features present in the matrix
    # This handles the realistic scenario where the mapping table is larger than the matrix
    mapping_filtered <- mapping_table %>%
        dplyr::filter(!!rlang::sym(feature_col) %in% matrix_columns) %>%
        dplyr::select(!!rlang::sym(feature_col), !!rlang::sym(group_col))

    # Check for missing features
    missing_features <- setdiff(matrix_columns, mapping_filtered %>%
                                    dplyr::pull(!!rlang::sym(feature_col)))
    if (length(missing_features) > 0) {
        warning(paste(length(missing_features), "features in the matrix were not found in the mapping table and will be dropped:",
                      paste(head(missing_features, 5), collapse = ", ")))
        # Drop missing features from the matrix
        numeric_matrix <- numeric_matrix[, !matrix_columns %in% missing_features]
        matrix_columns <- colnames(numeric_matrix) # Update columns
    }

    # Align the mapping table to the exact order of the matrix columns
    mapping_aligned <- mapping_filtered[match(matrix_columns, mapping_filtered %>%
                                                  dplyr::pull(!!rlang::sym(feature_col))), ]

    # --- Core Matrix Aggregation ---

    # 1. Create the Grouping/Design Matrix (G)
    # Convert the group column to a factor for proper matrix generation
    group_factor <- factor(mapping_aligned %>% dplyr::pull(!!rlang::sym(group_col)))

    # sparse.model.matrix is the fastest way to create the binary indicator matrix G
    # G has dimensions: [N_features x N_groups]
    group_matrix <- Matrix::sparse.model.matrix(~ 0 + group_factor)

    # Rename columns to be the actual amplicon names
    colnames(group_matrix) <- gsub("group_factor", "", colnames(group_matrix))

    # 2. Calculate the Sums per Group (M x G)
    # Note: The alignment above ensures that numeric_matrix %*% group_matrix is correct
    # [N_samples x N_features] %*% [N_features x N_groups] = [N_samples x N_groups (Sums)]
    summed_matrix <- numeric_matrix %*% group_matrix

    # 3. Calculate the Denominators (Number of features per group)
    # This is the normalization vector to convert sums to means
    denominator_vector <- Matrix::colSums(group_matrix)

    # 4. Calculate the Mean (Division by Denominator)
    # Use sweep or standard vectorized division for speed
    aggregated_matrix <- summed_matrix / rep(denominator_vector, each = nrow(summed_matrix))

    return(aggregated_matrix)
}

#' Filter CNV Matrix Based on Amplicon Read Depth and Completeness
#'
#' This function filters a Copy Number Variation (CNV) matrix (`obj$amp.mtx`)
#' based on two main criteria applied to an associated read depth matrix (`obj$dp.mtx`):
#' 1. **Amplicon Completeness**: Filters cells where the proportion of amplicons
#' with read depth greater than a threshold (`amp_readDepth`) is below a specified
#' percentage (`amp_completeness`).
#' 2. **Mean Cell Read Depth**: Filters rows (amplicons) from the resulting CNV
#' matrix where the mean read depth across the selected cells is below a
#' threshold (`amp_meanCellRead`).
#'
#' It requires an auxiliary function `aggregate_matrix_by_amplicon` to first
#' aggregate the depth matrix by amplicon identifiers.
#'
#' @param obj A list or environment containing the required matrices and table:
#'   \itemize{
#'     \item{\strong{amp.mtx}: A numeric matrix (CNV matrix) with amplicons as columns
#'     and cells as rows (or vice-versa, but the current logic suggests amplicons as columns
#'     for filtering by column name).}
#'     \item{\strong{dp.mtx}: A numeric matrix (read depth matrix) with IDs
#'     (variant positions) as columns and cells as rows.}
#'     \item{\strong{dna_id_table}: A data frame or tibble used for mapping IDs
#'     to amplicons, required by \code{aggregate_matrix_by_amplicon}.}
#'   }
#' @param amp_completeness A numeric value (0-100) representing the minimum
#'   percentage of amplicons that must have a read depth greater than \code{amp_readDepth}
#'   for a cell to be kept. Default is 80.
#' @param amp_readDepth A numeric value representing the minimum read depth
#'   required for an amplicon to be considered "complete" in a cell. Default is 10.
#' @param amp_meanCellRead A numeric value representing the minimum mean read depth
#'   required for an amplicon (across selected cells) to be kept in the final
#'   CNV matrix. Default is 10.
#' @param aggregate_matrix_by_amplicon A function that aggregates a numeric matrix
#'   (like \code{obj$dp.mtx}) by a grouping column (like 'amplicon') using a
#'   mapping table. This function is expected to take the following arguments:
#'   \code{numeric_matrix}, \code{mapping_table}, \code{feature_column_name},
#'   \code{group_column_name}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item{\strong{filtered_cnv_mtx}: The CNV matrix (\code{obj$amp.mtx})
#'       filtered by both cell completeness and amplicon mean read depth.}
#'     \item{\strong{cell_selection}: A character vector of the cell IDs that
#'       passed the amplicon completeness filter.}
#'     \item{\strong{amplicon_meanCellRead_summary}: A table summarizing the
#'       number of amplicons that passed/failed the mean read depth filter.}
#'   }
#' @export
#'
#' @examples
#' # The 'aggregate_matrix_by_amplicon' function is not defined here.
#' # In a real package, you would ensure it is available or mock it for examples.
#' \dontrun{
#' # Assuming 'obj' is a list with 'amp.mtx', 'dp.mtx', 'dna_id_table' and
#' # 'aggregate_matrix_by_amplicon' is a defined function.
#' filtered_data <- filter_cnv_matrix_by_completeness(
#'   obj = my_obj,
#'   amp_completeness = 80,
#'   amp_readDepth = 10,
#'   amp_meanCellRead = 10,
#'   aggregate_matrix_by_amplicon = my_aggregate_func
#' )
#' head(filtered_data$filtered_cnv_mtx)
#' }
filter_cnv_matrix_by_completeness <- function(
        obj,
        amp_completeness = 80,
        amp_readDepth = 10,
        amp_meanCellRead = 10,
        aggregate_fun = aggregate_matrix_by_mappingTable
) {
    # 1. Extraction BioConductor (Lignes = Features, Colonnes = Cellules)
    amp_mtx <- SummarizedExperiment::assay(obj$mae[["amplicons"]], "counts")
    dp_mtx <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "dp")
    dna_id_table <- as.data.frame(SummarizedExperiment::rowData(obj$mae[["dna_variants"]]))

    # 2. Transposition pour l'agrégation
    # aggregate_fun attend : Lignes = Samples (Cells), Colonnes = Features (Variants)
    tmp_dp_mtx <- t(as.matrix(dp_mtx))

    tmp_dp_mtx_aggregated <- aggregate_fun(
        numeric_matrix = tmp_dp_mtx,
        mapping_table = dna_id_table,
        feature_column_name = dna_id,
        group_column_name = amplicon
    ) |> as.matrix()

    tmp_dp_mtx_aggregated <- tmp_dp_mtx_aggregated[, colnames(tmp_dp_mtx_aggregated) != "NA", drop = FALSE]

    ## --- Filter according to amplicon completeness (Cell-level filter) ---
    # CORRECTION CRITIQUE : tmp_dp_mtx_aggregated est [Cells x Amplicons].
    # Pour évaluer chaque cellule, on itère sur les LIGNES (MARGIN = 1).
    cell_ampCompleteness_filter <- apply(
        X = tmp_dp_mtx_aggregated,
        MARGIN = 1,
        FUN = function(x) sum(x > amp_readDepth) / length(x)
    ) >= (amp_completeness / 100)

    cell_ampCompleteness_selected <- names(cell_ampCompleteness_filter)[cell_ampCompleteness_filter]

    if (length(cell_ampCompleteness_selected) == 0) {
        warning("No cells passed the amplicon completeness filter. Returning empty results.")
        return(list(
            filtered_cnv_mtx = amp_mtx[0, 0, drop = FALSE],
            cell_selection = character(0),
            amplicon_meanCellRead_summary = NULL
        ))
    }

    # amp_mtx est [Amplicons x Cells]. On filtre les colonnes.
    tmp_cnv_mtx <- amp_mtx[, cell_ampCompleteness_selected, drop = FALSE]

    ## --- Filter according to minimum mean read depth per amplicon ---
    amplicon_meanCellRead_filter <- rowMeans(tmp_cnv_mtx, na.rm = TRUE) >= amp_meanCellRead
    meanCellRead_summary <- table(amplicon_meanCellRead_filter)

    filtered_cnv_mtx <- tmp_cnv_mtx[amplicon_meanCellRead_filter, , drop = FALSE]

    return(list(
        filtered_cnv_mtx = filtered_cnv_mtx,
        cell_selection = cell_ampCompleteness_selected,
        amplicon_selection = names(amplicon_meanCellRead_filter)[amplicon_meanCellRead_filter],
        amplicon_meanCellRead_summary = meanCellRead_summary
    ))
}


#' @title Calculate Clonal CNV Profiles
#' @description Filters the CNV matrix, normalizes read counts, aggregates
#'   the normalized data by clone (using median), and calculates the final
#'   profile relative to a specified diploid reference clone.
#'
#' @param obj The primary data object containing the CNV matrix
#'   (expected to be used by filter_cnv_matrix_by_completeness).
#' @param dna_variant_clones A dataframe or list containing cell-to-clone assignments,
#'   specifically the 'clones' vector/factor.
#' @param diploid_ref The ID (string) of the clone group to be used as the 2N reference.
#' @param amp_completeness Minimum amplicon completeness (percentage).
#' @param amp_readDepth Minimum amplicon read depth.
#' @param amp_meanCellRead Minimum mean cell read depth.
#' @param exclude_clone A character string of a clone to exclude from the final profile
#'   (e.g., "small" or NA).
#'
#' @export
#'
#' @return A matrix where rows are genomic regions (amplicons) and columns are
#'   clonal profiles, normalized to the diploid reference.
#'
filter_cnv_profile <- function(obj,
                               dna_variant_clones,
                               amp_completeness = 80,
                               amp_readDepth = 10,
                               amp_meanCellRead = 10) {

    # Check if the filtering function exists (assuming it's external)
    if (!exists("filter_cnv_matrix_by_completeness")) {
        stop("Required function 'filter_cnv_matrix_by_completeness' not found.")
    }

    ## 1. Filter the raw data based on quality metrics
    filtered_data <- filter_cnv_matrix_by_completeness(
        obj = obj,
        amp_completeness = amp_completeness,
        amp_readDepth = amp_readDepth,
        amp_meanCellRead = amp_meanCellRead
    )

    return(filtered_data)
}


#' @title Process CNV Matrix to Diploid-Referenced Clonal Profiles
#' @description Normalizes read counts, aggregates the normalized data by clone
#'   (using median), and calculates the final profile relative to a specified
#'   diploid reference clone.
#' @import tibble
#'
#' @param filtered_data A list containing the CNV matrix under
#'   `$filtered_cnv_mtx` (amplicons x cells).
#' @param dna_variant_clones A named vector or factor where names are cell barcodes
#'   and values are clone IDs.
#' @param diploid_ref The ID (string) of the clone group to be used as the 2N reference.
#' @param exclude_clone A character string of a clone to exclude from the final profile
#'   (e.g., "small" or NA).
#'
#' @return A matrix where rows are genomic regions (amplicons) and columns are
#'   clonal profiles, normalized to the diploid reference.
#'
process_cnv_to_clonal_profile <- function(filtered_data,
                                          dna_variant_clones,
                                          diploid_ref = "2",
                                          exclude_clone = "small") {

    tmp_clones <- dna_variant_clones

    # --- INPUT CHECKS ---
    if (is.null(names(tmp_clones))) {
        warning("Clone assignments have no names/barcodes (cell IDs). Proceeding, but matching may fail.")
    }

    if (is.null(filtered_data$filtered_cnv_mtx)) {
        stop("filtered_data must contain '$filtered_cnv_mtx'.")
    }

    ## 1. Prepare and normalize the filtered read counts
    # Transpose the matrix to be cells x amplicons for easier column-wise normalization
    # tmp_dna_reandCounts_mtx <- t(as.matrix(filtered_data$filtered_cnv_mtx))
    tmp_dna_reandCounts_mtx <- as.matrix(filtered_data$filtered_cnv_mtx)

    # Keep only cells that were assigned to a clone (clusterised cells)
    cell_barcodes_to_keep <- intersect(colnames(tmp_dna_reandCounts_mtx), names(tmp_clones))

    if (length(cell_barcodes_to_keep) == 0) {
        stop("No overlapping cell barcodes found between CNV matrix and clone assignments.")
    }

    tmp_dna_reandCounts_mtx <- tmp_dna_reandCounts_mtx[, cell_barcodes_to_keep, drop = FALSE]

    # Calculate normalization factor (Total reads across all filtered cells)
    total_reads <- sum(tmp_dna_reandCounts_mtx, na.rm = TRUE)

    # Perform normalization (Counts per Million scaled by Total Reads)
    tmp_dna_reandCounts_mtx_norm <- apply(tmp_dna_reandCounts_mtx,
                                          2, # Apply across columns (cells)
                                          function(x) 10^6 * x / total_reads)

    ## 2. Aggregate normalized data by clone using the median
    clonal_profile <- tmp_dna_reandCounts_mtx_norm |>
        t() |> # Transpose back to amplicons x cells
        as_tibble(rownames = "cell_barcode") |>
        mutate(clone_id = tmp_clones[cell_barcode]) |>
        # Filter out cells with NA clone assignments
        filter(!is.na(clone_id)) |>
        group_by(clone_id) |>
        # Use MEDIAN for robust aggregation against outliers (e.g., PCR jackpots)
        summarise(across(where(is.numeric), \(x) median(x, na.rm = TRUE))) |>
        column_to_rownames("clone_id") |> t()

    # Remove specified clone from the profile (e.g., "small" or unassigned)
    if (exclude_clone %in% colnames(clonal_profile)) {
        clonal_profile <- clonal_profile[, colnames(clonal_profile) != exclude_clone, drop = FALSE]
    }

    ## 3. Normalize to the Diploid Reference Clone
    if (!(diploid_ref %in% colnames(clonal_profile))) {
        stop(paste("Diploid reference clone '", diploid_ref, "' not found in aggregated profiles."))
    }

    # Apply diploid normalization: 2 * Observed / Reference
    clonal_profile_diploid <- apply(clonal_profile, 2, function(x) {
        # Get the reference profile (amplicons x 1)
        ref <- clonal_profile[, diploid_ref]

        # Calculate the ratio and multiply by 2 (2N assumption for the reference)
        return(2 * x / ref)
    })

    return(clonal_profile_diploid)
}



#' Sort Chromosome Names in Genomic Order
#'
#' @description
#' This function takes a vector of chromosome names (e.g., "chr1", "1", "chrX")
#' and returns a sorted vector according to the standard genomic order (1, 2, ..., 22, X, Y).
#' Chromosomes not present in the input list are excluded from the result.
#'
#' @param chromosome_vector A character vector containing the names of the chromosomes.
#'
#' @return A character vector of unique, sorted chromosome names (with "chr" prefix).
#'
#' @examples
#' # Standard case
#' unsorted_list <- c("chr1", "chr10", "chr3", "chrX", "chr2")
#' sort_genomic_chromosomes(unsorted_list)
#' # [1] "chr1" "chr2" "chr3" "chr10" "chrX"
#'
#' # Handles the absence of chrX and chrY
#' no_sex_chroms <- c("chr11", "chr1", "chr2", "chr5")
#' sort_genomic_chromosomes(no_sex_chroms)
#' # [1] "chr1" "chr2" "chr5" "chr11"
#'
#' # Handles mixed formats and chrY
#' mixed_list <- c("1", "22", "chrY", "chr5", "chr13")
#' sort_genomic_chromosomes(mixed_list)
#' # [1] "chr1" "chr5" "chr13" "chr22" "chrY"
sort_genomic_chromosomes <- function(chromosome_vector) {

    # 1. Cleaning and Normalization: Remove "chr" prefix and normalize case for sex chromosomes.
    # We remove "chr" (case-insensitive) to extract the number or the letter.
    chrom_clean <- toupper(gsub("^(CHR|chr)", "", chromosome_vector, ignore.case = TRUE))

    # Work only with unique values
    chrom_unique <- unique(chrom_clean)

    # 2. Separate numerical and sex chromosomes
    chrom_nums <- chrom_unique[!chrom_unique %in% c("X", "Y")]
    chrom_sex <- chrom_unique[chrom_unique %in% c("X", "Y")]

    # 3. Sort numerical chromosomes correctly
    chrom_nums_sorted <- chrom_nums |>
        as.numeric() |>
        na.omit() |> # Removes elements that could not be converted to numbers (e.g., junk data)
        sort() |>
        as.character()

    # 4. Rebuild the final list in standard genomic order (1, 2, ..., 22, X, Y)
    final_order_list <- c(
        # Add the "chr" prefix back to the numbers
        paste0("chr", chrom_nums_sorted),

        # Add chrX and chrY only if they were present in the input
        if ("X" %in% chrom_sex) "chrX",
        if ("Y" %in% chrom_sex) "chrY"
    )

    # Return the sorted list
    return(final_order_list)
}


#' Rander annotation table
#'
render_annotation_table <- function(obj, ploidy_data){
    mat_data <- t(ploidy_data)

    # Extraction sécurisée des métadonnées
    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(obj$mae[["amplicons"]]))
    genome_v <- S4Vectors::metadata(obj$mae)$genome_version
    if (is.null(genome_v)) genome_v <- "hg19"

    tmp_var_table <- cnv_id_table |>
        dplyr::filter(dna_id %in% colnames(mat_data)) |>
        dplyr::arrange(as.numeric(chrom), as.numeric(start_pos)) |>
        dplyr::mutate(chr_lit = paste0("chr", chrom))

    mat_data <- mat_data[, tmp_var_table$dna_id, drop = FALSE]
    tmp_table <- tmp_var_table[match(colnames(mat_data), tmp_var_table$dna_id), ]

    annotate_genomic_regions(region_data = tmp_table, build = genome_v)
}

#' Plot CNV Heatmap
#'
#' @description
#' Generates a heatmap of Copy Number Variation (CNV) ploidy data.
#' The function aligns the provided ploidy data with chromosomal positions defined
#' in the ScIGMA object and visualizes it using `ComplexHeatmap`.
#'
#' @param obj A ScIGMA object containing a `cnv_id_table` with `dna_id`, `chrom`, and `start_pos` columns.
#' @param ploidy_data A matrix or data frame of ploidy values (result from `process_cnv_to_clonal_profile`).
#'
#' @return A `Heatmap` object (from the ComplexHeatmap package).
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar
#' @importFrom viridis viridis
#' @importFrom dplyr filter arrange mutate
#' @importFrom stats quantile setNames
#'
#' @export
plot_cnv_heatmap <- function(obj, ploidy_data, display_gene = FALSE) {

    mat_data <- t(ploidy_data)

    # Extraction sécurisée des métadonnées
    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(obj$mae[["amplicons"]]))
    genome_v <- S4Vectors::metadata(obj$mae)$genome_version
    if (is.null(genome_v)) genome_v <- "hg19"

    print("cnv_id_table")
    print(cnv_id_table)

    tmp_var_table <- cnv_id_table |>
        dplyr::filter(dna_id %in% colnames(mat_data)) |>
        dplyr::arrange(as.numeric(chrom), as.numeric(start_pos)) |>
        dplyr::mutate(chr_lit = paste0("chr", chrom))

    mat_data <- mat_data[, tmp_var_table$dna_id, drop = FALSE]
    tmp_split_table <- tmp_var_table[match(colnames(mat_data), tmp_var_table$dna_id), ]
    sorted_gen_levels <- sort_genomic_chromosomes(tmp_split_table$chrom)

    if (display_gene){
        tmp_split_vec <- annotate_genomic_regions(region_data = tmp_split_table, build = genome_v)
        split_vec <- factor(tmp_split_vec$symbol, levels = unique(tmp_split_vec$symbol))
    } else {
        split_vec <- factor(tmp_split_table$chr_lit, levels = sorted_gen_levels)
    }

    col_fun <- circlize::colorRamp2(
        breaks = c(quantile(mat_data, c(0, 0.25), na.rm=TRUE), 2, quantile(mat_data, c(0.75, 1), na.rm=TRUE)),
        colors = c("black", "#4575B4", "#F0F0F0", "#D73027", "#67001F")
    )

    group_colors <- setNames(viridis::viridis(nrow(mat_data)), nm = rownames(mat_data))

    left_ann <- ComplexHeatmap::rowAnnotation(
        df = data.frame(Group = rownames(mat_data)),
        col = list(Group = group_colors),
        show_legend = FALSE,
        simple_anno_size = grid::unit(1, "cm"),
        show_annotation_name = FALSE
    )

    ComplexHeatmap::Heatmap(
        mat_data,
        name = "Ploidy",
        col = col_fun,
        column_split = split_vec,
        row_split = rownames(mat_data),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        border = TRUE,
        row_gap = grid::unit(0, "mm"),
        column_gap = grid::unit(0, "mm"),
        border_gp = grid::gpar(col = "black", lwd = 1),
        column_title_side = "bottom",
        column_title_rot = 90,
        column_title_gp = grid::gpar(fontsize = 10),
        left_annotation = left_ann
    )
}


#' Annotate genomic regions using local Ensembl databases
#'
#' This function maps genomic coordinates to HGNC symbols using local
#' SQLite-backed Ensembl databases (EnsDb). It is specifically designed
#' for R Shiny applications to provide sub-second latency, zero network
#' dependency, and strict data privacy for clinical oncology data.
#'
#' @param region_data A data.frame containing genomic coordinates. It must
#'   contain the columns: 'chrom', 'start', and 'end'.
#' @param build A character string specifying the genome assembly.
#'   Must be either "hg38" (GRCh38) or "hg19" (GRCh37). Defaults to "hg38".
#'
#' @return A data.frame extending 'region_data' with 'hgnc_symbol' and
#'   'gene_biotype'. The output is filtered to retain only distinct
#'   protein-coding genes with valid HGNC symbols.
#'
#' @export
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom EnsDb.Hsapiens.v75 EnsDb.Hsapiens.v75
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom dplyr filter distinct
#' @importFrom ensembldb genes
annotate_genomic_regions <- function(region_data, build = "hg38") {
    # Selection logic for the local Ensembl database object
    # Ensures zero network latency for reactive Shiny environments
    if (build == "hg38") {
        edb <- EnsDb.Hsapiens.v86
    } else if (build == "hg19") {
        edb <- EnsDb.Hsapiens.v75
    } else {
        stop("Invalid build. Please use 'hg19' or 'hg38'.")
    }

    # Conversion to GRanges for high-performance spatial overlaps
    # This step leverages C-level optimization for genomic arithmetic
    query_gr <- GRanges(
        seqnames = region_data$chrom,
        ranges = IRanges(
            start = region_data$start_pos,
            end = region_data$end_pos
        )
    )

    # Extract gene-level metadata from the local SQLite backend
    # We restrict columns to minimize memory footprint in the Shiny session
    target_genes <- genes(
        edb,
        columns = c("symbol", "gene_biotype")
    )

    # C-optimized spatial join via findOverlaps
    overlaps <- findOverlaps(query_gr, target_genes)

    # Construct the result set using fast indexing
    # Prioritizing protein-coding genes for biological relevance in oncology
    tmp_y <- data.frame(symbol = target_genes$symbol[subjectHits(overlaps)],
                        gene_biotype = target_genes$gene_biotype[subjectHits(overlaps)],
                        dna_id = region_data[queryHits(overlaps), "dna_id"]) |>
        filter(gene_biotype == "protein_coding" & !duplicated(dna_id))

    result_df <- region_data |> left_join(y=tmp_y, by = "dna_id")

    # Replace NAs
    result_df[is.na(result_df)] <- "Unknown"

    # result_df <- data.frame(
    #     region_data[queryHits(overlaps), ],
    #     hgnc_symbol = target_genes$symbol[subjectHits(overlaps)],
    #     gene_biotype = target_genes$gene_biotype[subjectHits(overlaps)]
    # ) |>
    #     filter(gene_biotype == "protein_coding") |>
    #     distinct()


    return(result_df)
}


#' Plot Genomic CNV Profile (Ploidy) - Bold & Larger Labels
#' @export
plot_cnv_genome <- function(cnv_matrix,
                            gene_annotation = NULL,
                            sub_indices = NULL,
                            title = "Genomic Ploidy Profile",
                            max_points = 500,
                            lineplot_type = "Genes+amplicons") {

    # --- 1. Data Preparation ---
    if (!is.null(sub_indices)) {
        mat_sub <- cnv_matrix[sub_indices, , drop = FALSE]
    } else {
        mat_sub <- cnv_matrix
    }

    probes <- colnames(mat_sub)
    probe_means <- colMeans(mat_sub, na.rm = TRUE)
    # Base metadata
    plot_meta <- data.frame(
        Probe = probes,
        Mean_Ploidy = probe_means,
        stringsAsFactors = FALSE
    )
    # --- 2. Handling Order & X-Axis Logic ---

    # Valeurs par défaut
    plot_meta$x_index <- 1:nrow(plot_meta)
    x_tick_text <- paste0("<b>", probes, "</b>") # GRAS par défaut
    x_tick_vals <- plot_meta$x_index
    separators <- list()
    angle_val <- -90

    if (!is.null(gene_annotation)) {

        # A. Jointure
        plot_meta <- suppressMessages(
            dplyr::left_join(plot_meta, gene_annotation, by = "Probe")
        )

        # B. Tri Génomique
        plot_meta <- plot_meta %>%
            dplyr::mutate(
                chr_sort = suppressWarnings(as.numeric(Chromosome)),
                chr_sort = ifelse(Chromosome == "X", 98, chr_sort),
                chr_sort = ifelse(Chromosome == "Y", 99, chr_sort)
            ) %>%
            dplyr::arrange(chr_sort, Chrom_start)

        # C. Branchement Conditionnel
        if (lineplot_type == "Genes+amplicons") {
            # MODE STACKED (Gènes)
            angle_val <- -45

            plot_meta$group_label <- ifelse(is.na(plot_meta$Gene), plot_meta$Probe, plot_meta$Gene)
            unique_groups <- unique(plot_meta$group_label)

            plot_meta$x_index <- as.numeric(factor(plot_meta$group_label, levels = unique_groups))
            x_tick_vals <- 1:length(unique_groups)

            # Application du GRAS sur les labels des gènes
            x_tick_text <- paste0("<b>", unique_groups, "</b>")

            sep_pos <- seq(1.5, length(unique_groups) - 0.5, 1)

            # Recalcul des moyennes par gène
            tp_mean <- plot_meta |> dplyr::group_by(Gene) |> dplyr::summarize(Mean_Ploidy = mean(Mean_Ploidy, na.rm=TRUE))
            tp_mean_vec <- setNames(tp_mean$Mean_Ploidy, nm = tp_mean$Gene)
            plot_meta$Mean_Ploidy <- tp_mean_vec[plot_meta$Gene]

        } else {
            # MODE SPREAD (Positions)
            angle_val <- -45

            plot_meta$x_index <- 1:nrow(plot_meta)
            plot_meta$group_label <- plot_meta$Chrom_pos

            rle_res <- rle(as.character(plot_meta$group_label))
            end_indices <- cumsum(rle_res$lengths)
            sep_pos <- head(end_indices, -1) + 0.5
            start_indices <- c(1, head(end_indices, -1) + 1)
            x_tick_vals <- (start_indices + end_indices) / 2

            # Application du GRAS sur les labels des chromosomes
            x_tick_text <- paste0("<b>", rle_res$values, "</b>")
        }

        # D. Ré-ordonner
        ordered_probes <- plot_meta$Probe
        ordered_probes <- ordered_probes[ordered_probes %in% colnames(mat_sub)]
        mat_sub <- mat_sub[, ordered_probes, drop = FALSE]
        plot_meta <- plot_meta[match(ordered_probes, plot_meta$Probe), ]

        # E. Séparateurs
        separators <- lapply(sep_pos, function(p) {
            list(
                type = "line",
                x0 = p, x1 = p,
                y0 = 0, y1 = 1, yref = "paper",
                line = list(color = "#cccccc", width = 1, dash = "dot")
            )
        })
    }

    # --- 3. Prepare Scatter Data ---

    n_cells <- nrow(mat_sub)
    if (n_cells > max_points) {
        set.seed(42)
        mat_viz <- mat_sub[sample(n_cells, max_points), , drop = FALSE]
    } else {
        mat_viz <- mat_sub
    }

    probe_to_x <- setNames(plot_meta$x_index, plot_meta$Probe)
    x_indices_vector <- probe_to_x[colnames(mat_viz)]

    df_scatter <- data.frame(
        Probe = rep(colnames(mat_viz), each = nrow(mat_viz)),
        x_idx = rep(x_indices_vector, each = nrow(mat_viz)),
        Ploidy = as.vector(mat_viz)
    )

    # --- 4. Plotting ---

    p <- plotly::plot_ly() %>%
        plotly::layout(
            # Titre plus gros (20)
            title = list(text = paste0("<b>", title, "</b>"), font = list(size = 20)),
            template = "plotly_white",
            font = list(family = "Arial, sans-serif", size = 12),
            xaxis = list(
                title = "",
                tickmode = "array",
                tickvals = x_tick_vals,
                ticktext = x_tick_text, # Contient les tags <b>
                tickangle = angle_val,
                ticks = "outside",
                tickfont = list(size = 13), # Taille augmentée pour l'axe X (Gènes)
                automargin = TRUE,
                zeroline = FALSE,
                showgrid = FALSE
            ),
            yaxis = list(
                # Titre axe Y : GRAS + TAILLE 18
                title = list(text = "<b>Ploidy</b>", font = list(size = 18)),
                range = c(0, 5),
                zeroline = TRUE,
                zerolinewidth = 1,
                gridcolor = "#eeeeee",
                tickfont = list(size = 14) # Taille augmentée pour les chiffres Y
            ),
            legend = list(
                orientation = "h",
                xanchor = "center", x = 0.5, y = -0.2
            ),
            shapes = c(
                list(
                    list(
                        type = "rect",
                        fillcolor = "rgba(200, 200, 200, 0.25)",
                        line = list(width = 0),
                        xref = "paper", x0 = 0, x1 = 1,
                        y0 = 1.5, y1 = 2.5
                    )
                ),
                separators
            ),
            margin = list(b = 100)
        )

    # Add Cells (Blue)
    p <- p %>% plotly::add_trace(
        data = df_scatter,
        x = ~x_idx,
        y = ~Ploidy,
        type = "scatter",
        mode = "markers",
        marker = list(
            color = "#4a86e8",
            size = 6,
            opacity = 0.7,
            line = list(width = 0.5, color = "white")
        ),
        name = "Cells",
        hoverinfo = "text",
        text = ~paste("<b>Probe:</b>", Probe, "<br><b>Ploidy:</b>", round(Ploidy, 2))
    )

    # Add Means (Red)
    p <- p %>% plotly::add_trace(
        x = plot_meta$x_index,
        y = plot_meta$Mean_Ploidy,
        type = "scatter",
        mode = "markers",
        marker = list(
            color = "#cc0000",
            size = 10,
            symbol = "circle",
            line = list(width = 1, color = "black")
        ),
        name = "Mean",
        hoverinfo = "text",
        text = ~paste("<b>Probe:</b>", plot_meta$Probe,
                      "<br><b>Gene:</b>", if("Gene" %in% names(plot_meta)) plot_meta$Gene else "N/A",
                      "<br><b>Mean:</b>", round(plot_meta$Mean_Ploidy, 2))
    ) %>%
        plotly::toWebGL() %>%
        plotly::config(
            displaylogo = FALSE,
            toImageButtonOptions = list(
                format = "svg",
                filename = "genomic_ploidy_profile",
                height = 600,
                width = 1000,
                scale = 2
            )
        )

    return(p)
}
