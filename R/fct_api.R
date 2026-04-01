#' Normalize variant separators
#'
#' Replace `:` and `/` with `-` in variant strings. Useful to normalize
#' variant identifiers for URLs, file names, or APIs that do not accept
#' colon or slash characters.
#'
#' @param variant A character (or factor) vector of variant strings
#'   (e.g., `"chr1:115256598-T/C"`).
#'
#' @return A character vector with `:` and `/` replaced by `-`.
#'
#' @details
#' - Vectorized over `variant`.
#' - `NA` inputs remain `NA`.
#' - This function does not trim whitespace or collapse repeated hyphens.
#'
#' @examples
#' format_variant("chr1:115256598-T/C")
#' # "chr1-115256598-T-C"
#'
#' format_variant(c("chr2:200/A", NA))
#' # "chr2-200-A" NA
format_variant <- function(variant) {
    stopifnot(is.character(variant) || is.factor(variant))
    variant <- as.character(variant)
    chartr(":/", "--", variant)
}


#' Normalize a nested path specification
#'
#' Accepts either a single string in dot-notation (e.g., `"a.b.c"`)
#' or a character vector (e.g., `c("a", "b", "c")`) and returns
#' the path as a character vector of segments.
#'
#' @param path Character scalar in dot-notation or a character vector
#'   of path segments.
#'
#' @return A character vector of path segments.
#'
#' @examples
#' make_path("a.b.c")
#' make_path(c("a", "b", "c"))
make_path <- function(path) {
    if (is.character(path)){
        strsplit(path, "\\.", fixed = FALSE)[[1]]
    } else {
        path
    }
}

#' Safely pluck a nested value, returning NA for empty/missing inputs
#'
#' Uses \code{purrr::pluck()} to extract a nested element from a record.
#' If the value is \code{NULL}, an empty list, or a length-0 vector,
#' returns \code{NA}. Otherwise returns the extracted value as-is.
#'
#' @param x A list-like record (e.g., one parsed JSON object).
#' @param path Path to the element: dot-notation string like \code{"a.b.c"}
#'   or a character vector like \code{c("a","b","c")}.
#'
#' @return The extracted value, or \code{NA} if missing/empty.
#'
#' @details
#' Empty structures are normalized to \code{NA} so downstream binding
#' and tabular transformations are stable.
#'
#' @examples
#' rec <- list(a = list(b = list(c = 1)))
#' pluck_or_na(rec, "a.b.c")  # 1
#' pluck_or_na(rec, "a.b.d")  # NA
pluck_or_na <- function(x, path) {
    v <- purrr::pluck(x, !!!make_path(path), .default = NULL)
    # NA if NULL / empty list / length-0 vector
    if (is.null(v) || (is.list(v) && length(v) == 0) || (length(v) == 0)) {
        NA
    } else {
        v
    }
}

#' Extract multiple nested fields from a list of records
#'
#' Given a list of records (e.g., parsed JSON objects) and a named list of
#' paths, extracts each requested path into a tibble column. Each extracted
#' value is wrapped to produce a list-column, preserving heterogeneous
#' structures while enabling safe row-binding.
#'
#' @param records A list of records (each element is a list parsed from JSON).
#' @param paths A named list mapping output column names to paths.
#'   Paths can be dot-notation strings (e.g., \code{"a.b.c"}) or character
#'   vectors of segments (e.g., \code{c("a","b","c")}).
#'
#' @return A tibble with one row per input record and one column per requested
#'   path. Columns are list-columns to preserve nested/heterogeneous values.
#'
#' @examples
#' records <- list(
#'   list(id = 1, annotations = list(clinvar = list(value = list(sig = "pathogenic")))),
#'   list(id = 2, annotations = list())  # missing value
#' )
#' paths <- list(
#'   id = "id",
#'   clinvar_value = "annotations.clinvar.value"
#' )
#' extract_paths(records, paths)
#'
#' # Then, if desired, widen:
#' # extract_paths(records, paths) |> tidyr::unnest_wider(clinvar_value, names_sep = "_")
extract_paths <- function(records, paths) {
    # Example of 'paths': list(col_name = "a.b.c", other = "x.y")
    stopifnot(is.list(records))
    purrr::map(records, function(rec) {
        # Ensure list-columns by wrapping each plucked value in a list()
        vals <- purrr::imap(paths, ~ list(pluck_or_na(rec, .x)))  # guarantees list-columns
        tibble::as_tibble(vals)
    }) |>
        dplyr::bind_rows()
}

#' @import httr2
#' @import dplyr
#'
fetch_variants_batch_fields <- function(
        ids,
        paths,                                   # list(col_name = "a.b.c", ...)
        base_url = "https://api.missionbio.io/annotations/v1/variants",
        batch_size = 300,
        max_retries = 4
) {
    stopifnot(length(ids) >= 1L)
    ids_formated <- format_variant(ids)

    chunks <- split(ids_formated, ceiling(seq_along(ids_formated) / batch_size))

    message(paste0("Annotating ", length(chunks), " batches ..."))

    purrr::map_dfr(chunks, function(ch) {

        query_url <- paste0(base_url, "?ids=", paste0(ch, collapse = ","))

        res_query <- query_url |>
            httr2::request() |>
            httr2::req_retry(max_tries = max_retries, backoff = ~ min(60, 2^(.x - 1))) |> # Backoff & retry
            httr2::req_perform() |>
            httr2::resp_body_json()


        # colonne id + toutes les extractions demandées
        tibble::tibble(variant_id = ch) |>
            dplyr::bind_cols(extract_paths(res_query, paths)) |>
            dplyr::mutate("gene_function" = sapply(gene_function, paste0,collapse = ", "),
                   "clinvar" = sapply(clinvar, paste0,collapse = ", ")) |>
            dplyr::mutate(
                chromosome   = purrr::map_chr(chromosome, 1, .default = NA_character_),
                position     = purrr::map_chr(position, 1, .default = NA_character_),
                ref_allele   = purrr::map_chr(ref_allele, 1, .default = NA_character_),
                alt_allele   = purrr::map_chr(alt_allele, 1, .default = NA_character_),
                gene         = purrr::map_chr(gene, 1, .default = NA_character_),
                variant_type = purrr::map_chr(variant_type, 1, .default = NA_character_),
                impact = purrr::map_vec(impact,as.numeric,.default = NA_real_),
                clinvar = purrr::map_chr(clinvar, as.character ,.default = NA_character_),
                clinvar = data.table::fifelse(clinvar == "NA", NA, clinvar),
                variant_id = paste0(gene,":",sub("^([^-]+)-([^-]+)-([^-]+)-([^-]+)$", "\\1:\\2:\\3/\\4",variant_id))
            )
    })
}

#' Filter variants and fetch annotations
#'
#' @param obj ScIGMA_object instance.
#' @param paths List of API paths for fetch_variants_batch_fields.
#' @param min_dp Integer. Minimum depth.
#' @param min_gq Integer. Minimum genotype quality.
#' @param vaf_ref Numeric. Maximum VAF for reference call.
#' @param vaf_hom Numeric. Minimum VAF for homozygous alternate call.
#' @param vaf_het Numeric. Maximum VAF for heterozygous call.
#' @param min_cell_pt Numeric. Minimum percentage of cells covering a variant.
#' @param min_mut_cell_pt Numeric. Minimum percentage of mutated cells.
#' @param batch_size Integer. Batch size for API calls.
#' @return Filtered and annotated ScIGMA_object.
#' @export
filter_and_annotate_variants <- function(obj,
                                         paths,
                                         min_dp = 10,
                                         min_gq = 30,
                                         vaf_ref = 5,
                                         vaf_hom = 95,
                                         vaf_het = 30,
                                         min_cell_pt = 10,
                                         min_mut_cell_pt = 10,
                                         batch_size = 300) {

    # 1. Pipeline execution: Filtering
    obj <- filter_variant_ScIGMA_mae(
        obj = obj,
        min.dp = min_dp,
        min.gq = min_gq,
        vaf.ref = vaf_ref,
        vaf.hom = vaf_hom,
        vaf.het = vaf_het,
        min.cell.pt = min_cell_pt,
        min.mut.cell.pt = min_mut_cell_pt
    )

    message("Fetching and injecting annotations...")

    # 2. Strict ID extraction
    variant_ids <- rownames(obj$mae[["dna_variants"]])

    # Gatekeeper: Do not call API if filtering removed everything
    if (length(variant_ids) == 0) {
        warning("No variants left after filtering. Skipping annotation.")
        return(obj)
    }

    # 3. Protected API call
    annot_df <- tryCatch({
        fetch_variants_batch_fields(
            variant_ids,
            batch_size = batch_size,
            paths = paths
        )
    }, error = function(e) {
        stop(sprintf("API Error during variant annotation: %s", e$message))
    })

    if (!is.null(annot_df) && nrow(annot_df) > 0) {

        # Robust string cleaning (prevents truncating legitimate IDs)
        annot_df$query <- ifelse(
            grepl("^[A-Za-z0-9]+:chr", annot_df$variant_id),
            sub("^[^:]*:", "", annot_df$variant_id),
            annot_df$variant_id
        )



        # 4. Dimensional synchronization
        current_rowdata <- as.data.frame(
            SummarizedExperiment::rowData(obj$mae[["dna_variants"]])
        )
        current_rowdata$query_id <- rownames(current_rowdata)

        merged_rowdata <- merge(
            x = current_rowdata,
            y = annot_df,
            by.x = "query_id",
            by.y = "query",
            all.x = TRUE,
            sort = FALSE
        )

        # 5. MAE Topology Guardrail (CRITICAL)
        # Prevents 1-to-many API mappings from expanding the rowData and crashing the MAE
        if (nrow(merged_rowdata) > nrow(current_rowdata)) {
            warning("API returned multiple annotations per variant. Deduplicating to preserve MAE strict dimensions.")
            merged_rowdata <- merged_rowdata[!duplicated(merged_rowdata$query_id), ]
        }

        # 6. Canonical order restoration
        rownames(merged_rowdata) <- merged_rowdata$query_id
        merged_rowdata <- merged_rowdata[variant_ids, ]
        merged_rowdata$query_id <- NULL

        merged_rowdata$variant_id[is.na(merged_rowdata$variant_id)] <- paste0("Unmapped:",rownames(merged_rowdata))[is.na(merged_rowdata$variant_id)]

        print("merged_rowdata")
        print(merged_rowdata)

        # 7. In-place MAE injection
        SummarizedExperiment::rowData(obj$mae[["dna_variants"]]) <- S4Vectors::DataFrame(merged_rowdata)
        message("Annotation matrix successfully integrated into MAE rowData.")

    } else {
        warning("API returned an empty object. rowData remains unannotated.")
    }

    message("Calculating variant cell proportions...")

    # 1. Extraction du pointeur de la matrice (Aucune donnée chargée en RAM)
    vaf_mtx <- SummarizedExperiment::assay(obj$mae[["dna_variants"]], "vaf")

    # 2. Calcul vectorisé Out-Of-Core sur les lignes (Variants)
    # L'opérateur > 0 exclut nativement les Wild-Type (0) et les génotypes non-fiables (-1)
    mutated_cells_count <- DelayedMatrixStats::rowSums2(vaf_mtx > 0)

    # Le registre absolu des cellules dicte le dénominateur
    total_cells <- ncol(vaf_mtx)

    # 3. Injection directe et in-place dans le registre d'annotation
    SummarizedExperiment::rowData(obj$mae[["dna_variants"]])$cell_proportion <- mutated_cells_count / total_cells

    message("Cell proportions successfully added to rowData.")
    S4Vectors::metadata(obj$mae)$variant_filter <- "filtered"

    invisible(obj)
}
