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
    stopifnot(is.list(records), length(paths) >= 1L)

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
        obj,
        paths,                                   # list(col_name = "a.b.c", ...)
        base_url = "https://api.missionbio.io/annotations/v1/variants",
        batch_size = 300,
        max_retries = 4
) {
    ids <- format_variant(obj$data$variants)

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
                position     = purrr::map_chr(position, 1, .default = NA_character_),   # si numérique : map_dbl
                ref_allele   = purrr::map_chr(ref_allele, 1, .default = NA_character_),
                alt_allele   = purrr::map_chr(alt_allele, 1, .default = NA_character_),
                gene         = purrr::map_chr(gene, 1, .default = NA_character_),
                variant_type = purrr::map_chr(variant_type, 1, .default = NA_character_),
                impact = purrr::map_vec(impact,as.numeric,.default = NA_real_),
                clinvar = purrr::map_chr(clinvar, as.character ,.default = NA_character_),
                clinvar = data.table::fifelse(clinvar == "NA", NA, clinvar)
            )
    })
}
