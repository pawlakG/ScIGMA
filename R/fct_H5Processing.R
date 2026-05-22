
# ScIGMA — HDF5‑backed (DelayedArray/HDF5Array)
# ----------------------------------------------------
# Dependencies:
#   - BiocManager::install(c("DelayedArray", "HDF5Array", "BiocParallel", "rhdf5"))

suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
    library(BiocParallel)
    library(rhdf5)
    library(R6)
})

# ----------------------------------------------------
# HDF5 helpers
# ----------------------------------------------------
#' Safely read a 1-D character dataset from an HDF5 file
#'
#' This helper function attempts to read a one-dimensional dataset
#' (strings or integers) from an HDF5 file. The result is always
#' returned as a character vector. If the specified dataset does not
#' exist or an error occurs, a default value is returned instead.
#'
#' @param filepath Character string. Path to the HDF5 file.
#' @param path Character string. Internal path to the dataset within the HDF5 file.
#' @param default Character vector. Value to return if the dataset is missing
#'   or cannot be read. Defaults to an empty character vector.
#'
#' @return A character vector with the dataset values, or \code{default}
#'   if the dataset is not found or cannot be read.
#'
#' @examples
#' \dontrun{
#' # Suppose "sample_names" is stored under /metadata in file.h5
#' .h5read_char("file.h5", "/metadata/sample_names")
#'
#' # Return a default value if the dataset does not exist
#' .h5read_char("file.h5", "/invalid/path", default = "unknown")
#' }
#'
#' @importFrom rhdf5 h5read
#' @keywords internal
.h5read_char <- function(filepath, path, default = character()) {
    out <- tryCatch({
        v <- rhdf5::h5read(filepath, path)
        if (is.null(v)) default else as.character(v)
    }, error = function(e) default)
    out
}


#' Check if a path exists in an HDF5 file
#'
#' This helper function tests whether a given path (dataset or group)
#' exists within an HDF5 file. It safely opens the file, checks for
#' the path, and closes the file afterwards.
#'
#' @param filepath Character string. Path to the HDF5 file.
#' @param path Character string. Path to the dataset or group inside the HDF5 file.
#'
#' @return Logical value: \code{TRUE} if the path exists, \code{FALSE} otherwise.
#'
#' @examples
#' \dontrun{
#' # Check if a dataset exists inside the file
#' .h5_has_path("file.h5", "/assays/dna_variants/layers/AF")
#' }
#'
#' @importFrom rhdf5 H5Fopen H5Fclose H5Lexists
#' @keywords internal
.h5_has_path <- function(filepath, path) {
    fid <- rhdf5::H5Fopen(filepath)
    on.exit(rhdf5::H5Fclose(fid), add = TRUE)
    rhdf5::H5Lexists(fid, path)
}


#' Create a DelayedArray reference to an HDF5 dataset
#'
#' This helper function constructs an \code{HDF5Array} (a DelayedArray
#' backed by an HDF5 dataset) pointing to a specified dataset inside an HDF5 file.
#' If the dataset does not exist, an explicit error is raised.
#'
#' @param filepath Character string. Path to the HDF5 file.
#' @param path Character string. Path to the dataset within the HDF5 file.
#'
#' @return An \code{HDF5Array} object representing the dataset. This object is
#'   a DelayedArray and enables block-processed, out-of-memory operations
#'   without loading the full dataset into RAM.
#'
#' @examples
#' \dontrun{
#' # Access a dataset of variant allele frequencies
#' vaf <- .h5_delayed("file.h5", "/assays/dna_variants/layers/AF")
#' dim(vaf)
#' }
#'
#' @importFrom HDF5Array HDF5Array
#' @keywords internal
.h5_delayed <- function(filepath, path) {
    if (!.h5_has_path(filepath, path)) {
        stop(sprintf("Dataset not found: %s", path))
    }
    tryCatch({
        HDF5Array(filepath, path)
    }, error = function(e) {
        message("Error accessing dataset: ", path, "\n", conditionMessage(e), "\n",
                "Trying another method ,\nThis may be longer ...")
        tryCatch({
            h5 <- rhdf5::h5read(filepath, path)
            # Matrix transformed to a vector, we need to restore dimensions
            tmp <- as.logical(h5)
            dim(tmp) <- dim(h5)
            DelayedArray(tmp)
        }, error = function(e) {
            stop("Failed to read dataset: ", path, "\n", conditionMessage(e))
        })
    })
}


#' Read all datasets directly under an HDF5 group
#'
#' Lists the contents of a given group in an HDF5 file and returns a
#' named list of datasets found directly under that group (no recursion).
#' Each dataset is read with `rhdf5::h5read()`. Length-1 datasets are
#' coerced to a simple vector.
#'
#' @param file Character string. Path to the HDF5 file.
#' @param group_path Character string. Absolute path to the group inside the HDF5 file.
#' @return A named list where each element corresponds to a dataset under `group_path`.
#' @keywords internal
.h5_read_metadata_group <- function(file, group_path) {
    # List all entries in the HDF5 file
    all <- h5ls(file, recursive = TRUE, datasetinfo = FALSE)

    # Keep only direct children of the target group
    ls <- subset(all, group == group_path)
    if (nrow(ls) == 0L) stop("Group not found or empty: ", group_path)

    # Keep datasets only
    ds <- subset(ls, otype == "H5I_DATASET")
    if (nrow(ds) == 0L) stop("No datasets under: ", group_path)

    # Read each dataset and return a named list
    out <- setNames(
        lapply(ds$name, function(nm) {
            val <- h5read(file, paste0(group_path, "/", nm))
            if (length(val) == 1L) val <- as.vector(val)
            val
        }),
        ds$name
    )
    out
}




#' ScIGMA_object: HDF5-backed R6 class for single-cell multi-omics
#'
#' @description
#' `ScIGMA_object` is an R6 class designed to manage large single-cell
#' genomics/proteomics matrices using the **DelayedArray/HDF5Array** backend.
#' All heavy matrices are expected to be HDF5-backed `DelayedArray` objects,
#' allowing out-of-memory processing and block-wise computation without
#' materializing the full data in RAM.
#'
#' @details
#' The class stores cell- and feature-level metadata along with several
#' assay matrices (VAF, genotype, depth, genotype quality, amplicon counts,
#' proteins, and optional ploidy). A convenience method `realize_all()` can
#' persist any delayed computations to a single HDF5 file, with control over
#' chunking and compression.
#'
#' **Important:** matrix-like fields (those ending with `*.mtx`) must be
#' `DelayedArray`-compatible objects (typically `HDF5Array`). The class does
#' not implicitly coerce to in-memory matrices.
#'
#' @section Fields:
#' - `meta.data` (`character`): Sample-level metadata.
#' - `cell.ids` (`character`): Cell identifiers (barcodes).
#' - `cell.labels` (`character`): Human-readable cell labels.
#' - `variants` (`character`): Variant identifiers.
#' - `variant.filter` (`character`): Subset of `variants` after filtering.
#' - `vaf.mtx` (`DelayedArray`): Variant Allele Frequency matrix.
#' - `vaf.mtx.cells` (`character`): Column/row labels associated with `vaf.mtx`.
#' - `gt.mtx` (`DelayedArray`): Genotype matrix.
#' - `gt.mtx.cells` (`character`): Labels for `gt.mtx`.
#' - `dp.mtx` (`DelayedArray`): Read depth matrix.
#' - `dp.mtx.cells` (`character`): Labels for `dp.mtx`.
#' - `gq.mtx` (`DelayedArray`): Genotype quality matrix.
#' - `gq.mtx.cells` (`character`): Labels for `gq.mtx`.
#' - `amps` (`character`): Amplicon identifiers.
#' - `amp.normalize.method` (`character`): Normalization method used for amplicons.
#' - `amp.mtx` (`DelayedArray`): Amplicon read count matrix.
#' - `amp.mtx.cells` (`character`): Labels for `amp.mtx`.
#' - `ploidy.mtx` (`DelayedArray` or `NULL`): Optional ploidy matrix.
#' - `proteins` (`character`): Protein identifiers.
#' - `protein.normalize.method` (`character`): Normalization method for proteins.
#' - `protein.mtx` (`DelayedArray`): Protein expression/read count matrix.
#' - `protein.mtx.cells` (`character`): Labels for `protein.mtx`.
#' - `backing_files` (`list`): File bookkeeping (e.g., `original`, `realized`).
#'
#' @section Methods:
#' - `initialize(...)`: Constructor; validates HDF5-backed matrix fields.
#' - `print()`: Compact console summary of the object.
#' - `summary()`: Returns a list with counts and key matrix dimensions.
#' - `add_variant(variant)`: Add a new variant ID (de-duplicated).
#' - `filter_variants(filter_fun)`: Keep variants for which `filter_fun(id)` is TRUE.
#' - `realize_all(dir, file, chunkdim, level)`: Persist delayed matrices to a
#'    single HDF5 file with optional chunking/compression control.
#' )


# ----------------------------------------------------
# 100% HDF5‑backed import function
# ----------------------------------------------------
# Reads a .h5 file and constructs a ScIGMA_object whose matrices
# are HDF5Array/DelayedArray (no full in‑RAM copies).

#' Load an H5 file into an HDF5‑backed ScIGMA_object
#' @param filepath Path to the `.h5` file.
#' @param sample.name Sample name (stored in `meta.data`).
#' @param omic.type Either "DNA+protein" or "DNA".
#' @param block.size Target DelayedArray block size in bytes (e.g., `100e6`).
#' @return A `ScIGMA_object` with all matrices HDF5‑backed.
loadH5_HDF5 <- function(filepath, sample.name, omic.type = c("DNA+protein", "DNA"), block.size = 100e6){
    omic.type <- match.arg(omic.type)
    stopifnot(file.exists(filepath))

    options(DelayedArray.block.size = block.size)

    # ---- Protein (optional) ----
    if (omic.type == "DNA+protein"){
        proteins <- .h5read_char(filepath, "/assays/protein_read_counts/ca/id")
        protein_mtx <- t(.h5_delayed(filepath, "/assays/protein_read_counts/layers/read_counts"))

        if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
            protein_samples <- "/assays/protein_read_counts/ra/label"
        } else {
            protein_samples <- "/assays/protein_read_counts/ra/sample_name"
        }
        protein_id_path <- "/assays/protein_read_counts/ra/barcode"
        protein_barcodes <- .h5read_char(filepath, protein_id_path)
        pr_ra_names <- .h5read_char(filepath, protein_samples)

        if (length(protein_barcodes) > 0) rownames(protein_mtx) <- protein_barcodes
        if (length(proteins) > 0) colnames(protein_mtx) <- proteins

        protein.normalize.method <- "unnormalized"
        protein.cells <- setNames(pr_ra_names, nm = protein_barcodes)
    } else {
        message("DNA only: skip protein matrix")
        proteins <- "non-protein"
        protein.normalize.method <- "non-protein"
        protein_mtx <- NULL
        protein.cells <- character()
    }

    # ---- Check if label information exists ----
    ## dna_read_counts
    if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
        dna_read_count_samples <- "/assays/dna_read_counts/ra/label"
    } else {
        dna_read_count_samples <- "/assays/dna_read_counts/ra/sample_name"
    }

    amp_sample_ids <- .h5read_char(filepath, "/assays/dna_read_counts/ra/barcode") # Sample barcode
    amp_ca_ids <- .h5read_char(filepath, "/assays/dna_read_counts/ca/id") # Amplicon ids
    amp_sample_names <- setNames(amp_sample_ids, nm = .h5read_char(filepath, dna_read_count_samples))# Samples human readable names

    ## dna_variants
    if (.h5_has_path(filepath, "/assays/dna_variants/ra/label")){
        dna_variants_samples <- "/assays/dna_variants/ra/label"
    } else {
        dna_variants_samples <- "/assays/dna_variants/ra/sample_name"
    }
    dna_sample_name <- .h5read_char(filepath, dna_variants_samples)
    dna_sample_ids <- .h5read_char(filepath, "/assays/dna_variants/ra/barcode")
    dna_id <- .h5read_char(filepath, "/assays/dna_variants/ca/id")

    vaf_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/AF"))
    if (length(dna_sample_ids) > 0 || length(dna_id) > 0) dimnames(vaf_mtx) <- list(dna_sample_ids, dna_id)

    gt_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/NGT"))
    if (length(dna_sample_ids) > 0 || length(dna_id) > 0) dimnames(gt_mtx) <- list(dna_sample_ids, dna_id)

    dp_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/DP"))
    if (length(dna_sample_ids) > 0 || length(dna_id) > 0) dimnames(dp_mtx) <- list(dna_sample_ids, dna_id)

    gq_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/GQ"))
    if (length(dna_sample_ids) > 0 || length(dna_id) > 0) dimnames(gq_mtx) <- list(dna_sample_ids, dna_id)

    dna_id_table <- data.frame(dna_id  = dna_id,
                               amplicon = .h5read_char(filepath, "/assays/dna_variants/ca/amplicon"),
                               chrom = .h5read_char(filepath, "/assays/dna_variants/ ca/CHROM"),
                               pos = .h5read_char(filepath, "/assays/dna_variants/ca/POS"),
                               ado_gt_cells = .h5read_char(filepath, "/assays/dna_variants/ca/ado_gt_cells"),
                               ado_rate = .h5read_char(filepath, "/assays/dna_variants/ca/ado_rate"),
                               filtered = .h5read_char(filepath, "/assays/dna_variants/ca/filtered"))

    cnv_id_table <- data.frame(dna_id  = .h5read_char(filepath, "/assays/dna_read_counts/ca/id"),
                               chrom = .h5read_char(filepath, "/assays/dna_read_counts/ca/CHROM"),
                               start_pos = .h5read_char(filepath, "/assays/dna_read_counts/ca/start_pos") |>
                                   as.numeric(),
                               end_pos = .h5read_char(filepath, "/assays/dna_read_counts/ca/end_pos") |>
                                   as.numeric())

    dna_cell_table <- data.frame(dna_barcode = .h5read_char(filepath, "/assays/dna_variants/ra/barcode"),
                                 dna_filtered = .h5read_char(filepath, "/assays/dna_variants/ra/filtered"),
                                 dna_sample_name = .h5read_char(filepath, dna_read_count_samples))


    # ---- DNA read counts (amplicons) ----
    amp_mtx <- t(.h5_delayed(filepath, "/assays/dna_read_counts/layers/read_counts"))
    if (length(amp_sample_ids) > 0 || length(amp_ca_ids) > 0) dimnames(amp_mtx) <- list(amp_sample_ids, amp_ca_ids)

    amp.normalize.method <- "unnormalized"

    # ---- Metadata ----
    dna_variant_metadata <- .h5_read_metadata_group(filepath, "/assays/dna_variants/metadata")
    dna_read_counts_metadata <- .h5_read_metadata_group(filepath, "/assays/dna_read_counts/metadata")
    protein_read_counts_metadata <- .h5_read_metadata_group(filepath, "/assays/protein_read_counts/metadata")

    # ---- Cell labels (logique de fallback identique au code d'origine) ----
    if (.h5_has_path(filepath, dna_read_count_samples)){
        cell.labels <-  amp_sample_names
    } else if (.h5_has_path(dna_variants_samples)){
        cell.labels <- setNames(dna_sample_ids, nm = dna_sample_name)
    } else if (.h5_has_path(filepath, protein_samples)){
        cell.labels <- protein.cells
    } else if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/barcode")){
        barcode <- .h5read_char(filepath, "/assays/dna_read_counts/ra/barcode")
        cell.labels <- setNames(barcode, nm = rep("unassigned", length(barcode)))
    } else {
        cell.labels <- character()
    }

    # ---- Variant filter mask ----
    cnv_metadata <- list("genome_version" = .h5read_char(filepath, "/assays/dna_read_counts/metadata/genome_version")
    )

    # ---- Variant filter mask ----
    variant.filter.mask <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/FILTER_MASK"))
    if (length(dna_sample_name) > 0 || length(dna_id) > 0) dimnames(variant.filter.mask) <- list(dna_sample_ids, dna_id)

    # ---- Object assembly ----
    obj <- ScIGMA_object$new(
        meta.data = list("name" = sample.name,
                         "dna_variant" = dna_variant_metadata,
                         "dna_read_counts" = dna_read_counts_metadata,
                         "protein_read_counts" = protein_read_counts_metadata
        ),
        cell.ids = amp_sample_ids,
        cell.labels = cell.labels,
        variants = dna_id,
        variant.filter = "unfiltered",
        vaf.mtx = vaf_mtx,
        vaf.mtx.cells = dna_sample_name,
        gt.mtx = gt_mtx,
        gt.mtx.cells = dna_sample_name,
        dp.mtx = dp_mtx,
        dp.mtx.cells = dna_sample_name,
        gq.mtx = gq_mtx,
        gq.mtx.cells = dna_sample_name,
        amps = amp_ca_ids,
        amp.normalize.method = amp.normalize.method,
        amp.mtx = amp_mtx,
        amp.mtx.cells = amp_sample_names,
        ploidy.mtx = NULL,
        proteins = proteins,
        protein.normalize.method = protein.normalize.method,
        protein.mtx = protein_mtx,
        protein.mtx.cells = protein.cells,
        backing_files = list(original = filepath),
        dna.variant.filter.mask = variant.filter.mask,
        dna_id_table = dna_id_table,
        cnv_id_table = cnv_id_table,
        dna_cell_table = dna_cell_table,
        cnv_metadata = cnv_metadata
    )

    obj
}


#' Harmonize columns (features) across multiple matrices
#'
#' Given a list of matrices (typically `DelayedArray`/`HDF5Array`) and their
#' corresponding column-name vectors, this helper aligns features across files
#' according to a specified policy:
#' \itemize{
#'   \item \code{"identical"} — require that all feature sets (and their order)
#'         are exactly the same across files; otherwise, throw an error.
#'   \item \code{"intersect"} — keep only the intersection of features shared by
#'         all files, and reorder each matrix consistently.
#' }
#'
#' @param mtx_list List of matrix-like objects to be harmonized (each with
#'   column names; typically `DelayedArray`/`HDF5Array`).
#' @param cols_list List of character vectors; column names for each matrix in
#'   \code{mtx_list}, in the same order.
#' @param label Character scalar used in error messages to identify the assay
#'   (e.g., \code{"VAF"}, \code{"GT"}, \code{"AMP"}, \code{"PROTEIN"}).
#' @param feature_policy Character. Either `"identical"` or `"intersect"`.
#'
#' @return
#' A list of matrix-like objects with harmonized columns. Returns \code{NULL}
#' if \code{mtx_list} or \code{cols_list} is \code{NULL}. Throws an error when
#' feature sets violate the selected policy (e.g., no common features).
#'
#' @examples
#' \dontrun{
#' mtxs <- list(m1, m2, m3)
#' cols <- list(colnames(m1), colnames(m2), colnames(m3))
#'
#' # Require identical features
#' harmonize_cols(mtxs, cols, label = "VAF", feature_policy = "identical")
#'
#' # Keep only common features
#' harmonize_cols(mtxs, cols, label = "VAF", feature_policy = "intersect")
#' }
#'
#' @keywords internal
harmonize_cols <- function(mtx_list, cols_list, label, feature_policy = c("identical", "intersect")) {
    feature_policy <- match.arg(feature_policy)
    if (is.null(mtx_list) || is.null(cols_list)) return(NULL)

    if (feature_policy == "identical") {
        ref  <- cols_list[[1]]
        same <- vapply(cols_list, function(x) identical(x, ref), logical(1))
        if (!all(same)) {
            stop(sprintf("[%s] Feature sets differ across files. Use feature_policy='intersect' to proceed.", label))
        }
        mtx_list
    } else { # "intersect"
        common <- Reduce(intersect, cols_list)
        if (length(common) == 0) {
            stop(sprintf("[%s] No common features to merge.", label))
        }
        lapply(mtx_list, function(m) m[, match(common, colnames(m)), drop = FALSE])
    }
}


#' Load and merge multiple HDF5 files from a directory (HDF5-backed)
#'
#' Reads all `.h5` files from a directory, converts each into an
#' HDF5-backed `ScIGMA_object`, and **merges them by cells** using `rbind`.
#' Features (columns) are aligned beforehand via \code{\link{harmonize_cols}}
#' (policy: `"identical"` or `"intersect"`).
#'
#' @param dir Character. Directory containing `.h5` files.
#' @param pattern Character. Regex to match file names (default: `"\\.h5$"`).
#' @param recursive Logical. Search subdirectories (default: `FALSE`).
#' @param omic.type `"DNA+protein"` or `"DNA"`. Passed to `loadH5_HDF5()`.
#' @param block.size Numeric. DelayedArray block size in bytes (default `1e8`).
#' @param feature_policy Character. `"identical"` (require same features/order across files)
#'   or `"intersect"` (keep only common features; reorder consistently). Default `"identical"`.
#' @param verbose Logical. Print progress messages.
#'
#' @return A merged `ScIGMA_object` (all matrices HDF5-backed).
#'
#' @examples
#' \dontrun{
#' obj <- loadH5_dir_HDF5("data/h5/", omic.type = "DNA+protein",
#'                        feature_policy = "intersect")
#' obj$print()
#' }
loadH5_dir_HDF5 <- function(dir,
                            pattern = "\\.h5$",
                            recursive = FALSE,
                            omic.type = c("DNA+protein","DNA"),
                            block.size = 1e8,
                            feature_policy = c("identical","intersect"),
                            verbose = TRUE) {
    stopifnot(dir.exists(dir))
    omic.type <- match.arg(omic.type)
    feature_policy <- match.arg(feature_policy)

    files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = recursive)
    if (length(files) == 0) stop("No .h5 files found in: ", dir)

    if (verbose) message(sprintf("[ScIGMA] Found %d file(s). Loading as HDF5-backed…", length(files)))
    objs <- lapply(files, function(f) {
        if (verbose) message("  - ", basename(f))
        loadH5_HDF5(f, sample.name = tools::file_path_sans_ext(basename(f)),
                    omic.type = omic.type, block.size = block.size)
    })

    pull <- function(slot) lapply(objs, function(o) o[[slot]])

    # --- Collect feature names per assay (columns) ---
    vaf_list <- pull("vaf.mtx")
    has_vaf  <- !sapply(vaf_list, is.null)
    if (!all(has_vaf)) stop("Some objects lack vaf.mtx; cannot merge.")
    v_cols <- lapply(vaf_list, colnames)

    amp_list <- pull("amp.mtx")
    has_amp  <- !sapply(amp_list, is.null)
    a_cols <- if (all(has_amp)) lapply(amp_list, colnames) else NULL

    prot_list <- pull("protein.mtx")
    has_prot  <- !sapply(prot_list, is.null)
    p_cols <- if (omic.type == "DNA+protein" && any(has_prot)) lapply(prot_list, colnames) else NULL

    # --- Harmonize columns (features) before row-binding (cells) ---
    vaf_list  <- harmonize_cols(vaf_list, v_cols, "VAF",      feature_policy = feature_policy)
    gt_list   <- harmonize_cols(pull("gt.mtx"),  lapply(pull("gt.mtx"),  colnames), "GT", feature_policy = feature_policy)
    dp_list   <- harmonize_cols(pull("dp.mtx"),  lapply(pull("dp.mtx"),  colnames), "DP", feature_policy = feature_policy)
    gq_list   <- harmonize_cols(pull("gq.mtx"),  lapply(pull("gq.mtx"),  colnames), "GQ", feature_policy = feature_policy)
    amp_list  <- harmonize_cols(amp_list, a_cols, "AMP",      feature_policy = feature_policy)
    prot_list <- if (omic.type == "DNA+protein") harmonize_cols(prot_list, p_cols, "PROTEIN", feature_policy = feature_policy) else NULL

    # --- Merge by cells (rows) ---
    vaf_all <- do.call(rbind, vaf_list)
    gt_all  <- do.call(rbind, gt_list)
    dp_all  <- do.call(rbind, dp_list)
    gq_all  <- do.call(rbind, gq_list)
    amp_all <- if (!is.null(amp_list))  do.call(rbind, amp_list)  else NULL
    prot_all<- if (!is.null(prot_list)) do.call(rbind, prot_list) else NULL

    # --- Merge metadata / labels ---
    concat_char <- function(slot) unlist(lapply(objs, function(o) o[[slot]]), use.names = FALSE)
    cell.ids    <- concat_char("cell.ids")
    cell.labels <- concat_char("cell.labels")

    variants <- colnames(vaf_all)
    amps     <- if (!is.null(amp_all)) colnames(amp_all) else character()
    proteins <- if (!is.null(prot_all)) colnames(prot_all) else if (omic.type=="DNA") "non-protein" else character()

    vaf_cells <- rownames(vaf_all)
    amp_cells <- if (!is.null(amp_all))  rownames(amp_all)  else character()
    prot_cells<- if (!is.null(prot_all)) rownames(prot_all) else character()

    amp.normalize.method     <- if (!is.null(amp_all)) "unnormalized" else character()
    protein.normalize.method <- if (!is.null(prot_all)) "unnormalized" else if (omic.type=="DNA") "non-protein" else character()

    meta.data <- vapply(files, function(f) tools::file_path_sans_ext(basename(f)), character(1))

    if (verbose) message("[ScIGMA] Row-bind merge done (HDF5-backed, lazy).")

    ScIGMA_object$new(
        meta.data = meta.data,
        cell.ids = cell.ids,
        cell.labels = cell.labels,
        variants = variants,
        variant.filter = "unfiltered",
        vaf.mtx = vaf_all,
        vaf.mtx.cells = vaf_cells,
        gt.mtx = gt_all,
        gt.mtx.cells = vaf_cells,
        dp.mtx = dp_all,
        dp.mtx.cells = vaf_cells,
        gq.mtx = gq_all,
        gq.mtx.cells = vaf_cells,
        amps = amps,
        amp.normalize.method = amp.normalize.method,
        amp.mtx = amp_all,
        amp.mtx.cells = amp_cells,
        ploidy.mtx = NULL,
        proteins = proteins,
        protein.normalize.method = protein.normalize.method,
        protein.mtx = prot_all,
        protein.mtx.cells = prot_cells,
        backing_files = list(original = files)
    )
}


loadH5_HDF5_biocond <- function(filepath, sample_name, omic_type = c("DNA+protein", "DNA"), block.size = 100e6){
    omic_type <- match.arg(omic_type)
    stopifnot(file.exists(filepath))

    options(DelayedArray.block.size = block.size)

    if (omic_type == "DNA+protein"){
        proteins <- .h5read_char(filepath, "/assays/protein_read_counts/ca/id")

        # SUPPRESSION de t() : la matrice sort nativement en Features x Cells
        protein_mtx <- .h5_delayed(filepath, "/assays/protein_read_counts/layers/read_counts")

        if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
            protein_samples <- "/assays/protein_read_counts/ra/label"
        } else {
            protein_samples <- "/assays/protein_read_counts/ra/sample_name"
        }
        protein_id_path <- "/assays/protein_read_counts/ra/barcode"
        protein_barcodes <- .h5read_char(filepath, protein_id_path)
        pr_ra_names <- .h5read_char(filepath, protein_samples)

        # INVERSION STRICTE : Lignes = Protéines, Colonnes = Cellules
        if (length(proteins) > 0 && length(protein_barcodes) > 0) {
            dimnames(protein_mtx) <- list(proteins, protein_barcodes)
        }

        protein.normalize.method <- "unnormalized"
        protein.cells <- setNames(pr_ra_names, nm = protein_barcodes)
    } else {
        message("DNA only: skip protein matrix")
        proteins <- "non-protein"
        protein.normalize.method <- "non-protein"
        protein_mtx <- NULL
        protein.cells <- character()
    }

    # ---- Check if label information exists ----
    ## dna_read_counts
    if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
        dna_read_count_samples <- "/assays/dna_read_counts/ra/label"
    } else {
        dna_read_count_samples <- "/assays/dna_read_counts/ra/sample_name"
    }

    amp_sample_ids <- .h5read_char(filepath, "/assays/dna_read_counts/ra/barcode") # Sample barcode
    amp_ca_ids <- .h5read_char(filepath, "/assays/dna_read_counts/ca/id") # Amplicon ids
    amp_sample_names <- setNames(amp_sample_ids, nm = .h5read_char(filepath, dna_read_count_samples))# Samples human readable names

    ## dna_variants
    if (.h5_has_path(filepath, "/assays/dna_variants/ra/label")){
        dna_variants_samples <- "/assays/dna_variants/ra/label"
    } else {
        dna_variants_samples <- "/assays/dna_variants/ra/sample_name"
    }
    dna_sample_name <- .h5read_char(filepath, dna_variants_samples)
    dna_sample_ids <- .h5read_char(filepath, "/assays/dna_variants/ra/barcode")
    dna_id <- .h5read_char(filepath, "/assays/dna_variants/ca/id")

    # ---- 1. Extraction et Forçage Dimensionnel (Features x Cells) ----
    # Tapestri stocke souvent en Cellules x Variants. t() transpose en Variants x Cellules.

    vaf_mtx <- .h5_delayed(filepath, "/assays/dna_variants/layers/AF")
    if (length(dna_id) > 0 && length(dna_sample_ids) > 0) dimnames(vaf_mtx) <- list(dna_id, dna_sample_ids)

    gt_mtx <- .h5_delayed(filepath, "/assays/dna_variants/layers/NGT")
    if (length(dna_id) > 0 && length(dna_sample_ids) > 0) {
        dimnames(gt_mtx) <- list(dna_id, dna_sample_ids)
    }

    dp_mtx <- .h5_delayed(filepath, "/assays/dna_variants/layers/DP")
    if (length(dna_id) > 0 && length(dna_sample_ids) > 0) {
        dimnames(dp_mtx) <- list(dna_id, dna_sample_ids)
    }

    gq_mtx <- .h5_delayed(filepath, "/assays/dna_variants/layers/GQ")
    if (length(dna_id) > 0 && length(dna_sample_ids) > 0) {
        dimnames(gq_mtx) <- list(dna_id, dna_sample_ids)
    }

    variant_filter_mask <- .h5_delayed(filepath, "/assays/dna_variants/layers/FILTER_MASK")
    if (length(dna_id) > 0 && length(dna_sample_ids) > 0) {
        dimnames(variant_filter_mask) <- list(dna_id, dna_sample_ids)
    }

    # Assurez-vous que amp_mtx est bien extrait avant cette ligne si ce n'est pas déjà le cas
    amp_mtx <- .h5_delayed(filepath, "/assays/dna_read_counts/layers/read_counts")
    if (length(amp_ca_ids) > 0 && length(amp_sample_ids) > 0) dimnames(amp_mtx) <- list(amp_ca_ids, amp_sample_ids)

    # ---- 2. Validation stricte des dimensions ----
    expected_dim <- c(length(dna_id), length(dna_sample_ids))
    raw_assays <- list(
        vaf = vaf_mtx,
        gt = gt_mtx,
        dp = dp_mtx,
        gq = gq_mtx,
        variant_filter_mask = variant_filter_mask
    )

    assays_list <- Filter(
        function(x) {
            !is.null(x) && dim(x)[1] == expected_dim[1] && dim(x)[2] == expected_dim[2]
        },
        raw_assays
    )

    if (length(assays_list) == 0) stop("Critical failure: DNA assays dimension mismatch.")

    # ---- 3. Tables S4Vectors (Lignes = Identifiants de la dimension) ----
    dna_id_table <- S4Vectors::DataFrame(
        dna_id = dna_id,
        amplicon = .h5read_char(filepath, "/assays/dna_variants/ca/amplicon"),
        chrom = .h5read_char(filepath, "/assays/dna_variants/ca/CHROM"),
        pos = .h5read_char(filepath, "/assays/dna_variants/ca/POS"),
        ado_gt_cells = .h5read_char(filepath, "/assays/dna_variants/ca/ado_gt_cells"), # RE-ADDED : Métrique QC
        ado_rate = .h5read_char(filepath, "/assays/dna_variants/ca/ado_rate"),         # RE-ADDED : Métrique QC
        filtered = .h5read_char(filepath, "/assays/dna_variants/ca/filtered"),
        row.names = dna_id
    )

    dna_cell_table <- S4Vectors::DataFrame(
        dna_barcode = .h5read_char(filepath, "/assays/dna_variants/ra/barcode"),
        dna_filtered = .h5read_char(filepath, "/assays/dna_variants/ra/filtered"),
        dna_sample_name = dna_sample_name, # FIX : Utilise la variable extraite du H5, pas l'argument global
        row.names = dna_sample_ids
    )

    cnv_id_table <- S4Vectors::DataFrame(
        dna_id = amp_ca_ids,
        chrom = .h5read_char(filepath, "/assays/dna_read_counts/ca/CHROM"),
        start_pos = as.numeric(.h5read_char(filepath, "/assays/dna_read_counts/ca/start_pos")),
        end_pos = as.numeric(.h5read_char(filepath, "/assays/dna_read_counts/ca/end_pos")),
        # Add here gene info
        row.names = amp_ca_ids
    )

    genome_v <- .h5read_char(filepath, "/assays/dna_read_counts/metadata/genome_version")
    if (length(genome_v) == 0) genome_v <- "hg19"

    # RE-ADDED : Extraction des métadonnées de pipeline (Tapestri, etc.)
    dna_variant_metadata <- .h5_read_metadata_group(filepath, "/assays/dna_variants/metadata")
    dna_read_counts_metadata <- .h5_read_metadata_group(filepath, "/assays/dna_read_counts/metadata")

    if (omic_type == "DNA+protein" && .h5_has_path(filepath, "/assays/protein_read_counts/metadata")) {
        protein_read_counts_metadata <- .h5_read_metadata_group(filepath, "/assays/protein_read_counts/metadata")
    } else {
        protein_read_counts_metadata <- list()
    }

    # ---- 4. Assemblage BioConductor ----
    sce_variants <- SingleCellExperiment::SingleCellExperiment(
        assays = assays_list,
        rowData = dna_id_table
    )

    experiment_list <- list(dna_variants = sce_variants)

    if (!is.null(amp_mtx)) {
        experiment_list$amplicons <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = amp_mtx),
            rowData = cnv_id_table
        )
    }

    if (omic_type == "DNA+protein" && exists("protein_mtx") && !is.null(protein_mtx)) {
        experiment_list$proteins <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = protein_mtx)
        )
    }

    # ---- 5. Fédérateur MAE ----
    mae_main <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experiment_list,
        colData = dna_cell_table,
        metadata = list(
            name = sample_name,
            variant_filter = "unfiltered",
            genome_version = genome_v,
            dna_variant_meta = dna_variant_metadata,         # RE-ADDED
            dna_read_counts_meta = dna_read_counts_metadata, # RE-ADDED
            protein_meta = protein_read_counts_metadata      # RE-ADDED
        )
    )
    # ---- 6. Instanciation R6 ----
    ScIGMA_object$new(
        mae = mae_main,
        backing_files = list(original = filepath),
        filetype = omic_type
    )
}

#' Sanitize MAE string metadata to remove Windows carriage returns and whitespaces
#'
#' @param mae MultiAssayExperiment object
#' @return Sanitized MultiAssayExperiment
#' @export
sanitize_mae_strings <- function(mae) {

    # Fonction interne ultra-rapide pour nettoyer un vecteur
    clean_char_vector <- function(x) {
        if (is.character(x)) {
            # Pulvérise \r et supprime les espaces invisibles aux extrémités
            return(trimws(gsub("\r", "", x), whitespace = "[\\h\\v]"))
        }
        return(x) # Ignore les vecteurs numériques/logiques
    }

    # Fonction interne pour nettoyer un DataFrame entier
    clean_df <- function(df) {
        if (is.null(df) || nrow(df) == 0) return(df)

        # Nettoyage du contenu des colonnes
        df[] <- lapply(df, clean_char_vector)

        # Nettoyage des index
        if (!is.null(rownames(df))) rownames(df) <- clean_char_vector(rownames(df))
        if (!is.null(colnames(df))) colnames(df) <- clean_char_vector(colnames(df))

        return(df)
    }

    message("Sanitizing MAE object (Removing \\r and phantom whitespaces)...")

    # 1. Purge du colData global (Métadonnées cliniques)
    SummarizedExperiment::colData(mae) <- clean_df(as.data.frame(SummarizedExperiment::colData(mae)))

    # 2. Purge itérative de chaque modalité (DNA, Protéines, Amplicons)
    for (exp_name in names(mae)) {
        se <- mae[[exp_name]]

        # Purge des rowData (là où se trouve ton erreur)
        SummarizedExperiment::rowData(se) <- clean_df(as.data.frame(SummarizedExperiment::rowData(se)))

        # Purge des rownames/colnames des matrices d'assay
        rownames(se) <- clean_char_vector(rownames(se))
        colnames(se) <- clean_char_vector(colnames(se))

        mae[[exp_name]] <- se
    }

    return(mae)
}

#' Sanitize protein marker names for Seurat compatibility
#'
#' @description Normalizes raw protein names by converting Greek characters
#'   and HTML entities into their Latin equivalents via C-level vectorization.
#'   Removes incompatible symbols to ensure complete Seurat interoperability.
#'
#' @param marker_names Character vector. Raw protein marker names.
#'
#' @return Character vector of sanitized alphanumeric names.
#'
#' @importFrom stringr str_replace_all
#' @export
sanitize_protein_markers <- function(marker_names) {
    greek_map <- c(
        "\u03B1" = "alpha",   "\u0391" = "Alpha",
        "&#945;" = "alpha",   "&#913;" = "Alpha",   "&alpha;" = "alpha",
        "\u03B2" = "beta",    "\u0392" = "Beta",
        "&#946;" = "beta",    "&#914;" = "Beta",    "&beta;"  = "beta",
        "\u03B3" = "gamma",   "\u0393" = "Gamma",
        "&#947;" = "gamma",   "&#915;" = "Gamma",   "&gamma;" = "gamma",
        "\u03B4" = "delta",   "\u0394" = "Delta",
        "&#948;" = "delta",   "&#916;" = "Delta",   "&delta;" = "delta",
        "\u03B5" = "epsilon", "\u0395" = "Epsilon",
        "&#949;" = "epsilon", "&#917;" = "Epsilon", "&epsilon;" = "epsilon",
        "\u03BA" = "kappa",   "\u039A" = "Kappa",
        "&#954;" = "kappa",   "&#922;" = "Kappa",   "&kappa;" = "kappa",
        "\u03BC" = "mu",      "\u039C" = "Mu",
        "&#956;" = "mu",      "&#924;" = "Mu",      "&mu;"    = "mu"
    )

    # Pipeline: Tidyverse standard, strictly avoiding for loops
    clean_names <- marker_names |>
        stringr::str_replace_all(greek_map) |>
        stringr::str_replace_all("[ /\\-]", "_") |>
        stringr::str_replace_all("[^a-zA-Z0-9_]", "") |>
        stringr::str_replace_all("_", "-")

    return(clean_names)
}
