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
#'
#' @examples
#' \dontrun{
#' # Minimal construction with HDF5-backed matrices (placeholders shown):
#' library(HDF5Array)
#' vaf <- HDF5Array("data.h5", "assays/dna_variants/layers/AF")
#' obj <- ScIGMA_object$new(
#'   meta.data = "SampleA",
#'   cell.ids = c("cell1","cell2"),
#'   cell.labels = c("cell1","cell2"),
#'   variants = c("var1","var2"),
#'   vaf.mtx = vaf
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
        dna.variant.filter.mask = variant.filter.mask
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

#' Protein matrix normalization
#'
#' Normalizes the protein matrix within a `ScIGMA_object` using the
#' centered log‑ratio (CLR) transform, inspired by the analogous function
#' in the optima package.
#'
#' @param ScIGMA_object A `ScIGMA_object` (R6).
#' @import compositions
#' @return The input `ScIGMA_object` with `protein.mtx` normalized and
#'  `protein.normalize.method` set to "normalized".
#' @keywords ScIGMA protein normalization
#' @export
#' @examples
#' \dontrun{
#' scigma <- normalizeProtein(scigma)
#' }

normalizeProtein <- function(ScIGMA_object) {
    # extract count matrix
    inputMatrix <- ScIGMA_object$protein.mtx.filtered |> as.matrix()
    # apply normalization CLR method
    ret <- (compositions::clr(inputMatrix + 1))

    ScIGMA_object$protein.mtx.filtered.normalized <- as.matrix(ret)
    ScIGMA_object$protein.normalize.method <- "CLR normalized"
    return(ScIGMA_object)
}


#' CNV normalization function
#'
#' Normalizes the CNV/amplicon count matrix to correct for cell‑wise and
#' amplicon‑wise effects, and updates `amp.normalize.method` from
#' "unnormalized" to "normalized". Adapted from the optima R package.
#'
#' @param ScIGMA_object A `ScIGMA_object` (R6).
#' @return The input `ScIGMA_object` with CNV matrix normalized and
#'   `amp.normalize.method` set to "normalized".
#' @keywords ScIGMA CNV normalization
#' @export
#' @examples
#' \dontrun{
#' scigma <- normalizeCNV(scigma)
#' }

normalizeCNV <- function(ScIGMA_object){
    cnv.mtx <- ScIGMA_object$amp.mtx

    # conduct normalize
    rowsum.threshold <- sort(rowSums(cnv.mtx), decreasing = TRUE)[10] / 10
    keep.tf <- rowSums(cnv.mtx) > rowsum.threshold
    cnv.mtx <- cnv.mtx / (apply(cnv.mtx, 1, mean) + 1)
    cnv.mtx <- t(t(cnv.mtx) / (apply(cnv.mtx[keep.tf,], 2, median) + 0.05))
    cnv.mtx <- cnv.mtx * 2

    # update object
    ScIGMA_object$amp.mtx <- cnv.mtx
    ScIGMA_object$amp.normalize.method <- "normalized"

    return(ScIGMA_object)
}

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

    print("keep_cells")
    print(length(keep_cells))
    print("obj$amp.mtx")
    print(dim(obj$amp.mtx))

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



#' Fonction to generate the dna_variant heatmap
#'
#' @import ComplexHeatmap
#' @import colorBlindness
# generate_dna_variant_heatmap <- function(obj, selected_variants, n_cluster = 6, min_prop_cluster = 0.01){
generate_dna_variant_heatmap <- function(obj, selected_variants_df, n_cluster = 6, min_prop_cluster = 0.01,
                                         heatmap_include_all_samples = TRUE){

    selected_variants <- sub(x = selected_variants_df$variant_id, pattern = "^([^:]+:)|^:", "")

    # Get genotype matrix for selected variants
    gt <- obj$gt.mtx[,selected_variants, drop = FALSE] |> as.matrix()
    # Get NGT_MASK matrix for selected variants
    msk <- obj$dna.variant.filter.mask.filtered[, selected_variants, drop = FALSE] |> as.matrix() != 0

    # Align rows and columns
    common_rows <- intersect(rownames(gt),
                             rownames(msk))
    common_cols <- intersect(colnames(gt),
                             colnames(msk))
    if (length(common_rows) == 0L || length(common_cols) == 0L) {
        stop("No overlap (rows/columns) between gt.mtx and mask.")
    }

    # Apply NGT_MASK matrix
    gt_filtered <- gt[common_rows,]

    # Set value 3 for filtered genotypes
    gt_filtered[cbind(row(msk)[msk],
                      col(msk)[msk])] <- 3

    # Take a subset where no row has values equal to 3 : meaning no missing info
    tmp_heamtap_matrix_filtered_noMissing <- gt_filtered[rowSums(gt_filtered == 3) == 0,]

    # Take rows with any values equal to 3
    tmp_heamtap_matrix_filtered_withMissing <- gt_filtered[rowSums(gt_filtered == 3) > 0,]
    tmp_heamtap_matrix_filtered_withMissing[tmp_heamtap_matrix_filtered_withMissing == 3] <- NA

    # perform clustering on matrix with no value equal to 3
    ## Hamming distance
    sample_dist_matrix <- proxy::dist(tmp_heamtap_matrix_filtered_noMissing, method = function(x,y) sum(x != y)) |>
        as.matrix()
    variant_dist_matrix <- proxy::dist(t(tmp_heamtap_matrix_filtered_noMissing), method = function(x,y) sum(x != y)) |>
        as.matrix()

    ## reorder rows and columns
    ### Get clusters
    hc_sample_noMissing <- hclust(as.dist(sample_dist_matrix), method = "ward.D2") |> cutree(k = sqrt(nrow(sample_dist_matrix))) |> sort()
    hc_variant_noMissing <- hclust(as.dist(variant_dist_matrix), method = "ward.D2")|> cutree(k = n_cluster) |> sort()
    ### Reassign too small samples clusters (compared to total samples) to a "small" group
    res_table_clusters <- table(hc_sample_noMissing)
    too_small_clusters <- names(res_table_clusters)[res_table_clusters/nrow(gt) < min_prop_cluster]
    hc_sample_noMissing[hc_sample_noMissing %in% too_small_clusters] <- "small"
    hc_sample_noMissing <- hc_sample_noMissing |> sort() |> as.factor()
    ### Reorder
    tmp_heamtap_matrix_filtered_noMissing_ordered <- tmp_heamtap_matrix_filtered_noMissing[names(hc_sample_noMissing),names(hc_variant_noMissing)]

    # rbind with samples with missing info if heatmap_include_all_samples is TRUE
    if (heatmap_include_all_samples){
        tmp_heamtap_matrix_filtered_complete_ordered <- rbind(tmp_heamtap_matrix_filtered_noMissing_ordered,
                                                              tmp_heamtap_matrix_filtered_withMissing)
        # Set labels / annotations
        heatmap_split_vector <- c(as.character(hc_sample_noMissing),
                                  rep("missing", nrow(tmp_heamtap_matrix_filtered_withMissing)))
    } else {
        tmp_heamtap_matrix_filtered_complete_ordered <- tmp_heamtap_matrix_filtered_noMissing_ordered
        heatmap_split_vector <- as.character(hc_sample_noMissing)
        # Remove "small" cluster
        tmp_heamtap_matrix_filtered_complete_ordered <- tmp_heamtap_matrix_filtered_complete_ordered[heatmap_split_vector != "small",]
        heatmap_split_vector <- heatmap_split_vector[heatmap_split_vector != "small"]
    }


    # plot heatmap
    ## Color palette
    dna_variant_colorPalette <- setNames(c("#E9E8EC", "#BAB7D0", "#3C2692"), nm = c("0","1","2"))
    ### Legend
    ## Annotation
    heatmap_true_levels <- levels(hc_sample_noMissing)[levels(hc_sample_noMissing)!="small"]
    annotationColor <- list(Cluster = setNames(c(colorBlindness::paletteMartin[-1][1:length(heatmap_true_levels)],"grey","#333333"),
                                               nm = c(heatmap_true_levels, "missing","small")))
    dna_variant_annotation <- rowAnnotation(Cluster = heatmap_split_vector,
                                            col =annotationColor,
                                            show_annotation_name = FALSE,
                                            show_legend = FALSE)

    # ---------------------------- #
    # Verify colnames before plotting
    idx_colnames <- sapply(colnames(tmp_heamtap_matrix_filtered_complete_ordered), function(x) {
        which(grepl(x, selected_variants_df$variant_id))

    }, simplify = FALSE, USE.NAMES = FALSE) |> unlist()

    new_colnames <- selected_variants_df$variant_id[idx_colnames]

    ## Heatmap
    print("colnames tmp_heamtap_matrix_filtered_complete_ordered")
    print(colnames(tmp_heamtap_matrix_filtered_complete_ordered))

    Heatmap(tmp_heamtap_matrix_filtered_complete_ordered,
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
                at = c("0","1","2"),
                labels = c("WT","HET","HOM"),
                grid_height = grid::unit(4, "cm")
            ))


}



#' Fonction to generate the dna_variant heatmap
#'
#' @import ComplexHeatmap
#' @import colorBlindness
# generate_dna_variant_heatmap <- function(obj, selected_variants, n_cluster = 6, min_prop_cluster = 0.01){
get_no_missing_dna_variant <- function(obj, selected_variants_df, n_cluster = 6, min_prop_cluster = 0.01 ){

    selected_variants <- sub(x = selected_variants_df$variant_id, pattern = "^([^:]+:)|^:", "")

    # Get genotype matrix for selected variants
    gt <- obj$gt.mtx[,selected_variants, drop = FALSE] |> as.matrix()
    # Get NGT_MASK matrix for selected variants
    msk <- obj$dna.variant.filter.mask.filtered[, selected_variants, drop = FALSE] |> as.matrix() != 0

    # Align rows and columns
    common_rows <- intersect(rownames(gt),
                             rownames(msk))
    common_cols <- intersect(colnames(gt),
                             colnames(msk))
    if (length(common_rows) == 0L || length(common_cols) == 0L) {
        stop("No overlap (rows/columns) between gt.mtx and mask.")
    }

    # Apply NGT_MASK matrix
    gt_filtered <- gt[common_rows,]

    # Set value 3 for filtered genotypes
    gt_filtered[cbind(row(msk)[msk],
                      col(msk)[msk])] <- 3

    # Take a subset where no row has values equal to 3 : meaning no missing info
    tmp_heamtap_matrix_filtered_noMissing <- gt_filtered[rowSums(gt_filtered == 3) == 0,]

    # Take rows with any values equal to 3
    tmp_heamtap_matrix_filtered_withMissing <- gt_filtered[rowSums(gt_filtered == 3) > 0,]
    tmp_heamtap_matrix_filtered_withMissing[tmp_heamtap_matrix_filtered_withMissing == 3] <- NA

    # perform clustering on matrix with no value equal to 3
    ## Hamming distance
    sample_dist_matrix <- proxy::dist(tmp_heamtap_matrix_filtered_noMissing, method = function(x,y) sum(x != y)) |>
        as.matrix()
    variant_dist_matrix <- proxy::dist(t(tmp_heamtap_matrix_filtered_noMissing), method = function(x,y) sum(x != y)) |>
        as.matrix()

    ## reorder rows and columns
    ### Get clusters
    hc_sample_noMissing <- hclust(as.dist(sample_dist_matrix), method = "ward.D2") |> cutree(k = sqrt(nrow(sample_dist_matrix))) |> sort()
    hc_variant_noMissing <- hclust(as.dist(variant_dist_matrix), method = "ward.D2")|> cutree(k = n_cluster) |> sort()
    ### Reassign too small samples clusters (compared to total samples) to a "small" group
    res_table_clusters <- table(hc_sample_noMissing)
    too_small_clusters <- names(res_table_clusters)[res_table_clusters/nrow(gt) < min_prop_cluster]
    hc_sample_noMissing[hc_sample_noMissing %in% too_small_clusters] <- "small"
    ## Remove "small variants"
    hc_sample_noMissing <- hc_sample_noMissing[hc_sample_noMissing != "small"]
    hc_sample_noMissing <- hc_sample_noMissing |> sort() |> as.factor()
    ### Reorder
    tmp_heamtap_matrix_filtered_noMissing_ordered <- tmp_heamtap_matrix_filtered_noMissing[names(hc_sample_noMissing),names(hc_variant_noMissing)]

    tmp_heamtap_matrix_filtered_noMissing_ordered

}


# --------------------------------------------------------------- #
#' Render ridge plot with ggplot2 and plotly
#'
#' @import ggplot2
#' @import ggridges
#' @import plotly
#' @import tidyr
render_protein_ridge_plot <- function(obj){
    print("render_protein_ridge_plot called")
    print("dim obj$protein.mtx.filtered.normalized")
    print(dim(obj$protein.mtx.filtered.normalized))
    print(head(obj$protein.mtx.filtered.normalized))
    tmp_data <- obj$protein.mtx.filtered.normalized |>
        as_tibble() |>
        pivot_longer(everything())
    tmp_plot <- tmp_data |>
        ggplot(aes(x=value, y=name, fill = name)) +
        geom_density_ridges() +
        theme_ridges() +
        theme(legend.position = "none")
    print("render plotlied plot")
    # print("tmp_plot")
    # print(tmp_plot)
    print(tmp_plot) |> plotly::ggplotly()
}

#' Génère un barplot de la proportion de protéines
#'
#' @import MatrixGenerics
#'
#' @param obj Un objet contenant la matrice protéique (obj$protein.mtx)
#' @param title Titre optionnel du graphique
#'
#' @return Un objet ggplot représentant la distribution relative des protéines
#' @export
#'
#' @examples
#' plot_protein_barplot(obj)
plot_protein_barplot <- function(obj, title = "Protein Percentage Distribution") {

    print("obj$protein.mtx")
    print(dim(obj$protein.mtx))
    print(head(obj$protein.mtx, 2))
    print(class(obj$protein.mtx))

    # Vérifications basiques
    if (is.null(obj$protein.mtx)) {
        stop("Object do not have a protein matrix (obj$protein.mtx).")
    }

    # Calcul des pourcentages
    protein_barplot_df <- data.frame(
        protein = colnames(obj$protein.mtx),
        percent = round(colSums2(obj$protein.mtx) / sum(obj$protein.mtx), 5) * 100
    )

    # Génération du plot
    protein_barplot <- protein_barplot_df |>
        ggplot(aes(x = protein, y = percent)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        xlab("Antibody") +
        ylab("Percentage") +
        ggtitle(title) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )

    print(protein_barplot) |> plotly::ggplotly()
}


# nsp_transform.R --------------------------------------------------------------
# Ré-implémentation R "NSP.transform()" d'après la doc officielle Mosaic
# Référence : https://missionbio.github.io/mosaic/pages/missionbio.demultiplex.protein.nsp.NSP.html
# Hypothèses clé (doc) :
#  - Le background et le signal sont linéairement dépendants du total de lectures par cellule
#  - On ajuste 2 modèles linéaires (fond & signal), puis on corrige : (x - y_bg_hat) / (y_sig_hat - y_bg_hat)
#  - Jitter gaussien (sd=jitter) centré 0 ; option scale pour oversequencing ; sous-échantillonnage (sample_size)
#
# Conventions :
#  - 'counts_mat' = matrice (cells x proteins) de lectures brutes (>=0)
#  - Noms snake_case, explicites, sans conflits avec fonctions natives
#  - Retourne une liste ; $normalized = matrice normalisée ; $models = coefficients fond/signal ; $scaling_factor

#----------------------------#
# Utilitaires internes
#----------------------------#

#' Estimation d'un facteur d'échelle global pour runs "oversequencés"
#' Heuristique : la doc NSP évoque 'scaling_factor(reads, jitter)' et une contrainte
#'   max_zero_read_cells : fraction max de cellules avec 0 reads après scaling.
#' Ici, on propose une stratégie robuste & simple, cohérente avec la doc.
estimate_scaling_factor <- function(total_reads,
                                    max_zero_read_cells = 0.05) {
    # total_reads : vecteur (longueur = nb cellules)
    # max_zero_read_cells : seuil (0–1) accepté de cellules à 0 après scaling
    #
    # 1) Si déjà peu de zéros, on n'impose pas de scaling.
    zero_frac <- mean(total_reads == 0)
    if (!is.finite(zero_frac) || zero_frac > max_zero_read_cells) {
        return(1.0)
    }
    # 2) Vise à réduire une longue traîne en utilisant un ratio médiane/p95
    #    (si p95 >> médiane, on rabaisse l'échelle)
    if (all(total_reads == 0)) return(1.0)
    tr_pos <- total_reads[total_reads > 0]
    if (length(tr_pos) < 10) return(1.0)

    med <- stats::median(tr_pos)
    p95 <- stats::quantile(total_reads, 0.95, names = FALSE)
    if (!is.finite(med) || !is.finite(p95) || p95 <= 0) return(1.0)

    # 3) Proposition : scale = med / p95 (<= 1 si p95 > med)
    proposal <- as.numeric(med / p95)
    proposal <- max(min(proposal, 1.0), 0.0)

    # 4) Vérifie que le scaling ne crée pas trop de zéros
    trial <- floor(total_reads * proposal)
    if (mean(trial == 0) <= max_zero_read_cells) {
        return(proposal)
    } else {
        return(1.0)
    }
}

#' Sélectionne indices "background" (bas) et "signal" (haut) pour un vecteur y
select_bg_sig_indices <- function(y, p_low = 0.10, p_high = 0.90,
                                  min_points = 20L) {
    # y : intensités d'une protéine sur un sous-ensemble de cellules
    # p_low/p_high : quantiles pour définir fond/signal
    # min_points : garde-fou si quantiles donnent trop peu de points
    q_low  <- stats::quantile(y, p_low,  names = FALSE, type = 7)
    q_high <- stats::quantile(y, p_high, names = FALSE, type = 7)

    bg_idx <- which(y <= q_low)
    sg_idx <- which(y >= q_high)

    # garde-fou : si trop peu de points, on prend extrêmes
    if (length(bg_idx) < min_points) {
        bg_idx <- order(y, decreasing = FALSE)[seq_len(min(min_points, length(y)))]
    }
    if (length(sg_idx) < min_points) {
        sg_idx <- order(y, decreasing = TRUE)[seq_len(min(min_points, length(y)))]
    }
    list(bg_idx = bg_idx, sg_idx = sg_idx, q_low = q_low, q_high = q_high)
}

#----------------------------#
# Fonction principale NSP
#----------------------------#

#' nsp_transform
#' Normalisation "Noise-corrected and Scaled Protein counts" (NSP) reconstituée
#'
#' @param counts_mat matrice (cells x proteins) de lectures brutes (>=0)
#' @param jitter sd du bruit gaussien ajouté avant modélisation (doc NSP)
#' @param scale facteur d'échelle global (<=1) ; si NULL -> estimation auto (doc NSP)
#' @param sample_size nb de cellules max pour ajuster les modèles (Inf = toutes)
#' @param random_state graine pseudo-aléatoire pour reproductibilité
#' @param p_low quantile bas pour "background" (par défaut 0.10 comme compromis robuste)
#' @param p_high quantile haut pour "signal" (par défaut 0.90)
#' @param max_zero_read_cells fraction max autorisée de cellules 0 après scaling auto (doc NSP)
#' @return list(normalized, models, scaling_factor)
#' @references Doc officielle NSP : missionbio.github.io/mosaic/pages/missionbio.demultiplex.protein.nsp.NSP.html
nsp_transform <- function(counts_mat,
                          jitter = 0.5,
                          scale = NULL,
                          sample_size = Inf,
                          random_state = NULL,
                          p_low = 0.10,
                          p_high = 0.90,
                          max_zero_read_cells = 0.05) {
    # Vérifs d'entrée
    if (!(is.matrix(counts_mat) || is.data.frame(counts_mat))) {
        stop("counts_mat doit être une matrice/data.frame (cells x proteins).")
    }
    counts_mat <- as.matrix(counts_mat)
    if (any(counts_mat < 0, na.rm = TRUE)) {
        stop("counts_mat doit contenir des valeurs >= 0.")
    }

    # Dimensions
    n_cells    <- nrow(counts_mat)  # nb de cellules
    n_proteins <- ncol(counts_mat)  # nb de protéines

    # Graine aléatoire (pour jitter & sampling)
    if (!is.null(random_state)) set.seed(as.integer(random_state))

    # Copie float
    x <- matrix(as.numeric(counts_mat), nrow = n_cells, ncol = n_proteins,
                dimnames = dimnames(counts_mat))

    # (1) Jitter gaussien centré 0 (doc NSP)
    if (isTRUE(jitter > 0)) {
        x <- x + matrix(stats::rnorm(n_cells * n_proteins, mean = 0, sd = jitter),
                        nrow = n_cells, ncol = n_proteins)
        x[x < 0] <- 0  # tronque à 0 (lectures ne peuvent pas être négatives)
    }

    # (2) Total reads par cellule (covariable principale des modèles)
    total_reads <- rowSums(x, na.rm = TRUE)

    # (3) Scaling global optionnel (doc NSP: scaling_factor pour oversequencing)
    if (is.null(scale)) {
        scale <- estimate_scaling_factor(total_reads, max_zero_read_cells = max_zero_read_cells)
    } else {
        # garde-fou : on force [0,1]
        scale <- max(min(as.numeric(scale), 1.0), 0.0)
    }
    if (!isTRUE(all.equal(scale, 1.0))) {
        # Applique l'échelle et recalcule total_reads
        x <- floor(x * scale)
        total_reads <- rowSums(x, na.rm = TRUE)
    }

    # (4) Sous-échantillonnage pour l'ajustement des modèles (doc NSP: sample_size)
    if (is.finite(sample_size) && sample_size < n_cells) {
        idx_fit <- sample.int(n_cells, size = sample_size)
    } else {
        idx_fit <- seq_len(n_cells)
    }
    t_fit <- total_reads[idx_fit]  # total reads des cellules utilisées pour fit

    # (5) Ajustement fond/signal par protéine + correction
    normalized <- matrix(NA_real_, nrow = n_cells, ncol = n_proteins,
                         dimnames = dimnames(counts_mat))
    model_list <- vector("list", n_proteins)
    names(model_list) <- colnames(counts_mat)

    eps <- 1e-6  # pour éviter division par ~0

    for (j in seq_len(n_proteins)) {
        # y = lectures (avec jitter/scale) pour la protéine j
        y       <- x[, j]
        y_fit   <- y[idx_fit]

        # 5a) Sélection des indices "background" (bas) et "signal" (haut) sur l'échantillon
        sel <- select_bg_sig_indices(y_fit, p_low = p_low, p_high = p_high, min_points = 20L)
        bg_idx_fit <- sel$bg_idx
        sg_idx_fit <- sel$sg_idx

        # 5b) Ajustements linéaires : y_bg ~ total_reads ; y_sig ~ total_reads (doc NSP)
        bg_lm <- stats::lm(y_fit[bg_idx_fit] ~ t_fit[bg_idx_fit])
        sg_lm <- stats::lm(y_fit[sg_idx_fit] ~ t_fit[sg_idx_fit])

        # 5c) Prédictions pour T de TOUTES les cellules
        y_bg_hat <- as.numeric(stats::coef(bg_lm)[1] + stats::coef(bg_lm)[2] * total_reads)
        y_sg_hat <- as.numeric(stats::coef(sg_lm)[1] + stats::coef(sg_lm)[2] * total_reads)

        # 5d) Garde-fous numériques (fond >= 0, signal >= fond + eps)
        y_bg_hat <- pmax(0, y_bg_hat)
        y_sg_hat <- pmax(y_bg_hat + eps, y_sg_hat)

        # 5e) Transformation NSP : (y − y_bg_hat) / (y_sg_hat − y_bg_hat), tronquée à >= 0
        num <- y - y_bg_hat
        den <- y_sg_hat - y_bg_hat
        val <- ifelse(den > eps, num / den, 0)
        val[val < 0] <- 0

        # 5f) Écrit la colonne normalisée
        normalized[, j] <- val

        # 5g) Stocke les modèles pour audit/reproductibilité
        model_list[[j]] <- list(
            protein      = colnames(counts_mat)[j],
            bg_coef      = stats::coef(bg_lm),
            sg_coef      = stats::coef(sg_lm),
            q_low        = sel$q_low,
            q_high       = sel$q_high,
            sample_size  = length(idx_fit),
            jitter       = jitter,
            scale        = scale
        )
    }

    # (6) Retour : matrice normalisée + méta
    list(
        normalized      = normalized,
        models          = model_list,
        scaling_factor  = scale
    )
}


