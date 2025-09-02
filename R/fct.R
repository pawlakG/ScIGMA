# ScIGMA — HDF5‑backed (DelayedArray/HDF5Array)
# ----------------------------------------------------
# Dépendances :
#   - BiocManager::install(c("DelayedArray", "HDF5Array", "BiocParallel", "rhdf5"))

suppressPackageStartupMessages({
    library(DelayedArray)
    library(HDF5Array)
    library(BiocParallel)
    library(rhdf5)
    library(R6)
})

# ----------------------------------------------------
# Helpers HDF5
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
    HDF5Array(filepath, path)
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
# Fonction d'import 100 % HDF5‑backed
# ----------------------------------------------------
# lit le fichier .h5 et construit un ScIGMA_object dont les matrices
# sont des HDF5Array/DelayedArray (aucune copie pleine en RAM).

#' Charger un fichier H5 en HDF5‑backed ScIGMA_object
#' @param filepath Chemin du fichier .h5
#' @param sample.name Nom d'échantillon (stocké dans meta.data)
#' @param omic.type "DNA+protein" ou "DNA"
#' @param block.size Taille cible des blocs (octets) pour DelayedArray (ex. 100e6)
#' @return Un objet ScIGMA_object (toutes matrices HDF5‑backed)
loadH5_HDF5 <- function(filepath, sample.name, omic.type = c("DNA+protein", "DNA"), block.size = 100e6){
    omic.type <- match.arg(omic.type)
    stopifnot(file.exists(filepath))

    options(DelayedArray.block.size = block.size)

    # ---- Protein (optionnel) ----
    if (omic.type == "DNA+protein"){
        proteins <- .h5read_char(filepath, "/assays/protein_read_counts/ca/id")
        protein_mtx <- t(.h5_delayed(filepath, "/assays/protein_read_counts/layers/read_counts"))
        if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
            protein_samples <- "/assays/protein_read_counts/ra/label"
        } else {
            protein_samples <- "/assays/protein_read_counts/ra/sample_name"
        }
        pr_ra_names <- .h5read_char(filepath, protein_samples)
        if (length(pr_ra_names) > 0) rownames(protein_mtx) <- pr_ra_names
        if (length(proteins) > 0) colnames(protein_mtx) <- proteins

        protein.normalize.method <- "unnormalized"
        protein.cells <- pr_ra_names
    } else {
        message("DNA only: skip protein matrix")
        proteins <- "non-protein"
        protein.normalize.method <- "non-protein"
        protein_mtx <- NULL
        protein.cells <- character()
    }


    # ---- Check is there is labels info ----
    ## dna_read_counts
    if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/label")){
        dna_read_count_samples <- "/assays/dna_read_counts/ra/label"
    } else {
        dna_read_count_samples <- "/assays/dna_read_counts/ra/sample_name"
    }
    ## dna_variants
    if (.h5_has_path(filepath, "/assays/dna_variants/ra/label")){
        dna_variants_samples <- "/assays/dna_variants/ra/label"
    } else {
        dna_variants_samples <- "/assays/dna_variants/ra/sample_name"
    }

    # ---- DNA variants ----
    dna_sample_name <- .h5read_char(filepath, dna_variants_samples)
    dna_id <- .h5read_char(filepath, "/assays/dna_variants/ca/id")

    vaf_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/AF"))
    if (length(dna_sample_name) > 0 || length(dna_id) > 0) dimnames(vaf_mtx) <- list(dna_sample_name, dna_id)

    gt_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/NGT"))
    if (length(dna_sample_name) > 0 || length(dna_id) > 0) dimnames(gt_mtx) <- list(dna_sample_name, dna_id)

    dp_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/DP"))
    if (length(dna_sample_name) > 0 || length(dna_id) > 0) dimnames(dp_mtx) <- list(dna_sample_name, dna_id)

    gq_mtx <- t(.h5_delayed(filepath, "/assays/dna_variants/layers/GQ"))
    if (length(dna_sample_name) > 0 || length(dna_id) > 0) dimnames(gq_mtx) <- list(dna_sample_name, dna_id)

    # ---- DNA read counts (amplicons) ----
    amp_mtx <- t(.h5_delayed(filepath, "/assays/dna_read_counts/layers/read_counts"))
    amp_ra_names <- .h5read_char(filepath, dna_read_count_samples)
    amp_ca_ids <- .h5read_char(filepath, "/assays/dna_read_counts/ca/id")
    if (length(amp_ra_names) > 0 || length(amp_ca_ids) > 0) dimnames(amp_mtx) <- list(amp_ra_names, amp_ca_ids)

    amp.normalize.method <- "unnormalized"

    # ---- Cell labels (logique de fallback identique au code d'origine) ----
    if (.h5_has_path(filepath, dna_read_count_samples)){
        cell.labels <- amp_ra_names
    } else if (length(dna_sample_name) == 0 ){
        cell.labels <- dna_sample_name
    } else if (.h5_has_path(filepath, protein_samples)){
        cell.labels <- pr_ra_names
    } else if (.h5_has_path(filepath, "/assays/dna_read_counts/ra/barcode")){
        barcode <- .h5read_char(filepath, "/assays/dna_read_counts/ra/barcode")
        cell.labels <- rep("unassigned", length(barcode))
    } else {
        cell.labels <- character()
    }

    # ---- Cell IDs (barcodes) ----
    cell.ids <- .h5read_char(filepath, "/assays/dna_read_counts/ra/barcode")

    # ---- Assemblage objet ----
    obj <- ScIGMA_object$new(
        meta.data = sample.name,
        cell.ids = cell.ids,
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
        amp.mtx.cells = amp_ra_names,
        ploidy.mtx = NULL,
        proteins = proteins,
        protein.normalize.method = protein.normalize.method,
        protein.mtx = protein_mtx,
        protein.mtx.cells = protein.cells,
        backing_files = list(original = filepath)
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





#' --------------------------------------------------------------- #
#' Function to get DNA clones
#'
#' @description Function inspired from optima R package.
#'
filter_variant <- function (ScIGMA_object, min.dp = 10, min.gq = 30, vaf.ref = 5,
                            vaf.hom = 95, vaf.het = 35, min.cell.pt = 50, min.mut.cell.pt = 1)
{
    vaf.mtx <- ScIGMA_object@vaf.mtx
    gt.mtx <- ScIGMA_object@gt.mtx
    dp.mtx <- ScIGMA_object@dp.mtx
    gq.mtx <- ScIGMA_object@gq.mtx
    dp.tf <- dp.mtx < min.dp
    gq.tf <- gq.mtx < min.gq
    vaf.ref.tf <- (vaf.mtx > vaf.ref) & (gt.mtx == 0)
    vaf.hom.tf <- (vaf.mtx < vaf.hom) & (gt.mtx == 2)
    vaf.het.tf <- (vaf.mtx < vaf.het) & (gt.mtx == 1)
    keep = !(dp.tf | gq.tf | vaf.ref.tf | vaf.hom.tf | vaf.het.tf)
    gt.mtx[!keep] <- 3
    vaf.mtx[gt.mtx == 3] <- -1
    num.cells <- nrow(keep)
    num.variants <- ncol(keep)
    cell.num.keep.tf <- colSums(apply(gt.mtx, 2, function(x) {
        x %in% 0:2
    })) > num.cells * min.cell.pt/100
    mut.cell.num.keep.tf <- colSums(apply(gt.mtx, 2, function(x) {
        x %in% 1:2
    })) > num.cells * min.mut.cell.pt/100
    variant.keep.tf <- cell.num.keep.tf & mut.cell.num.keep.tf
    v.names <- ScIGMA_object@variants
    v.names[variant.keep.tf]
    cell.variants.keep.tf <- rowSums(gt.mtx != 3) > num.variants *
        min.cell.pt/100
    c.names <- ScIGMA_object@cell.ids[cell.variants.keep.tf]
    if (ScIGMA_object@proteins[1] == "non-protein") {
        print("non-protein")
        my.protein.mtx <- ScIGMA_object@protein.mtx
    }
    else {
        my.protein.mtx <- ScIGMA_object@protein.mtx[cell.variants.keep.tf,
        ]
    }
    filtered.obj <- new("optima", meta.data = ScIGMA_object@meta.data,
                        cell.ids = ScIGMA_object@cell.ids[cell.variants.keep.tf],
                        cell.labels = ScIGMA_object@cell.labels[cell.variants.keep.tf],
                        variants = ScIGMA_object@variants[variant.keep.tf], variant.filter = "filtered",
                        vaf.mtx = ScIGMA_object@vaf.mtx[cell.variants.keep.tf, variant.keep.tf],
                        gt.mtx = ScIGMA_object@gt.mtx[cell.variants.keep.tf, variant.keep.tf],
                        dp.mtx = ScIGMA_object@dp.mtx[cell.variants.keep.tf, variant.keep.tf],
                        gq.mtx = ScIGMA_object@gq.mtx[cell.variants.keep.tf, variant.keep.tf],
                        amps = ScIGMA_object@amps, amp.normalize.method = ScIGMA_object@amp.normalize.method,
                        amp.mtx = ScIGMA_object@amp.mtx[cell.variants.keep.tf, ],
                        proteins = ScIGMA_object@proteins, protein.normalize.method = ScIGMA_object@protein.normalize.method,
                        protein.mtx = my.protein.mtx)
    cat("Number of cells removed: ")
    cat(length(cell.variants.keep.tf) - sum(cell.variants.keep.tf))
    cat("\nNumber of variants removed: ")
    cat(length(variant.keep.tf) - sum(variant.keep.tf))
    cat("\n")
    return(filtered.obj)
}

#' Protein matrix normalization
#'
#' The function normalizes protein matrix within an ScIGMA_object object using
#' CLR method inspired from the same fonction from optima package
#'
#' @param ScIGMA_object optima object.
#' @import compositions
#' @return An optima object with protein matrix being normalized and
#'  protein.normalize.method label updated to "normalized".
#' @keywords optima.obj
#' @export
#' @examples normalizeProtein(optima.object)

normalizeProtein <- function(ScIGMA_object) {
    # extract count matrix
    inputMatrix <- ScIGMA_object$protein.mtx
    # apply normalization CLR method
    ret <- (compositions::clr(inputMatrix + 1))

    ScIGMA_object$protein.mtx <- as.matrix(ret)
    ScIGMA_object$protein.normalize.method <- "normalized"
    return(ScIGMA_object)
}


#' CNV normalization function
#'
#' The function normalizes the CNV matrix to correct for column-wise and row-wise variation and
#' updates the ScIGMA_object object amp.normalize.method from "unnormalized" to "normalized".
#' This function is an adaptation of the same fonction from the optima R package
#'
#' @param optima.obj optima object.
#' @return optima object with normalized CNV and amp.normalize.method updated to "normalized".
#' @keywords optima.obj
#' @export
#' @examples normalizeCNV(optima.obj)

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

    filtered$amp.mtx.filtered      <- obj$amp.mtx[keep_cells, , drop = FALSE]

    filtered$protein.mtx.filtered  <- my.protein.mtx

    print("filtered$protein.mtx")
    print(dim(filtered$protein.mtx))

    # Summary
    removed_cells    <- length(cell.variants.keep.tf) - sum(cell.variants.keep.tf)
    removed_variants <- length(variant.keep.tf) - sum(variant.keep.tf)
    message("Number of cells removed: ", removed_cells)
    message("Number of variants removed: ", removed_variants)

    invisible(filtered)
}


