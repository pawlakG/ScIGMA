#' @importFrom R6 R6Class
#' @importFrom DelayedArray DelayedArray realize
#' @importFrom HDF5Array HDF5Array writeHDF5Array setHDF5DumpDir setHDF5DumpFile
#' @keywords classes
ScIGMA_object <- R6::R6Class(
    classname = "ScIGMA_object",
    public = list(
        #' @field data Named list holding all payload fields
        #' (matrices/vectors/metadata).
        data = NULL,

        #' @description
        #' Create a new `ScIGMA_object`.
        #'
        #' @param meta.data List of sample-level and assays metadata.
        #' @param cell.ids Character vector of cell IDs (barcodes).
        #' @param cell.ids.filtered Character vector of cell IDs (barcodes)
        #' after variant filtering.
        #' @param cell.labels Named character vector of human-readable cell
        #' labels.
        #' @param cell.labels.filtered Character vector of human-readable cell
        #' labels after variant filtering.
        #' @param variants Character vector of variant IDs.
        #' @param variants.filtered Character vector of variant ID after variant
        #' filterings.
        #' @param variant.filter Character vector of filtered variant IDs.
        #' @param dna.variant.filter.mask matrix of cell and variant filtering,
        #' extracted from
        #' sample.dna.layers.FILTER_MASK, sample.dna.col_attrs.filtered,
        #' sample.dna.row_attrs.filtered,
        #' n_passing_variants, n_passing_cells, and n_passing_variants_per_cell.
        #' @param dna.variant.filter.mask.filtered matrix of cell and variant
        #' filtering after variant filtering according to user parameters
        #' @param vaf.mtx `DelayedArray` (e.g., `HDF5Array`) for VAF data
        #' (or `NULL`).
        #' @param vaf.mtx.filtered `DelayedArray` (e.g., `HDF5Array`) for VAF
        #' data (or `NULL`) after variant filtering.
        #' @param vaf.mtx.cells Character vector of labels associated with
        #' `vaf.mtx`.
        #' @param gt.mtx `DelayedArray` for genotype data (or `NULL`).
        #' @param gt.mtx.filtered `DelayedArray` for genotype data (or `NULL`)
        #' after variant filtering.
        #' @param gt.mtx.cells Character vector of labels associated with
        #' `gt.mtx`.
        #' @param dp.mtx `DelayedArray` for read depth (or `NULL`).
        #' @param dp.mtx.filtered `DelayedArray` for read depth (or `NULL`)
        #' after variant filtering.
        #' @param dp.mtx.cells Character vector of labels associated with
        #' `dp.mtx`.
        #' @param gq.mtx `DelayedArray` for genotype quality (or `NULL`).
        #' @param gq.mtx.filtered `DelayedArray` for genotype quality
        #' (or `NULL`) after variant filtering.
        #' @param gq.mtx.cells Character vector of labels associated with
        #' `gq.mtx`.
        #' @param amps Character vector of amplicon IDs.
        #' @param amp.normalize.method Character scalar describing
        #' normalization used for amplicons.
        #' @param amp.mtx `DelayedArray` for amplicon read counts (or `NULL`).
        #' @param amp.mtx.filtered `DelayedArray` for amplicon read counts
        #' (or `NULL`) after variant filtering.
        #' @param amp.mtx.normalized `DelayedArray` for amplicon read
        #' normalized counts (or `NULL`).
        #' @param amp.mtx.cells Character vector of labels associated with
        #' `amp.mtx`.
        #' @param ploidy.mtx Optional `DelayedArray` for ploidy data
        #' (or `NULL`).
        #' @param proteins Character vector of protein IDs.
        #' @param protein.normalize.method Character scalar describing protein
        #' normalization.
        #' @param protein.mtx `DelayedArray` for protein expression/read counts
        #' (or `NULL`).
        #' @param protein.mtx.filtered `DelayedArray` for protein
        #' expression/read counts (or `NULL`) after variant filtering.
        #' @param protein.mtx.cells Character vector of labels associated
        #' with `protein.mtx`.
        #' @param backing_files List for bookkeeping HDF5 files
        #' (e.g., `list(original=..., realized=...)`).
        #' @param variant.annotation data.frame of annotated variants.
        #' @param selected.variants Selected variant during varaint selection step
        initialize = function(meta.data = list(),
                              cell.ids = character(),
                              cell.ids.filtered = character(),
                              cell.labels = character(),
                              cell.labels.filtered = character(),
                              variants = character(),
                              variants.filtered = character(),
                              variant.filter = character(),
                              dna.variant.filter.mask = matrix(),
                              dna.variant.filter.mask.filtered = matrix(),
                              vaf.mtx = NULL,
                              vaf.mtx.filtered = NULL,
                              vaf.mtx.cells = character(),
                              gt.mtx = NULL,
                              gt.mtx.filtered = NULL,
                              gt.mtx.cells = character(),
                              dp.mtx = NULL,
                              dp.mtx.filtered = NULL,
                              dp.mtx.cells = character(),
                              gq.mtx = NULL,
                              gq.mtx.filtered = NULL,
                              gq.mtx.cells = character(),
                              amps = character(),
                              amp.normalize.method = character(),
                              amp.mtx = NULL,
                              amp.mtx.filtered = NULL,
                              amp.mtx.normalized = NULL,
                              amp.mtx.cells = character(),
                              ploidy.mtx = NULL,
                              proteins = character(),
                              protein.normalize.method = character(),
                              protein.mtx = NULL,
                              protein.mtx.filtered = NULL,
                              protein.mtx.cells = character(),
                              backing_files = list(),
                              variant.annotation = NULL,
                              selected.variants = NULL) {
            # simple validator for HDF5-backed matrices
            .chk <- function(x, nm) {
                if (!is.null(x) && !inherits(x, "DelayedArray")) {
                    stop(sprintf("%s must be a DelayedArray/HDF5Array (HDF5-backed) or NULL.", nm),
                        call. = FALSE
                    )
                }
            }
            .chk(vaf.mtx, "vaf.mtx")
            .chk(gt.mtx, "gt.mtx")
            .chk(dp.mtx, "dp.mtx")
            .chk(gq.mtx, "gq.mtx")
            .chk(amp.mtx, "amp.mtx")
            .chk(ploidy.mtx, "ploidy.mtx")
            .chk(protein.mtx, "protein.mtx")

            # pack everything into the single `data` slot
            self$data <- list(
                meta.data = meta.data,
                cell.ids = cell.ids,
                cell.labels = cell.labels,
                cell.ids.filtered = cell.ids.filtered,
                cell.labels.filtered = cell.labels.filtered,
                variants = variants,
                variants.filtered = variants.filtered,
                variant.filter = variant.filter,
                dna.variant.filter.mask = dna.variant.filter.mask,
                dna.variant.filter.mask.filtered = dna.variant.filter.mask.filtered,
                vaf.mtx = vaf.mtx,
                vaf.mtx.filtered = vaf.mtx.filtered,
                vaf.mtx.cells = vaf.mtx.cells,
                gt.mtx = gt.mtx,
                gt.mtx.filtered = gt.mtx.filtered,
                gt.mtx.cells = gt.mtx.cells,
                dp.mtx = dp.mtx,
                dp.mtx.filtered = dp.mtx.filtered,
                dp.mtx.cells = dp.mtx.cells,
                gq.mtx = gq.mtx,
                gq.mtx.filtered = gq.mtx.filtered,
                gq.mtx.cells = gq.mtx.cells,
                amps = amps,
                amp.normalize.method = amp.normalize.method,
                amp.mtx = amp.mtx,
                amp.mtx.filtered = amp.mtx.filtered,
                amp.mtx.normalized = amp.mtx.normalized,
                amp.mtx.cells = amp.mtx.cells,
                ploidy.mtx = ploidy.mtx,
                proteins = proteins,
                protein.normalize.method = protein.normalize.method,
                protein.mtx = protein.mtx,
                protein.mtx.filtered = protein.mtx.filtered,
                protein.mtx.cells = protein.mtx.cells,
                backing_files = backing_files,
                variant.annotation = variant.annotation,
                selected.variants = selected.variants
            )
        },

        #' @description
        #' Print a compact summary of the object.
        print = function(...) {
            cat("ScIGMA_object (HDF5-backed) — data slot layout\n")
            cat("------------------------------------------------\n")
            cat("Metadata entries  : ", length(self$data$meta.data), "\n")
            cat("Cells filtered    : ", length(self$data$cell.ids.filtered), "\n")
            cat("Cells             : ", length(self$data$cell.ids), "\n")
            cat("Variants          : ", length(self$data$variants), "\n")
            cat("Variants filtered : ", length(self$data$variants.filtered), "\n")
            cat("Proteins          : ", length(self$data$proteins), "\n")
            cat("Proteins filtered : ", length(self$data$proteins.filtered), "\n")
            invisible(self)
        },

        #' @description
        #' Realize (persist) delayed matrices to a single HDF5 file.
        #'
        #' @param dir Character. Target directory.
        #' @param file Character. HDF5 file name (e.g., "ScIGMA_store.h5").
        #' @param chunkdim Optional integer vector for chunking (passed to writeHDF5Array()).
        #' @param level Integer compression level (0–9). Default 6.
        #' @return `self` (invisible).
        realize_all = function(dir, file = "ScIGMA_store.h5", chunkdim = NULL, level = 6) {
            HDF5Array::setHDF5DumpDir(dir)
            HDF5Array::setHDF5DumpFile(file)

            .real <- function(x, name) {
                if (is.null(x)) {
                    return(NULL)
                }
                if (!is.null(chunkdim)) {
                    return(HDF5Array::writeHDF5Array(
                        x,
                        filepath = file.path(dir, file), name = name,
                        chunkdim = chunkdim, level = level
                    ))
                } else {
                    return(DelayedArray::realize(x))
                }
            }

            self$data$vaf.mtx <- .real(self$data$vaf.mtx, "vaf_mtx")
            self$data$gt.mtx <- .real(self$data$gt.mtx, "gt_mtx")
            self$data$dp.mtx <- .real(self$data$dp.mtx, "dp_mtx")
            self$data$gq.mtx <- .real(self$data$gq.mtx, "gq_mtx")
            self$data$amp.mtx <- .real(self$data$amp.mtx, "amp_mtx")
            self$data$amp.mtx.normalized <- .real(self$data$amp.mtx.normalized, "amp_mtx")
            self$data$ploidy.mtx <- .real(self$data$ploidy.mtx, "ploidy_mtx")
            self$data$protein.mtx <- .real(self$data$protein.mtx, "protein_mtx")

            self$data$backing_files$realized <- file.path(dir, file)
            invisible(self)
        }
    ),

    # ---- Backward compatible access via active bindings ----
    active = list(
        meta.data = function(value) {
            if (missing(value)) {
                return(self$data$meta.data)
            }
            self$data$meta.data <- value
        },
        cell.ids = function(value) {
            if (missing(value)) {
                return(self$data$cell.ids)
            }
            self$data$cell.ids <- value
        },
        cell.labels = function(value) {
            if (missing(value)) {
                return(self$data$cell.labels)
            }
            self$data$cell.labels <- value
        },
        cell.ids.filtered = function(value) {
            if (missing(value)) {
                return(self$data$cell.ids.filtered)
            }
            self$data$cell.ids.filtered <- value
        },
        cell.labels.filtered = function(value) {
            if (missing(value)) {
                return(self$data$cell.labels.filtered)
            }
            self$data$cell.labels.filtered <- value
        },
        variants = function(value) {
            if (missing(value)) {
                return(self$data$variants)
            }
            self$data$variants <- value
        },
        variants.filtered = function(value) {
            if (missing(value)) {
                return(self$data$variants.filtered)
            }
            self$data$variants.filtered <- value
        },
        variant.filter = function(value) {
            if (missing(value)) {
                return(self$data$variant.filter)
            }
            self$data$variant.filter <- value
        },
        dna.variant.filter.mask = function(value) {
            if (missing(value)) {
                return(self$data$dna.variant.filter.mask)
            }
            self$data$dna.variant.filter.mask <- value
        },
        dna.variant.filter.mask.filtered = function(value) {
            if (missing(value)) {
                return(self$data$dna.variant.filter.mask.filtered)
            }
            self$data$dna.variant.filter.mask.filtered <- value
        },
        vaf.mtx = function(value) {
            if (missing(value)) {
                return(self$data$vaf.mtx)
            }
            self$data$vaf.mtx <- value
        },
        vaf.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$vaf.mtx.filtered)
            }
            self$data$vaf.mtx.filtered <- value
        },
        vaf.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$vaf.mtx.cells)
            }
            self$data$vaf.mtx.cells <- value
        },
        gt.mtx = function(value) {
            if (missing(value)) {
                return(self$data$gt.mtx)
            }
            self$data$gt.mtx <- value
        },
        gt.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$gt.mtx.filtered)
            }
            self$data$gt.mtx.filtered <- value
        },
        gt.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$gt.mtx.cells)
            }
            self$data$gt.mtx.cells <- value
        },
        dp.mtx = function(value) {
            if (missing(value)) {
                return(self$data$dp.mtx)
            }
            self$data$dp.mtx <- value
        },
        dp.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$dp.mtx.filtered)
            }
            self$data$dp.mtx.filtered <- value
        },
        dp.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$dp.mtx.cells)
            }
            self$data$dp.mtx.cells <- value
        },
        gq.mtx = function(value) {
            if (missing(value)) {
                return(self$data$gq.mtx)
            }
            self$data$gq.mtx <- value
        },
        gq.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$gq.mtx.filtered)
            }
            self$data$gq.mtx.filtered <- value
        },
        gq.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$gq.mtx.cells)
            }
            self$data$gq.mtx.cells <- value
        },
        amps = function(value) {
            if (missing(value)) {
                return(self$data$amps)
            }
            self$data$amps <- value
        },
        amp.normalize.method = function(value) {
            if (missing(value)) {
                return(self$data$amp.normalize.method)
            }
            self$data$amp.normalize.method <- value
        },
        amp.mtx = function(value) {
            if (missing(value)) {
                return(self$data$amp.mtx)
            }
            self$data$amp.mtx <- value
        },
        amp.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$amp.mtx.filtered)
            }
            self$data$amp.mtx.filtered <- value
        },
        amp.mtx.normalized = function(value) {
            if (missing(value)) {
                return(self$data$amp.mtx.normalized)
            }
            self$data$amp.mtx.normalized <- value
        },
        amp.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$amp.mtx.cells)
            }
            self$data$amp.mtx.cells <- value
        },
        ploidy.mtx = function(value) {
            if (missing(value)) {
                return(self$data$ploidy.mtx)
            }
            self$data$ploidy.mtx <- value
        },
        proteins = function(value) {
            if (missing(value)) {
                return(self$data$proteins)
            }
            self$data$proteins <- value
        },
        protein.normalize.method = function(value) {
            if (missing(value)) {
                return(self$data$protein.normalize.method)
            }
            self$data$protein.normalize.method <- value
        },
        protein.mtx = function(value) {
            if (missing(value)) {
                return(self$data$protein.mtx)
            }
            self$data$protein.mtx <- value
        },
        protein.mtx.filtered = function(value) {
            if (missing(value)) {
                return(self$data$protein.mtx.filtered)
            }
            self$data$protein.mtx.filtered <- value
        },
        protein.mtx.cells = function(value) {
            if (missing(value)) {
                return(self$data$protein.mtx.cells)
            }
            self$data$protein.mtx.cells <- value
        },
        backing_files = function(value) {
            if (missing(value)) {
                return(self$data$backing_files)
            }
            self$data$backing_files <- value
        },
        variant.annotation = function(value) {
            if (missing(value)) {
                return(self$data$variant.annotation)
            }
            self$data$variant.annotation <- value
        },
        selected.variants = function(value) {
            if (missing(value)) {
                return(self$data$selected.variants)
            }
            self$data$selected.variants <- value
        }
    )
)
