#' @importFrom R6 R6Class
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom DelayedArray DelayedArray realize
#' @importFrom HDF5Array setHDF5DumpDir setHDF5DumpFile
#' @keywords classes
ScIGMA_object <- R6::R6Class(
    classname = "ScIGMA_object",
    public = list(
        mae = NULL,             # Dataset actif
        mae_raw = NULL,         # Dataset d'origine (immuable)
        backing_files = list(),
        filetype = NULL,
        dna.clones = NULL,
        dna_clones_renamed = NULL,
        protein.filtered = NULL,
        variants.filtered = NULL,
        ploidy.mtx = NULL,
        cnv_dp_filtered = NULL,
        seurat_object = NULL,
        protein_gating_tree = list(),
        umaps = list(),
        cnv.active.clones = NULL,
        protein.gating_tree = NULL,
        dna.clones_pre_compass = NULL,
        dna_clone_colors = NULL,


        #' @description
        #' Renomme un clone et synchronise sa couleur
        update_dna_clone_names = function(old_name, new_name) {
            # 1. Mise à jour des niveaux du facteur dna.clones
            if (!is.null(self$dna.clones)) {
                levels(self$dna.clones)[levels(self$dna.clones) == old_name] <- new_name
                self$dna_clones_renamed <- self$dna.clones
            }

            # 2. Mise à jour de la palette de couleurs
            # FIX CRITIQUE : Suppression de la double référence (self$self)
            if (!is.null(self$dna_clone_colors)) {
                names(self$dna_clone_colors)[names(self$dna_clone_colors) == old_name] <- new_name
            }

            # Signal de mise à jour pour les graphiques réactifs
            gargoyle::trigger("dna_clones_renamed")
            invisible(self)
        },


        initialize = function(mae = NULL,
                              mae_raw = NULL,
                              backing_files = list(),
                              filetype = NULL,
                              seurat_object = NULL) {
            if (!is.null(mae) && !inherits(mae, "MultiAssayExperiment")) {
                stop("mae must be a MultiAssayExperiment object.")
            }
            self$mae <- mae
            self$mae_raw <- if (is.null(mae_raw)) mae else mae_raw
            self$backing_files <- backing_files
            self$filetype <- filetype
            self$seurat_object <- seurat_object
        },

        print = function(...) {
            cat("ScIGMA_object (HDF5-backed MultiAssayExperiment wrapper)\n")
            if (!is.null(self$mae)) {
                print(self$mae)
            } else {
                cat("Empty object.\n")
            }
            invisible(self)
        },

        reset_analysis = function() {
            if (!is.null(self$mae_raw)) {
                self$mae <- self$mae_raw
            }

            self$dna.clones <- NULL
            self$dna_clones_renamed <- NULL
            self$protein.filtered <- NULL
            self$variants.filtered <- NULL
            self$ploidy.mtx <- NULL
            self$cnv_dp_filtered <- NULL
            self$seurat_object <- NULL
            self$protein_gating_tree <- list()
            self$umaps <- list()
            self$cnv.active.clones <- NULL
            self$protein.gating_tree <- NULL
            dna.clones_pre_compass <- NULL
            self$dna_clone_colors <- NULL

            # Sécurité additionnelle sur les métadonnées
            if (!is.null(self$mae)) {
                S4Vectors::metadata(self$mae) <- list()
            }
            invisible(self)
        },

        realize_all = function(dir,
                               file = "scigma_store.h5",
                               chunkdim = NULL,
                               level = 6) {
            if (is.null(self$mae)) stop("No MAE object to realize.")
            HDF5Array::setHDF5DumpDir(dir)
            HDF5Array::setHDF5DumpFile(file)

            self$mae <- S4Vectors::endoapply(
                self$mae,
                DelayedArray::realize
            )
            self$backing_files$realized <- file.path(dir, file)
            invisible(self)
        }
    )
)
