#' lensR6
#'
#' @description Create LENS R6 object
#'
#' @return An R6 object
#'
#' @noRd
ScIGMA_object <- R6::R6Class(
    classname = "ScIGMA_object",
    public = list(
        meta.data = NULL,
        cell.ids = NULL,
        cell.labels = NULL,
        variants = NULL,
        variant.filter = NULL,
        vaf.mtx = NULL,
        vaf.mtx.cells = NULL,
        gt.mtx = NULL,
        gt.mtx.cells = NULL,
        dp.mtx = NULL,
        dp.mtx.cells = NULL,
        gq.mtx = NULL,
        gq.mtx.cells = NULL,
        amps = NULL,
        amp.normalize.method = NULL,
        amp.mtx = NULL,
        amp.mtx.cells = NULL,
        ploidy.mtx = NULL,
        proteins = NULL,
        protein.normalize.method = NULL,
        protein.mtx = NULL,
        protein.mtx.cells = NULL,
        # ---------------------------- #
        # methods
        initialize = function(data = NULL){
            self$data <- data
        },
        print = function(...){
            print(self$data)
        }
    )
)
