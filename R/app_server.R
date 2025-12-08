#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import gargoyle
#' @noRd
app_server <- function(input, output, session) {
    # Your application server logic
    options(shiny.maxRequestSize=2048*1024^2)
    # --------------------------------------------------------------- #
    # Source config
    source("R/utils-config.R")
    # --------------------------------------------------------------- #
    #Initialize lens object
    # ScIGMA_data <- ScIGMA_object$new()
    ScIGMA_data <- ScIGMA_object$new()
    # --------------------------------------------------------------- #
    # Initialize watchers
    init("initApp")
    init("dataLoaded")
    init("dnaVariant_filtered")
    # --------------------------------------------------------------- #
    # Analysis
    # ---------------------------- #
    # Overview
    mod_analysis_overview_server("analysis_overview_1", ScIGMA_data)
    # ---------------------------- #
    # Analysis right up
    mod_analysis_DNA_server("analysis_DNA_1", ScIGMA_data)

    # --------------------------------------------------------------- #
    # Protein
    mod_analysis_Protein_server("analysis_Protein_1", ScIGMA_data)

}

