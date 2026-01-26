#' The application server-side
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import gargoyle
#' @import future
#' @import promises
#' @noRd
app_server <- function(input, output, session) {
    # Your application server logic
    options(shiny.maxRequestSize=2048*1024^2)

    # Setup async
    plan(multisession)

    # --------------------------------------------------------------- #
    # Source config
    source("R/utils-config.R")
    # --------------------------------------------------------------- #
    #Initialize lens object
    # ScIGMA_data <- ScIGMA_object$new()
    ScIGMA_data <- ScIGMA_object$new()
    # --------------------------------------------------------------- #
    # Initialize watchers
    init("initApp",
         "dataLoaded",
         "dnaVariant_filtered",
         "dnaVariant_selected",
         "CNV_filtered",
         "CNV_ploidy_computed",
         "CNV_ui_cnv_plot_parameters_rendered",
         "CNV_ui_cnv_plot_additionalParameters_rendered",
         "launch_umap",
         "umap_computed")
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

    # --------------------------------------------------------------- #
    # CNV
    mod_analysis_CNV_server("analysis_CNV_1", ScIGMA_data)

}

