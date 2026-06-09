future::plan(future::multisession, workers = 4)

#' The application server-side
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import gargoyle
#' @import future
#' @import promises
# #' @import bslib
#' @noRd
app_server <- function(input, output, session) {
    # Your application server logic
    options(shiny.maxRequestSize = 2048 * 1024^2)

    # Setup async
    # plan(multisession)
    plan(sequential)

    # --------------------------------------------------------------- #
    # --------------------------------------------------------------- #
    nav_hide(id = "tabs_analysis_inner", target = "tab_protein")
    nav_hide(id = "tabs_analysis_inner", target = "tab_multi_omics")

    # --------------------------------------------------------------- #
    # Source config
    source("R/utils-config.R")

    # --------------------------------------------------------------- #
    # Initialize lens object
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
         "umap_computed",
         "dna_clones_renamed",
         "compass_completed",
         "multiomics_annotated",
         "clusters_computed",
         "gating_updated")
    print(paste("App Server R6 ID:", data.table::address(ScIGMA_data)))
    # --------------------------------------------------------------- #
    # --------------------------------------------------------------- #
    on("dataLoaded", {

        # scenario_check <- isTRUE(ScIGMA_data$has_proteomics)
        # or
        # scenario_check <- ScIGMA_data$experiment_type == "multi_omics"

        show_advanced_tabs <- FALSE
        print("ScIGMA_data$filetype app_server")
        print(ScIGMA_data$filetype )
        if (ScIGMA_data$filetype == "DNA+protein") {
                show_advanced_tabs <- TRUE
        } else {
            show_advanced_tabs <- FALSE
        }

        if (show_advanced_tabs) {
            nav_show(id = "tabs_analysis_inner", target = "tab_protein")
            nav_show(id = "tabs_analysis_inner", target = "tab_multi_omics")
        } else {
            nav_hide(id = "tabs_analysis_inner", target = "tab_protein")
            nav_hide(id = "tabs_analysis_inner", target = "tab_multi_omics")
        }
    })

    # [ NODE_ACCESS : Analysis ]
    # ----------------------------------------------------- _
    # >> Overview _
    mod_analysis_overview_server("analysis_overview_1", ScIGMA_data)

    # >> DNA _
    mod_analysis_DNA_server("analysis_DNA_1", ScIGMA_data)

    # >> Protein _
    mod_analysis_Protein_server("analysis_Protein_1", ScIGMA_data)

    # >> CNV _
    mod_analysis_CNV_server("analysis_CNV_1", ScIGMA_data)

    # >> MultiOmics _
    mod_analysis_multiomics_server("analysis_multiomics_1", ScIGMA_data)

    # [ NODE_ACCESS : DOWWLOAD ]
    # ----------------------------------------------------- _
    # >> Download data _
    mod_download_panel_server("download_panel_1", ScIGMA_data)

    # [ NODE_ACCESS : HELP ]
    # ----------------------------------------------------- _
    mod_help_panel_server("help_panel_1", ScIGMA_data)

    # [ NODE_ACCESS : ABOUT ]
    # ----------------------------------------------------- _
    mod_about_panel_server("about_panel_1", ScIGMA_data)

}
