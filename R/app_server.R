#' The application server-side
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import gargoyle
#' @import future
#' @import promises
#' @import bslib
#' @noRd
app_server <- function(input, output, session) {
    # Your application server logic
    options(shiny.maxRequestSize = 2048 * 1024^2)

    # Setup async
    # plan(multisession)
    plan(sequential)

    # --------------------------------------------------------------- #
    # 1. Masquage initial des onglets (BSLIB)
    # --------------------------------------------------------------- #
    # On cache immédiatement les onglets avancés au lancement.
    # L'ID "tabs_analysis_inner" doit correspondre à celui défini dans navset_tab() (voir app_ui).
    nav_hide(id = "tabs_analysis_inner", target = "tab_protein")
    nav_hide(id = "tabs_analysis_inner", target = "tab_multi_omics")

    # --------------------------------------------------------------- #
    # Source config
    source("R/utils-config.R")

    # --------------------------------------------------------------- #
    # Initialize lens object
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
         "umap_computed",
         "dna_clones_renamed")
    print(paste("App Server R6 ID:", data.table::address(ScIGMA_data)))
    # --------------------------------------------------------------- #
    # 2. Logique de révélation conditionnelle (GARGOYLE)
    # --------------------------------------------------------------- #
    # On écoute l'événement déclenché par l'onglet Overview (ex: "dataLoaded" ou un nouveau trigger)
    on("dataLoaded", {
        # CONDITION : À adapter selon la structure exacte de ton objet R6.
        # Exemple : Si le type d'analyse contient "proteomics" ou si une checkbox est cochée.

        # scenario_check <- isTRUE(ScIGMA_data$has_proteomics)
        # ou
        # scenario_check <- ScIGMA_data$experiment_type == "multi_omics"

        # Placeholder pour la logique :
        show_advanced_tabs <- FALSE
        # Exemple hypothétique basé sur ta demande :
        if (ScIGMA_data$filetype == "DNA+protein") {
                show_advanced_tabs <- TRUE
        } else {
            show_advanced_tabs <- FALSE
        }

        if (show_advanced_tabs) {
            nav_show(id = "tabs_analysis_inner", target = "tab_protein")
            nav_show(id = "tabs_analysis_inner", target = "tab_multi_omics")
        } else {
            # Optionnel : recacher si l'utilisateur change d'avis et recharge des données simples
            nav_hide(id = "tabs_analysis_inner", target = "tab_protein")
            nav_hide(id = "tabs_analysis_inner", target = "tab_multi_omics")
        }
    })

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

    # --------------------------------------------------------------- #
    # MultiOmics
    mod_analysis_multiOmics_server("analysis_multiOmics_1", ScIGMA_data)

}
