#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import bslib
#' @import gridlayout
#' @importFrom shinyWidgets radioGroupButtons actionBttn sliderTextInput virtualSelectInput pickerInput pickerOptions updatePickerInput
#' @noRd
app_ui <- function(request) {
    tagList(
        golem_add_external_resources(),

        use_waiter(),

        page_navbar(
            title = tagList(
                tags$img(src = "www/logo.png", height = "70px", style = "vertical-align: bottom; margin-right: -23px; margin-bottom: -5px;"),
                "cIGMA portal"
            ),
            selected = "welcome_tab",

            # 1. CORRECTION DU THEME
            # On utilise font_scale pour grossir le texte proprement
            # sans casser la hiérarchie visuelle de Lux.
            theme = bs_theme(
                bootswatch = "lux",
                font_scale = 1 # Augmente la taille globale de 10% (ajustez selon besoin)
            ),

            # 2. NETTOYAGE
            # J'ai retiré le bloc tags$style violent qui écrasait le thème.
            # Si vous avez besoin de CSS spécifique, ciblez des classes précises sans toucher aux h1-h6 globalement.

            navbar_options = navbar_options(collapsible = TRUE),

            nav_panel(
                title = HTML("Welcome"),
                value = "welcome_tab",
                mod_welcome_ui("welcome_1")
            ),

            nav_panel(
                title = HTML("Analysis"),
                value = "analysis_tab",

                navset_card_underline(
                    id = "tabs_analysis_inner",
                    nav_panel(
                        title = "Load and preprocess",
                        value = "tab_overview",
                        mod_analysis_overview_ui("analysis_overview_1")
                    ),
                    nav_panel(
                        title = "SNV",
                        value = "tab_dna",
                        mod_analysis_DNA_ui("analysis_DNA_1")
                    ),
                    nav_panel(
                        title = "CNV",
                        value = "tab_cnv",
                        mod_analysis_CNV_ui("analysis_CNV_1")
                    ),
                    nav_panel(
                        title = "Protein",
                        value = "tab_protein",
                        mod_analysis_Protein_ui("analysis_Protein_1")
                    ),
                    nav_panel(
                        title = "Multi-omics",
                        value = "tab_multi_omics",
                        mod_analysis_multiomics_ui("analysis_multiomics_1")
                    )
                )
            ),
            nav_panel(
                title = HTML("Download"),
                value = "download_tab",
                mod_download_panel_ui("download_panel_1")
            ),
            nav_panel(
                title = HTML("Help"),
                value = "help_tab",
                mod_help_panel_ui("help_panel_1")
            ),
            nav_panel(
                title = HTML("About"),
                value = "about_tab",
                mod_about_panel_ui("about_panel_1")
            )
        )
    )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
    add_resource_path(
        "www",
        app_sys("app/www")
    )

    tags$head(
        favicon(ext = "png"),
        bundle_resources(
            path = app_sys("app/www"),
            app_title = "ScIGMA"
        )
    )
}
