font_size_global <- "12px"

#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import bslib
#' @import gridlayout
#' @import shinythemes
#' @import shinyWidgets
#' @import shinydashboard
#' @noRd
app_ui <- function(request) {
    tagList(
        # Leave this function for adding external resources
        golem_add_external_resources(),

        # Your application UI logic
        page_navbar(
            fluid = TRUE,
            theme = bs_theme(bootswatch = "lux"),

            tags$head(
                tags$style(HTML(sprintf("
      body, .container, .well, .navbar, .form-control,
      h1, h2, h3, h4, h5, h6, p, label, .btn {
        font-size: %s !important;
      }
    ", font_size_global)))
            ),

            title = HTML("ScIGMA portal"),
            selected = "welcome_tab", # Utilisation d'ID explicites
            navbar_options = navbar_options(collapsible = TRUE),

            nav_panel(
                title = HTML("Welcome"),
                value = "welcome_tab",
                mod_welcome_ui("welcome_1")
            ),

            nav_panel(
                title = HTML("Analysis"),
                value = "analysis_tab",
                grid_container(
                    layout = c("analysis"),
                    row_sizes = c("1fr"),

                    grid_card(
                        area = "analysis",
                        card_body(
                            # Remplacement de tabsetPanel par navset_tab pour compatibilité bslib
                            # L'ID 'tabs_analysis_inner' est crucial pour le contrôle serveur
                            navset_tab(
                                id = "tabs_analysis_inner",

                                nav_panel(
                                    title = "Overview.Preprocess",
                                    value = "tab_overview",
                                    mod_analysis_overview_ui("analysis_overview_1")
                                ),

                                nav_panel(
                                    title = "DNA",
                                    value = "tab_dna",
                                    mod_analysis_DNA_ui("analysis_DNA_1")
                                ),

                                nav_panel(
                                    title = "CNV",
                                    value = "tab_cnv",
                                    mod_analysis_CNV_ui("analysis_CNV_1")
                                ),

                                # Ces onglets seront masqués au démarrage via le serveur
                                nav_panel(
                                    title = "Protein",
                                    value = "tab_protein",
                                    mod_analysis_Protein_ui("analysis_Protein_1")
                                ),

                                nav_panel(
                                    title = "Multi-omics",
                                    value = "tab_multi_omics",
                                    mod_analysis_multiOmics_ui("analysis_multiOmics_1")
                                )
                            )
                        )
                    )
                )
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
        favicon(),
        bundle_resources(
            path = app_sys("app/www"),
            app_title = "LENSapp"
        )
    )
}
