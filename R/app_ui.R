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
            # theme = bs_theme(bootswatch = "simplex"),
            # title = HTML("<h2><align=center>LENS portal</h2>"),
            title = HTML("ScIGMA portal"),
            selected = HTML("Welcome"),
            navbar_options = navbar_options(
                collapsible = TRUE),
            nav_panel(
                title =  HTML("Welcome"),
                mod_welcome_ui("welcome_1")
            ),
            nav_panel(
                title = HTML("Analysis"),
                grid_container(
                    layout = c(
                        "analysis_left analysis_right"
                    ),
                    row_sizes = c(
                        "1fr"
                    ),
                    col_sizes = c(
                        "0.40fr",
                        "1.60fr"
                    ),
                    gap_size = "10px",
                    grid_card(
                        area = "analysis_left",
                        full_screen = FALSE,
                        card_header(HTML("Start here:")),
                        mod_analysis_left_ui("analysis_left_1")
                    )
                    ,
                    grid_card(
                        area = "analysis_right",
                        card_body(
                            tabsetPanel(
                                nav_panel(title = "Overview.Preprocess",
                                          mod_analysis_right_overview_ui("analysis_right_overview_1")
                                ),
                                nav_panel(title = "DNA",
                                          mod_analysis_right_DNA_ui("analysis_right_DNA_1"),
                                ),
                                nav_panel(title = "CNV",
                                          mod_analysis_right_CNV_ui("analysis_right_CNV_1")),
                                nav_panel(title = "Protein",
                                          mod_analysis_right_Protein_ui("analysis_right_Protein_1")),
                                nav_panel(title = "Multi-omics",
                                          mod_analysis_right_multiOmics_ui("analysis_right_multiOmics_1"))
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
        # Add here other external resources
        # for example, you can add shinyalert::useShinyalert()
    )
}
