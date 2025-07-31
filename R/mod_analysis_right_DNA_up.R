#' analysis_right_DNA_up UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom DT DTOutput datatable renderDT
mod_analysis_right_DNA_up_ui <- function(id) {
    ns <- NS(id)
    tagList(
        accordion(
            id = ns("acc"),
            open = FALSE,

            accordion_panel(
                "Select DNA variants",
                DTOutput(ns("variant_selection")),
                br(),
                actionButton("btn_filtrer", "Filter",
                             class = "btn-primary")
            ),
            accordion_panel(
                title = "Section B",
                icon = bsicons::bs_icon("sliders"),
                "Section B content"
            )
        )
    )



}

#' analysis_right_DNA_up Server Functions
#'
#' @noRd
mod_analysis_right_DNA_up_server <- function(id, lensObject){
    moduleServer(id, function(input, output, session){
        ns <- session$ns


        df_filtrage <- as.data.frame(matrix(1:25, 5, 5))

        # Render DNA variants dataframe

        # Afficher la table de sélection
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            datatable(lensObject$variantAnnotation,
                      selection = 'multiple',
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })


    })
}

## To be copied in the UI
# mod_analysis_right_DNA_up_ui("analysis_right_DNA_up_1")

## To be copied in the server
# mod_analysis_right_DNA_up_server("analysis_right_DNA_up_1", lensObject)
