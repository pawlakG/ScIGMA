#' analysis_right_DNA UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom DT DTOutput datatable renderDT
mod_analysis_DNA_ui <- function(id) {
    ns <- NS(id)
    tagList(
        accordion(
            id = ns("acc"),
            open = FALSE,

            accordion_panel(
                "Select DNA variants",
                DTOutput(ns("variant_selection")),
                br(),

                fluidRow(
                    column(width = 6,
                           materialSwitch(
                               inputId = ns("heatmap_include_all_samples"),
                               label = "Include all samples ?",
                               value = TRUE,
                               status = "success"
                           )
                    ),
                    column(
                        width = 6,
                        actionButton(ns("btn_filtrer"), "Filter",
                                     class = "btn-primary")
                    )
                )
            ),

            fluidRow(
                plotOutput(ns("heatmap"), height = "1200px")
            )
        )
    )

}

#' analysis_right_DNA Server Functions
#'
#' @noRd
#'
#' @import InteractiveComplexHeatmap
mod_analysis_DNA_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
        # Render DNA variants dataframe
        # Afficher la table de sélection
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            datatable(ScIGMA_data$variant.annotation,
                      selection = 'multiple',
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })


        # Récupérer les lignes sélectionnées seulement quand l'utilisateur clique
        observeEvent(input$btn_filtrer, {
            sel <- input$variant_selection_rows_selected   # vecteur des index sélectionnéss
            heatmap_include_all_samples <- input$heatmap_include_all_samples
            print("heatmap_include_all_samples")
            print(heatmap_include_all_samples)
            if (length(sel) > 0) {
                # récupérer les données correspondantes
                print("test1")
                ScIGMA_data$variants.filtered <- ScIGMA_data$variant.annotation[sel, ]

                tmp_selected_variant <- sub(x = ScIGMA_data$variants.filtered$variant_id, pattern = "^([^:]+:)|^:", "")

                # make heatmap
                ht <- generate_dna_variant_heatmap(obj = ScIGMA_data,
                                                   selected_variants_df = ScIGMA_data$variants.filtered,
                                                   n_cluster = 3,
                                                   heatmap_include_all_samples = heatmap_include_all_samples)
                ht <- draw(ht)

                output$heatmap <- renderPlot({
                    ht
                })

            }
        })

    })


}

## To be copied in the UI
# mod_analysis_DNA_ui("analysis_right_DNA_1")

## To be copied in the server
# mod_analysis_DNA_server("analysis_right_DNA_1", ScIGMA_data)
