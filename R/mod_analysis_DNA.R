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
                        actionButton(ns("btn_filtrer"), "Apply",
                                     class = "btn-primary")
                    )
                )
            ),

            accordion_panel(
                "Rename clusters",
                uiOutput(ns("rename_cluster_ui")
                )
            ),
            br(),
            fluidRow(
                div(
                    plotOutput(ns("heatmap"),
                               height = "600px", width = "700px"),
                    align = "center")
            )
        )
    )
}

#' analysis_right_DNA Server Functions
#'
#' @noRd
#'
#'
#' @import InteractiveComplexHeatmap
#' @importFrom ComplexHeatmap draw
mod_analysis_DNA_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
        # Render DNA variants dataframe
        # Afficher la table de sélection
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            datatable(ScIGMA_data$variant.annotation|> arrange(desc(cell_proportion), desc(impact)),
                      selection = 'multiple',
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })
        # Récupérer les lignes sélectionnées seulement quand l'utilisateur clique
        observeEvent({input$btn_filtrer
            watch("dna_clones_renamed")}, {
                print("Rendering DNA heatmap")
                sel <- input$variant_selection_rows_selected
                heatmap_include_all_samples <- input$heatmap_include_all_samples
                if (length(sel) > 0) {
                    # récupérer les données correspondantes
                    ScIGMA_data$variants.filtered <- ScIGMA_data$variant.annotation[sel, ]
                    tmp_selected_variant <- sub(x = ScIGMA_data$variants.filtered$variant_id, pattern = "^([^:]+:)|^:", "")
                    # make heatmap
                    ht_res <- generate_dna_variant_heatmap(obj = ScIGMA_data,
                                                           selected_variants_df = ScIGMA_data$variants.filtered,
                                                           heatmap_include_all_samples = heatmap_include_all_samples)
                    ht <- ht_res$heatmap
                    ht <- draw(ht)
                    print("new heatmap rendered")
                    if(is.null(ScIGMA_data$dna_clones_renamed)){ # Initiate dna.clones if not present
                        # Set dna.clones
                        ScIGMA_data$dna.clones <- ht_res$clones
                    }
                    print("ScIGMA_data$dna.clones")
                    print(levels(ScIGMA_data$dna.clones))
                    # Trigger event
                    trigger("dnaVariant_selected")
                    # render heatmap
                    output$heatmap <- renderPlot({
                        ht
                    })
                }
            })


        observeEvent(watch("dnaVariant_selected"),
                     {
                         req(ScIGMA_data$dna.clones)
                         output$rename_cluster_ui <-  renderUI({
                             tagList(
                                 p("Here you can rename a cluster, select a cluster name on drop list on the left, write its new name in right box and click on 'Apply New Labels' button."),
                                 fluidRow(column(6,
                                                 pickerInput(
                                                     inputId = ns("rename_cluster_ui_oldName"),
                                                     label = "Style : primary",
                                                     choices = levels(ScIGMA_data$dna.clones),
                                                     options = pickerOptions(container = "body",
                                                                             style = "btn-outline-primary"),
                                                     width = "100%"
                                                 )
                                 ),
                                 column(6,
                                        textInput(ns("rename_cluster_ui_newName"),
                                                  "New name")
                                 )
                                 ),
                                 div(
                                     actionButton(
                                         inputId = ns("btn_update_cluster_labels"),
                                         label = "Apply New Labels",
                                         icon = icon("check"),
                                         class = "btn-primary w-100" # w-100 pour prendre toute la largeur
                                     ), style = "margin-top:10px;"
                                 )
                             )
                         })
                     })

        observeEvent(input$btn_update_cluster_labels,
                     {
                         req(ScIGMA_data$dna.clones)
                         # update dna.clones labels
                         oldName <- input$rename_cluster_ui_oldName
                         newName <- input$rename_cluster_ui_newName
                         levels(ScIGMA_data$dna.clones)[levels(ScIGMA_data$dna.clones) == oldName] <- newName

                         ScIGMA_data$dna_clones_renamed <- ScIGMA_data$dna.clones
                         trigger("dna_clones_renamed")
                     })


    })
}

## To be copied in the UI
# mod_analysis_DNA_ui("analysis_right_DNA_1")

## To be copied in the server
# mod_analysis_DNA_server("analysis_right_DNA_1", ScIGMA_data)
