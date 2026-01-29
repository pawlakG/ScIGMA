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
        navset_card_underline(
            nav_panel(
                "Variant selection",
                accordion(
                    id = ns("acc"),
                    open = FALSE,
                    accordion_panel(
                        "Select DNA variants",
                        DTOutput(ns("variant_selection")),
                        br(),
                        fluidRow(
                            actionButton(ns("btn_filtrer"), "Apply",
                                         class = "btn-primary")
                        )
                    ),
                    accordion_panel(
                        "DNA variant heatmap",
                        fluidRow(
                            div(
                                plotOutput(ns("dna_variant_heatmap"),
                                           height = "600px", width = "900px"),
                                align = "center"),
                            fluidRow(
                                column(6,
                                       materialSwitch(
                                           inputId = ns("heatmap_include_all_samples"),
                                           label = "Show missings ?",
                                           value = TRUE,
                                           status = "success"
                                       )
                                ),
                                column(6,
                                       actionButton(ns("btn_dna_variant_download"), "Download plot",
                                                    class = "btn-primary")
                                )
                            )
                        )
                    ),
                    accordion_panel(
                        "Rename clusters",
                        uiOutput(ns("rename_cluster_ui")
                        )
                    ),
                )
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
            print("ScIGMA_data$variant.annotation |> arrange(desc(cell_proportion), desc(impact)")
            print(ScIGMA_data$variant.annotation)
            datatable(ScIGMA_data$variant.annotation |> arrange(desc(cell_proportion), desc(impact)),
                      selection = 'multiple',
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })
        # Récupérer les lignes sélectionnées seulement quand l'utilisateur clique
        observeEvent({input$btn_filtrer
            input$heatmap_include_all_samples
            watch("dna_clones_renamed")}, {
                print("Rendering DNA heatmap")
                sel <- input$variant_selection_rows_selected
                heatmap_include_all_samples <- input$heatmap_include_all_samples
                if (length(sel) > 0) {
                    print(ScIGMA_data$variant.annotation)
                    # récupérer les données correspondantes
                    # ScIGMA_data$variants.filtered <- ScIGMA_data$variant.annotation[sel, ]
                    ScIGMA_data$variants.filtered <- ScIGMA_data$variant.annotation[ScIGMA_data$variant.annotation$row_id %in% sel, ]
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
                    # Trigger event
                    trigger("dnaVariant_selected")
                    # render heatmap
                    output$dna_variant_heatmap <- renderPlot({
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
                         levels <- oldName
                         names(levels) <- newName
                         # ScIGMA_data$dna.clones <- fct_recode(ScIGMA_data$dna.clones, newName = oldName)
                         ScIGMA_data$dna.clones <- fct_recode(ScIGMA_data$dna.clones, !!!levels)
                         ScIGMA_data$dna_clones_renamed <- ScIGMA_data$dna.clones
                         trigger("dna_clones_renamed")
                     })


    })
}

## To be copied in the UI
# mod_analysis_DNA_ui("analysis_right_DNA_1")

## To be copied in the server
# mod_analysis_DNA_server("analysis_right_DNA_1", ScIGMA_data)
