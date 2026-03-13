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
#' @importFrom forcats fct_recode
mod_analysis_DNA_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        # UPDATED
        # File: R/mod_analysis_overview.R (ou fichier contenant ce module)

        ns <- session$ns

        # 1. Render DNA variants dataframe
        output$variant_selection <- renderDT({
            watch("dnaVariant_filtered")
            req(ScIGMA_data$mae) # Sécurité

            # Extraction et tri (La "Vue")
            tmp_variant_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
                as.data.frame() |>
                dplyr::select(variant_id, gene, variant_type, gene_function, impact, clinvar, cell_proportion) |>
                dplyr::arrange(desc(cell_proportion), desc(impact))

            datatable(tmp_variant_annotation,
                      selection = 'multiple',
                      rownames = FALSE, # Désactivé car variant_id est déjà présent
                      options = list(pageLength = 5,
                                     lengthMenu = c(5, 10, 15)))
        })

        # 2. Récupérer les lignes sélectionnées
        observeEvent({
            input$btn_filtrer
            input$heatmap_include_all_samples
            watch("dna_clones_renamed")
        }, {
            print("Rendering DNA heatmap")
            sel_indices <- input$variant_selection_rows_selected
            heatmap_include_all_samples <- input$heatmap_include_all_samples

            if (length(sel_indices) > 0) {

                # RECONSTRUCTION DE LA VUE : Indispensable pour mapper les index de l'UI (sel_indices)
                # avec les véritables identifiants biologiques, car arrange() a mélangé les lignes.
                sorted_annotation <- SummarizedExperiment::rowData(ScIGMA_data$mae[["dna_variants"]]) |>
                    as.data.frame() |>
                    dplyr::select(variant_id, gene, variant_type, gene_function, impact, clinvar, cell_proportion) |>
                    dplyr::arrange(desc(cell_proportion), desc(impact))

                # Extraction sécurisée des variants sélectionnés
                selected_df <- sorted_annotation[sel_indices, , drop = FALSE]

                # Mise à jour de l'objet global (au cas où d'autres modules l'utilisent)
                ScIGMA_data$variants.filtered <- selected_df

                # Génération de la Heatmap
                ht_res <- generate_dna_variant_heatmap(
                    obj = ScIGMA_data,
                    selected_variants_df = selected_df,
                    heatmap_include_all_samples = heatmap_include_all_samples
                )

                ht <- ComplexHeatmap::draw(ht_res$heatmap)
                print("New heatmap rendered")

                # Gestion des clones
                if (is.null(ScIGMA_data$dna_clones_renamed)) {
                    ScIGMA_data$dna.clones <- ht_res$clones
                }

                # Déclenchement des événements avals
                trigger("dnaVariant_selected")

                # Affichage
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
