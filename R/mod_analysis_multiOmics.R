# NEW
# File: R/mod_analysis_multiomics.R

#' analysis_multiomics UI Function
#' @importFrom shiny NS tagList uiOutput
mod_analysis_multiomics_ui <- function(id) {
    ns <- NS(id)
    tagList(
        # Réceptacle 100% R (Identique au module Protéine)
        uiOutput(ns("multiomics_main_ui"))
    )
}

#' analysis_multiomics Server Functions
#' @importFrom bslib card navset_card_underline nav_panel
#' @importFrom plotly renderPlotly plot_ly layout config
mod_analysis_multiomics_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ---------------------------------------------------------
        # 1. Contrôleur d'affichage 100% R
        # ---------------------------------------------------------
        is_ready_flag <- shiny::reactiveVal(FALSE)

        observeEvent({
            watch("umap_computed")
            watch("dnaVariant_filtered")
            watch("dataLoaded")
        }, {
            # Condition : L'UMAP doit exister ET les clones ADN doivent être définis
            has_seurat <- !is.null(ScIGMA_data$seurat_object)
            has_umap <- has_seurat && !is.null(ScIGMA_data$seurat_object@reductions$umap)
            has_clones <- !is.null(ScIGMA_data$dna.clones)

            if (has_umap && has_clones) {
                is_ready_flag(TRUE)
            } else {
                is_ready_flag(FALSE)
            }
        }, ignoreNULL = FALSE, ignoreInit = FALSE)

        output$multiomics_main_ui <- shiny::renderUI({
            if (!isTRUE(is_ready_flag())) {
                # Cas 1 : Bloqué
                card(
                    br(), br(),
                    h3("Please complete DNA filtering and Protein UMAP first.",
                       style = "text-align: center; color: #7f8c8d;"),
                    br(), br()
                )
            } else {
                # Cas 2 : Autorisé
                navset_card_underline(
                    nav_panel(
                        "Multi-Omic Projection",
                        fluidRow(
                            column(3,
                                   grid_card(
                                       area = "sidebar",
                                       h3("Controls"),
                                       shinyWidgets::materialSwitch(
                                           inputId = ns("use_compass"),
                                           label = "Use COMPASS clones (if available)",
                                           value = FALSE,
                                           status = "primary"
                                       ),
                                       hr(),
                                       actionButton(ns("btn_auto_annotate"), "Auto-Annotate Cells",
                                                    icon = icon("magic"), class = "btn-info w-100"),
                                       helpText("Predicts cell types based on canonical surface markers.")
                                   )
                            ),
                            column(9,
                                   grid_card(
                                       area = "main",
                                       h3("Protein UMAP x DNA Clones"),
                                       plotlyOutput(ns("multi_umap"), height = "600px")
                                   )
                            )
                        )
                    ),
                    nav_panel(
                        "Clonal Distribution",
                        fluidRow(
                            column(12,
                                   grid_card(
                                       h3("DNA Clone distribution per Protein Phenotype"),
                                       plotlyOutput(ns("multi_barplot"), height = "500px")
                                   )
                            )
                        )
                    )
                )
            }
        })

        # ---------------------------------------------------------
        # 2. Préparation des données (Jointure Multi-Omique)
        # ---------------------------------------------------------
        reactive_multi_df <- reactive({
            req(is_ready_flag())
            watch("multiomics_annotated") # Écoute les changements de labels

            # A. UMAP & Clusters Protéines
            umap_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |> as.data.frame()
            umap_df$barcode <- rownames(umap_df)

            if ("predicted_celltype" %in% colnames(ScIGMA_data$seurat_object@meta.data)) {
                umap_df$ptn_cluster <- ScIGMA_data$seurat_object$predicted_celltype
            } else if ("seurat_clusters" %in% colnames(ScIGMA_data$seurat_object@meta.data)) {
                umap_df$ptn_cluster <- paste("Cluster", ScIGMA_data$seurat_object$seurat_clusters)
            } else {
                umap_df$ptn_cluster <- "Unclustered"
            }

            # B. Clones ADN
            if (isTRUE(input$use_compass) && !is.null(ScIGMA_data$cnv.active.clones)) {
                dna_vec <- ScIGMA_data$cnv.active.clones
            } else {
                dna_vec <- ScIGMA_data$dna.clones
            }

            dna_df <- data.frame(
                barcode = names(dna_vec),
                dna_clone = as.character(dna_vec)
            )

            # C. Fusion sécurisée (Left Join : Protège contre les dropouts ADN)
            final_df <- merge(umap_df, dna_df, by = "barcode", all.x = TRUE)
            final_df$dna_clone[is.na(final_df$dna_clone)] <- "Missing/Filtered"

            return(final_df)
        })

        # ---------------------------------------------------------
        # 3. Moteur de Rendu Graphique
        # ---------------------------------------------------------
        output$multi_umap <- renderPlotly({
            df <- reactive_multi_df()

            plot_ly(
                data = df,
                x = ~umap_1,
                y = ~umap_2,
                color = ~dna_clone,
                type = 'scattergl',
                mode = 'markers',
                marker = list(size = 6, opacity = 0.8),
                # FIX : Infobulles Q1
                text = ~paste("<b>Barcode:</b>", barcode,
                              "<br><b>Protein Phenotype:</b>", ptn_cluster,
                              "<br><b>DNA Genotype:</b>", dna_clone),
                hoverinfo = "text"
            ) %>%
                layout(
                    xaxis = list(title = "UMAP 1", zeroline = FALSE),
                    yaxis = list(title = "UMAP 2", zeroline = FALSE),
                    legend = list(title = list(text = "<b>DNA Clone</b>"))
                ) %>%
                config(displaylogo = FALSE)
        })

        output$multi_barplot <- renderPlotly({
            df <- reactive_multi_df()

            # Comptage des proportions
            count_df <- as.data.frame(table(df$ptn_cluster, df$dna_clone))
            colnames(count_df) <- c("Protein_Cluster", "DNA_Clone", "Count")
            count_df <- count_df[count_df$Count > 0, ]

            plot_ly(
                data = count_df,
                x = ~Protein_Cluster,
                y = ~Count,
                color = ~DNA_Clone,
                type = 'bar'
            ) %>%
                layout(
                    barmode = 'stack',
                    xaxis = list(title = "Protein Phenotype"),
                    yaxis = list(title = "Number of Cells"),
                    legend = list(title = list(text = "<b>DNA Clone</b>"))
                ) %>%
                config(displaylogo = FALSE)
        })

        # ---------------------------------------------------------
        # 4. Action : Auto-Annotation
        # ---------------------------------------------------------
        observeEvent(input$btn_auto_annotate, {
            req(ScIGMA_data$seurat_object)
            req("seurat_clusters" %in% colnames(ScIGMA_data$seurat_object@meta.data))

            assay_to_use <- ifelse("clr" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "clr", "counts")

            shiny::showNotification("Predicting cell types based on protein signatures...", id = "annot_notif")

            new_labels <- auto_annotate_clusters(ScIGMA_data$seurat_object, assay = assay_to_use)

            current_clusters <- as.character(ScIGMA_data$seurat_object$seurat_clusters)
            ScIGMA_data$seurat_object$predicted_celltype <- unname(new_labels[current_clusters])

            shiny::removeNotification("annot_notif")
            shiny::showNotification("Auto-annotation complete!", type = "message")

            # Force la mise à jour du data.frame réactif
            gargoyle::trigger("multiomics_annotated")
        })
    })
}

## To be copied in app_ui.R :
# mod_analysis_multiomics_ui("analysis_multiomics_1")
## To be copied in app_server.R :
# mod_analysis_multiomics_server("analysis_multiomics_1", ScIGMA_data)
