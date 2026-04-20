#' analysis_multiomics UI Function
#' @importFrom shiny NS tagList uiOutput plotOutput
#' @importFrom bslib card navset_card_underline nav_panel grid_card
mod_analysis_multiomics_ui <- function(id) {
    ns <- shiny::NS(id)

    shiny::tagList(
        shiny::uiOutput(ns("multiomics_main_ui"))
    )
}

#' analysis_multiomics Server Functions
mod_analysis_multiomics_server <- function(id, ScIGMA_data) {
    shiny::moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # [ NODE_ACCESS : "Prism-like" Plotly Specification]
        # ----------------------------------------------------- _
        prism_axis_style <- list(
            titlefont = list(size = 16, color = "black", family = "Arial"),
            tickfont = list(size = 14, color = "black", family = "Arial"),
            showline = TRUE,
            linewidth = 2,
            linecolor = "black",
            mirror = FALSE,
            ticks = "outside",
            tickwidth = 2,
            ticklen = 6,
            tickcolor = "black",
            showgrid = FALSE,
            zeroline = FALSE
        )

        # [ NODE_ACCESS : UI GENERATION ]
        # ----------------------------------------------------- _
        output$multiomics_main_ui <- shiny::renderUI({
            watch("umap_computed")
            # Check if UMAP has been calculated
            is_umap_ready <- !is.null(ScIGMA_data$seurat_object) && !is.null(ScIGMA_data$seurat_object[["umap"]])
            bslib::navset_card_underline(
                bslib::nav_panel(
                    "Protein UMAP x DNA",
                    if (any(!is_umap_ready, is.null(ScIGMA_data$dna.clones))){
                        tagList(fluidRow(h2("Please compute UMAP and DNA clones first")))
                    } else {
                        accordion(
                            id = ns("ptnUMAP_DNA_acc"),
                            open = FALSE,
                            accordion_panel(
                                "DNA clones",
                                shiny::fluidRow(
                                    shiny::column(3,
                                                  bslib::card(
                                                      shiny::h3("Projection Controls"),
                                                      shinyWidgets::materialSwitch(
                                                          inputId = ns("ptnUMAP_DNA_acc_gtMtx_choice"),
                                                          label = "Use COMPASS imputation",
                                                          value = FALSE,
                                                          status = "success"
                                                      )
                                                  )
                                    ),
                                    shiny::column(9,
                                                  bslib::card(
                                                      shiny::h3("Orthogonal Projection"),
                                                      plotlyOutput(ns("ptnUMAP_DNA_acc_dnaClones_plot"), height = "600px")
                                                  )
                                    )
                                )
                            ),
                            accordion_panel(
                                "DNA Variants",
                                shiny::fluidRow(
                                    shiny::column(3,
                                                  bslib::card(
                                                      shiny::h4("Variant Projection"),
                                                      shinyWidgets::pickerInput(
                                                          inputId = ns("selected_variant"),
                                                          label = "Targeted Variant",
                                                          # choices = rownames(ScIGMA_data$mae[["dna_variants"]]),,
                                                          choices = ScIGMA_data$variants.filtered$label,
                                                          options = list(`live-search` = TRUE)
                                                      ),
                                                      shinyWidgets::materialSwitch(
                                                          inputId = ns("use_compass_variant"),
                                                          label = "Use COMPASS imputation",
                                                          value = FALSE,
                                                          status = "success"
                                                      )
                                                  )
                                    ),
                                    shiny::column(9,
                                                  bslib::card(
                                                      plotly::plotlyOutput(ns("ptnUMAP_DNA_acc_dnaVariants_plot"), height = "600px")
                                                  )
                                    )
                                )
                            )
                        )
                    }
                ),
                bslib::nav_panel(
                    "UMAP Unsupervised Clusters x DNA",
                    if (any(!is_umap_ready, !"seurat_clusters" %in% colnames(ScIGMA_data$seurat_object@meta.data))){
                        tagList(fluidRow(h2("Please compute UMAP unsupervised clusters and DNA clones first")))
                    } else {
                        accordion(
                            id = ns("UMAPclusters_DNA_acc"),
                            open = FALSE,
                            accordion_panel(
                                "DNA clones",
                                bslib::card(
                                    shiny::plotOutput(ns("UMAPclusters_DNAclones_acc_barplot"), height = "600px")
                                )
                            ),
                            accordion_panel(
                                "DNA Variants",
                                bslib::card(
                                    shiny::plotOutput(ns("UMAPclusters_DNAvariants_acc_barplot"), height = "600px")
                                )
                            )
                        )
                    }
                ),
                bslib::nav_panel(
                    "Bi-plot Gates x DNA",
                    if (length(ScIGMA_data$protein_gating_tree) == 0) {
                        shiny::tagList(shiny::fluidRow(shiny::h2("Please define Biplot gates first")))
                    } else {
                        shiny::fluidRow(
                            shiny::column(3,
                                          bslib::card(
                                              shiny::h4("Population Profiling"),
                                              shinyWidgets::pickerInput(
                                                  inputId = ns("selected_biplot_pop"),
                                                  label = "Target Population",
                                                  choices = names(ScIGMA_data$protein_gating_tree$gates_list),
                                                  options = list(`live-search` = TRUE)
                                              ),
                                              shinyWidgets::materialSwitch(
                                                  inputId = ns("use_compass_biplot"),
                                                  label = "Use COMPASS imputation",
                                                  value = FALSE,
                                                  status = "success"
                                              )
                                          )
                            ),
                            shiny::column(9,
                                          bslib::card(
                                              plotly::plotlyOutput(ns("biplot_dna_distribution_plot"), height = "600px")
                                          )
                            )
                        )
                    }
                )
            )
        })

        # [ NODE_ACCESS : COMPUTATIONS ]
        # ----------------------------------------------------- _
        # >> Protein UMAP x DNA _
        # [!] Clones
        observeEvent(
            list(watch("umap_computed"),
                 watch("dnaVariant_selected")), {
                     output$ptnUMAP_DNA_acc_dnaClones_plot <- renderPlotly({
                         w <- Waiter$new(
                             id = ns("ptnUMAP_DNA_acc_dnaClones_plot"),
                             html = spin_3(),
                             color = transparent(0.5)
                         )
                         w$show()

                         # [!] Set parameters
                         # Clones
                         dna_clones <- ScIGMA_data$dna.clones
                         # Colors
                         n_clusters <- length(unique(ScIGMA_data$dna.clones[!ScIGMA_data$dna.clones %in% c("Missing","small")]))
                         dna_clones_colors <- setNames(c(viridis::viridis(n_clusters), NA, "grey"),
                                                       nm = c(unique(ScIGMA_data$dna.clones[!ScIGMA_data$dna.clones %in% c("Missing","small")]),
                                                              "Missing", "small"))
                         # [!] Set df
                         ptnUmap_dna_df <-ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                             as.data.frame()
                         ptnUmap_dna_df$dna_clones = ScIGMA_data$dna.clones[rownames(ptnUmap_dna_df)]

                         p <- plot_ly(data = ptnUmap_dna_df,
                                      x = ~umap_1,
                                      y = ~umap_2,
                                      type = 'scattergl',
                                      mode = 'markers',
                                      color = ~dna_clones,
                                      colors = dna_clones_colors,
                                      marker = list(size = 5, opacity = 0.8)) %>%
                             layout(
                                 plot_bgcolor = "white",
                                 paper_bgcolor = "white",
                                 xaxis = list(visible = FALSE),
                                 yaxis = list(visible = FALSE),
                                 legend = list(
                                     title = list(text = "<b>Cluster</b>", font = list(family = "Arial", color = "black")),
                                     font = list(family = "Arial", size = 12, color = "black")
                                 ),
                                 margin = list(l = 60, r = 30, b = 10, t = 30)
                             ) %>%
                             config(displaylogo = FALSE)

                         w$hide()
                         return(p)
                     })

                 },
            ignoreInit = TRUE)
        # [!] Variants _
        observeEvent(watch("dataLoaded"), {
            shiny::req(ScIGMA_data$mae)
            print("ScIGMA_data$variants.filtered")
            print(ScIGMA_data$variants.filtered)
            # available_variants <- rownames(ScIGMA_data$mae[["dna_variants"]])
            available_variants <- ScIGMA_data$variants.filtered$label
            shinyWidgets::updatePickerInput(
                session = session,
                inputId = "selected_variant",
                choices = available_variants
            )
        }, ignoreInit = TRUE)


        output$ptnUMAP_DNA_acc_dnaVariants_plot <- plotly::renderPlotly({
            # 1. Évaluation paresseuse stricte
            shiny::req("DNA Variants" %in% input$ptnUMAP_DNA_acc)
            watch("umap_computed")
            shiny::req(ScIGMA_data$seurat_object, input$selected_variant)

            # 2. Extraction des coordonnées
            umap_df <- as.data.frame(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings)
            umap_df$Barcode <- rownames(umap_df)

            print("input$selected_variant")
            print(input$selected_variant)
            print("input$use_compass_variant")
            print(input$use_compass_variant)

            # Retrieve variant id
            tmp_selected_variant <- rownames(ScIGMA_data$variants.filtered)[ScIGMA_data$variants.filtered$label == input$selected_variant]

            # 3. Extraction du génotype
            geno_df <- extract_variant_genotypes(
                mae_data = ScIGMA_data$mae,
                # variant_id = input$selected_variant,
                variant_id = tmp_selected_variant,
                use_compass = input$use_compass_variant
            )
            # 4. Jointure asymétrique
            plot_df <- merge(umap_df, geno_df, by = "Barcode", all.x = TRUE)

            # 5. Nomenclature clinique
            geno_map <- c(
                "0" = "WT",
                "1" = "HET",
                "2" = "HOM",
                "3" = "Missing/ADO",
                "Missing/ADO" = "Missing/ADO"
            )

            # plot_df$Variant_Genotype <- geno_map[as.character(plot_df$Variant_Genotype)]
            plot_df$Variant_Genotype[is.na(plot_df$Variant_Genotype)] <- "Unknown"

            plot_df$Variant_Genotype <- factor(
                plot_df$Variant_Genotype,
                levels = c("WT", "HET", "HOM", "Missing/ADO", "Unknown")
            )

            # 6. Verrouillage strict de la palette colorimétrique (Inspirée de Viridis)
            # Permet de figer les couleurs indépendamment des dropouts présents dans le lot
            variant_colors <- c(
                "WT" = "#440154",          # Violet profond
                "HET" = "#21918c",         # Sarcelle
                "HOM" = "#fde725",         # Jaune brillant
                "Missing/ADO" = "#e0e0e0", # Gris clair (Inerte)
                "Unknown" = "#000000"      # Noir (Anomalie)
            )

            # 7. Rendu Plotly (Hardware Accelerated via scattergl)
            p <- plotly::plot_ly(
                data = plot_df,
                x = ~umap_1,
                y = ~umap_2,
                type = 'scattergl', # CRITIQUE : Performance Single-Cell
                mode = 'markers',
                color = ~Variant_Genotype,
                colors = variant_colors,
                marker = list(size = 5, opacity = 0.8),
                # Tooltip interactif propre
                text = ~paste("Barcode:", Barcode, "<br>Genotype:", Variant_Genotype),
                hoverinfo = "text"
            ) |>
                plotly::layout(
                    plot_bgcolor = "white",
                    paper_bgcolor = "white",
                    xaxis = list(visible = FALSE),
                    yaxis = list(visible = FALSE),
                    legend = list(
                        title = list(text = "<b>Genotype</b>", font = list(family = "Arial", color = "black")),
                        font = list(family = "Arial", size = 12, color = "black")
                    ),
                    margin = list(l = 60, r = 30, b = 10, t = 30)
                ) |>
                plotly::config(displaylogo = FALSE)

            return(p)
        })

        # >> UMAP Unsupervised Clusters x DNA _
        # observeEvent({})

        # >> Bi-plot Gates x DNA Distribution _
        output$biplot_dna_distribution_plot <- plotly::renderPlotly({
            # Déclencheur : l'onglet doit être actif
            shiny::req(input$selected_biplot_pop, ScIGMA_data$variants.filtered)


            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            assay_colnames <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use) |> colnames()
            cel_barcode <- assay_colnames[ScIGMA_data$protein_gating_tree$gates_list[[input$selected_biplot_pop]]]

            print("cel_barcode")
            print(length(cel_barcode))
            print(head(cel_barcode))

            # Extraction des données via le helper
            plot_df <- compute_population_genotype_distribution(
                mae_data = ScIGMA_data$mae,
                # variant_ids = ScIGMA_data$variants.filtered$variant_id,
                variant_ids = rownames(ScIGMA_data$variants.filtered),
                # cell_barcodes = ScIGMA_data$protein_gating_tree$gates_list[[input$selected_biplot_pop]],
                cell_barcodes = cel_barcode,
                use_compass = input$use_compass_biplot
            ) |> as.data.frame()

            # Palette clinique (Gris pour ADO)
            variant_colors <- c(
                "WT" = "#440154", "HET" = "#21918c",
                "HOM" = "#fde725", "Missing/ADO" = "#e0e0e0"
            )


            # Merge plot_df with variant ScIGMA_data$variants.filtered

            variants.filtered_tmp <- ScIGMA_data$variants.filtered |> select(-variant_id) |>
                rownames_to_column("Variant_ID")

            plot_df_joined <- left_join(plot_df,variants.filtered_tmp,
                             by = "Variant_ID")

            plot_df_joined$Variant <- paste(plot_df_joined$protein, plot_df_joined$cdna, sep = ", ")


            plotly::plot_ly(
                data = plot_df_joined,
                x = ~Variant,
                y = ~Percentage,
                color = ~Variant_Genotype,
                colors = variant_colors,
                type = 'bar',
                text = ~paste0(round(Percentage, 1), "% (n=", Count, ")"),
                textposition = 'outside',
                textfont = list(size = 15, color = "black", family = "Arial"),
                constraintext = 'none',
                hoverinfo = 'text'
            ) |>
                plotly::layout(
                    barmode = 'group',
                    xaxis = c(list(title = "<b>DNA Variants</b>"), prism_axis_style),
                    yaxis = c(list(title = "<b>Frequency (%)</b>", range = c(0, 115)), prism_axis_style),
                    legend = list(title = list(text = "<b>Genotype</b>")),
                    margin = list(b = 100)
                ) |>
                plotly::config(displaylogo = FALSE)
        })


    })
}
