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
            # 1. Écoute active des signaux (Gargoyle)
            watch("umap_computed")
            watch("dnaVariant_selected")
            watch("dna_clones_renamed")
            watch("compass_completed")
            # Ajoute ces signaux s'ils existent dans tes modules Protéines
            watch("clusters_computed")
            watch("gating_updated")

            # 2. Vérification stricte des états (Ground Truth)
            # A. Les clones ADN existent-ils ?
            has_dna <- !is.null(ScIGMA_data$dna.clones_pre_compass) || !is.null(ScIGMA_data$dna.clones)

            # B. La UMAP protéique est-elle faite ?
            has_umap <- !is.null(ScIGMA_data$seurat_object) && !is.null(ScIGMA_data$seurat_object@reductions$umap)

            # C. Le clustering non-supervisé est-il fait ?
            # Vérifie si la colonne des clusters existe dans les métadonnées (adapte "seurat_clusters" si tu l'as nommée autrement)
            has_clusters <- has_umap && ("seurat_clusters" %in% colnames(ScIGMA_data$seurat_object@meta.data))

            # D. Le gating manuel est-il fait ?
            # Vérifie que la liste des sous-populations (gates) n'est pas vide
            has_gates <- !is.null(ScIGMA_data$protein_gating_tree) && (length(ScIGMA_data$protein_gating_tree) > 0)

            # 3. Condition par défaut : Aucun pré-requis n'est rempli
            if (!has_dna || (!has_umap && !has_clusters && !has_gates)) {
                return(
                    shiny::div(class = "text-center mt-5",
                               shiny::h4("Please compute UMAP unsupervised clusters and DNA clones first", class = "text-muted")
                    )
                )
            }

            # 4. Construction dynamique des onglets
            tabs <- list()

            compass_exists <- !is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)

            # Panel 1 : Protein UMAP x DNA
            if (has_dna && has_umap) {
                # 2. Construction conditionnelle de l'interrupteur
                compass_switch_ui_ptnUMAP_DNA <- shinyWidgets::materialSwitch(
                    inputId = ns("ptnUMAP_DNA_acc_gtMtx_choice"),
                    label = if(compass_exists) {
                        "Use COMPASS imputed matrix ?"
                    } else {
                        "Run MCMC first"
                    },
                    value = compass_exists,
                    status = "success"
                )

                compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants <- shinyWidgets::materialSwitch(
                    inputId =  ns("use_compass_variant"),
                    label = if(compass_exists) {
                        "Use COMPASS imputed matrix ?"
                    } else {
                        "Run MCMC first"
                    },
                    value = compass_exists,
                    status = "success"
                )

                if (!compass_exists) {
                    compass_switch_ui_ptnUMAP_DNA <- shiny::div(
                        style = "pointer-events: none; opacity: 0.5;",
                        compass_switch_ui_ptnUMAP_DNA
                    )

                    compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants <- shiny::div(
                        style = "pointer-events: none; opacity: 0.5;",
                        compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants
                    )
                } else {
                    compass_switch_ui_ptnUMAP_DNA <- shiny::div(
                        compass_switch_ui_ptnUMAP_DNA
                    )

                    compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants <- shiny::div(
                        compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants
                    )
                }

                tabs[[length(tabs) + 1]] <- bslib::nav_panel(
                    title = "Protein UMAP x DNA",


                    accordion(
                        id = ns("ptnUMAP_DNA_acc"),
                        open = FALSE,
                        accordion_panel(
                            "DNA clones",
                            shiny::fluidRow(
                                shiny::column(3,
                                              bslib::card(
                                                  shiny::h3("Projection Controls"),
                                                  shiny::div(compass_switch_ui_ptnUMAP_DNA, align = "left")
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

                                                  shiny::div(compass_switch_ui_ptnUMAP_DNA_acc_dnaVariants,
                                                             align = "left")
                                                  # shinyWidgets::materialSwitch(
                                                  #     inputId = ns("use_compass_variant"),
                                                  #     label = "Use COMPASS imputation",
                                                  #     value = FALSE,
                                                  #     status = "success"
                                                  # )
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
                )
            }

            # Panel 2 : UMAP Unsupervised Clusters x DNA
            if (has_dna && has_clusters) {
                tabs[[length(tabs) + 1]] <- bslib::nav_panel(
                    title = "UMAP Unsupervised Clusters x DNA",
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
                )
            }

            # Panel 3 : Bi-plot Gates x DNA
            if (has_dna && has_gates) {

                # 2. Construction conditionnelle de l'interrupteur
                compass_switch_ui_use_compass_biplot <- shinyWidgets::materialSwitch(
                    inputId = ns("use_compass_biplot"),
                    label = if(compass_exists) {
                        "Use COMPASS imputed matrix ?"
                    } else {
                        "Run MCMC first"
                    },
                    value = compass_exists,
                    status = "success"
                )

                if (!compass_exists) {
                    compass_switch_ui_use_compass_biplot <- shiny::div(
                        style = "pointer-events: none; opacity: 0.5;",
                        compass_switch_ui_use_compass_biplot
                    )
                } else {
                    compass_switch_ui_use_compass_biplot <- shiny::div(
                        compass_switch_ui_use_compass_biplot
                    )
                }

                tabs[[length(tabs) + 1]] <- bslib::nav_panel(
                    title = "Bi-plot Gates x DNA",
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
                                          shiny::div(compass_switch_ui_use_compass_biplot, align = "left")
                                      )
                        ),
                        shiny::column(9,
                                      bslib::card(
                                          plotly::plotlyOutput(ns("biplot_dna_distribution_plot"), height = "600px")
                                      )
                        )
                    )
                )
            }

            # 5. Rendu final (Si au moins un onglet est valide, on l'affiche)
            if (length(tabs) > 0) {
                do.call(bslib::navset_card_underline, c(list(id = ns("multiomics_tabs")), tabs))
            } else {
                return(
                    shiny::div(class = "text-center mt-5",
                               shiny::h4("Please compute UMAP unsupervised clusters and DNA clones first", class = "text-muted")
                    )
                )
            }
        })

        # [ NODE_ACCESS : COMPUTATIONS ]
        # ----------------------------------------------------- _
        # >> Protein UMAP x DNA _
        # [!] Clones
        observeEvent(
            list(watch("umap_computed"),
                 watch("dna_clones_renamed"),
                 watch("compass_completed")), {

                     output$ptnUMAP_DNA_acc_dnaClones_plot <- renderPlotly({
                         w <- Waiter$new(id = ns("ptnUMAP_DNA_acc_dnaClones_plot"), html = spin_3(), color = transparent(0.5))
                         w$show()

                         # 1. Sélection dynamique du vecteur selon le slider
                         use_compass <- isTRUE(input$ptnUMAP_DNA_acc_gtMtx_choice)

                         print("ScIGMA_data$dna.clones_pre_compass")
                         print(table(ScIGMA_data$dna.clones_pre_compass))
                         print("dna_clones_to_use <- ScIGMA_data$dna.clones")
                         print(table(dna_clones_to_use <- ScIGMA_data$dna.clones))

                         if (use_compass) {
                             # Sécurité : si COMPASS n'a pas été lancé, on force le retour au brut
                             if (is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)) {
                                 shiny::showNotification("COMPASS non exécuté. Affichage des clones bruts.", type = "warning")
                                 print("using pre-compass clones")
                                 dna_clones_to_use <- ScIGMA_data$dna.clones_pre_compass
                             } else {
                                 print("using post-compass clones")
                                 dna_clones_to_use <- ScIGMA_data$dna.clones
                             }
                         } else {
                             print("using pre-compass clones")
                             dna_clones_to_use <- ScIGMA_data$dna.clones_pre_compass
                         }

                         shiny::req(dna_clones_to_use)

                         # 2. Préparation du DataFrame
                         ptnUmap_dna_df <- as.data.frame(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings)
                         ptnUmap_dna_df$dna_clones <- as.character(dna_clones_to_use[rownames(ptnUmap_dna_df)])
                         ptnUmap_dna_df$dna_clones[is.na(ptnUmap_dna_df$dna_clones)] <- "Missing"

                         # 3. Rendu avec palette synchronisée
                         p <- plot_ly(data = ptnUmap_dna_df,
                                      x = ~umap_1, y = ~umap_2,
                                      type = 'scattergl', mode = 'markers',
                                      color = ~dna_clones,
                                      colors = ScIGMA_data$dna_clone_colors, # Utilisation de la palette R6
                                      marker = list(size = 5, opacity = 0.8),
                                      text = ~paste("<b>DNA clones</b>:", dna_clones),
                                      hoverinfo = "text") %>%
                             layout(
                                 xaxis = c(list(title = "<b>UMAP 1</b>"), prism_axis_style),
                                 yaxis = c(list(title = "<b>UMAP 2</b>"), prism_axis_style),
                                 legend = list(title = list(text = "<b>Clones</b>"))
                             ) %>%
                             config(displaylogo = FALSE)

                         w$hide()
                         return(p)
                     })
                 }, ignoreInit = TRUE)
        # [!] Variants _
        observeEvent(watch("dataLoaded"), {
            shiny::req(ScIGMA_data$mae)
            # available_variants <- rownames(ScIGMA_data$mae[["dna_variants"]])
            available_variants <- ScIGMA_data$variants.filtered$label
            shinyWidgets::updatePickerInput(
                session = session,
                inputId = "selected_variant",
                choices = available_variants
            )
        }, ignoreInit = TRUE)

        observeEvent(
            list(watch("umap_computed"),
                 watch("dna_clones_renamed"),
                 watch("compass_completed")), {
                     output$ptnUMAP_DNA_acc_dnaVariants_plot <- plotly::renderPlotly({

                         w <- Waiter$new(id = ns("ptnUMAP_DNA_acc_dnaVariants_plot"), html = spin_3(), color = transparent(0.5))
                         w$show()


                         print("ptnUMAP_DNA_acc_dnaVariants_plot : input$use_compass_variant")
                         print(input$use_compass_variant)
                         # 1. Évaluation paresseuse stricte
                         shiny::req("DNA Variants" %in% input$ptnUMAP_DNA_acc)
                         watch("umap_computed")
                         shiny::req(ScIGMA_data$seurat_object, input$selected_variant)

                         # 2. Extraction des coordonnées
                         umap_df <- as.data.frame(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings)
                         umap_df$Barcode <- rownames(umap_df)

                         # Retrieve variant id
                         tmp_selected_variant <- rownames(ScIGMA_data$variants.filtered)[ScIGMA_data$variants.filtered$label == input$selected_variant]

                         # 3. Extraction du génotype
                         print("ScIGMA_data$mae")
                         print(ScIGMA_data$mae)
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
                             text = ~paste("<b>Genotype</b>:", Variant_Genotype),
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

                         w$hide()
                         return(p)
                     })
                 }, ignoreInit = TRUE)

        # >> UMAP Unsupervised Clusters x DNA _
        # observeEvent({})

        # >> Bi-plot Gates x DNA Distribution _
        output$biplot_dna_distribution_plot <- plotly::renderPlotly({
            # Déclencheur : l'onglet doit être actif
            shiny::req(input$selected_biplot_pop, ScIGMA_data$variants.filtered)

            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            assay_colnames <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use) |> colnames()
            cel_barcode <- assay_colnames[ScIGMA_data$protein_gating_tree$gates_list[[input$selected_biplot_pop]]]

            # Extraction des données via le helper
            plot_df <- compute_population_genotype_distribution(
                mae_data = ScIGMA_data$mae,
                variant_ids = rownames(ScIGMA_data$variants.filtered),
                cell_barcodes = cel_barcode,
                use_compass = input$use_compass_biplot
            ) |> as.data.frame()

            # Palette clinique (Gris pour ADO)
            variant_colors <- c(
                "WT" = "#440154", "HET" = "#21918c",
                "HOM" = "#fde725", "Missing/ADO" = "#e0e0e0"
            )

            # Merge plot_df with variant ScIGMA_data$variants.filtered

            variants.filtered_tmp <- ScIGMA_data$variants.filtered |>
                select(-variant_id) |>
                rownames_to_column("Variant_ID")

            plot_df_joined <- left_join(plot_df,variants.filtered_tmp,
                                        by = "Variant_ID")

            plot_df_joined$Variant <- paste(plot_df_joined$protein,
                                            plot_df_joined$cdna, sep = "<br>")

            plotly::plot_ly(
                data = plot_df_joined,
                x = ~Variant,
                y = ~Percentage,
                color = ~Variant_Genotype,
                colors = variant_colors,
                type = 'bar',
                text = ~paste0(round(Percentage, 1), "%<br>(n=", Count, ")"),
                textposition = 'outside',
                textfont = list(size = 12, color = "black", family = "Arial"),
                constraintext = 'none',
                hoverinfo = 'text'
            ) |>
                plotly::layout(
                    barmode = 'group',
                    # xaxis = c(list(title = "<b>DNA Variants</b>"), prism_axis_style),
                    xaxis = c(list(title = ""), prism_axis_style),
                    yaxis = c(list(title = "<b>Frequency (%)</b>", range = c(0, 125)), prism_axis_style),
                    legend = list(title = list(text = "<b>Genotype</b>")),
                    margin = list(b = 100, t = 50)
                ) |>
                plotly::config(displaylogo = FALSE)
        })


    })
}
