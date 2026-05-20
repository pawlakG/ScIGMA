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

            # Panel 1 : Protein UMAP x SNV
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
                    title = "Protein UMAP x SNV",


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

            # Panel 2 : UMAP Unsupervised Clusters x SNV
            if (has_dna && has_clusters) {


                compass_switch_ui_clusterPtn_DNA <- shinyWidgets::materialSwitch(
                    inputId = ns("clusterPtn_DNA_compass_choice"),
                    label = if(compass_exists) {
                        "Use COMPASS imputed matrix ?"
                    } else {
                        "Run MCMC first"
                    },
                    value = compass_exists,
                    status = "success"
                )



                if (!compass_exists) {
                    compass_switch_ui_clusterPtn_DNA <- shiny::div(
                        style = "pointer-events: none; opacity: 0.5;",
                        compass_switch_ui_clusterPtn_DNA
                    )
                } else {
                    compass_switch_ui_clusterPtn_DNA <- shiny::div(
                        compass_switch_ui_clusterPtn_DNA
                    )
                }

                tabs[[length(tabs) + 1]] <- bslib::nav_panel(
                    title = "UMAP Unsupervised Clusters x SNV",
                    accordion(
                        id = ns("UMAPclusters_DNA_acc"),
                        open = FALSE,
                        # --- SOUS-PANEL 1 : CLONES ---
                        bslib::accordion_panel(
                            title = "DNA clones",
                            shiny::fluidRow(
                                grid_card(
                                    area = "main",
                                    shiny::h3("Clonal Distribution"),
                                    plotly::plotlyOutput(ns("plot_clust_clones"), height = "500px")
                                )
                            ),
                            shiny::fluidRow(
                                shiny::div(compass_switch_ui_clusterPtn_DNA, align = "left")
                            )
                        ),
                        # --- SOUS-PANEL 2 : VARIANTS ---
                        bslib::accordion_panel(
                            title = "DNA Variants",
                            shiny::fluidRow(
                                grid_card(
                                    area = "main",
                                    shiny::h3("Variant Genotypes"),
                                    plotly::plotlyOutput(ns("plot_clust_variants"),
                                                         height = "720px"
                                                         # height = "600px"
                                    )
                                )
                            )
                        )
                    )
                )
            }

            # Panel 3 : Bi-plot Gates x SNV
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
                    title = "Biplot Gates x SNV",
                    shiny::fluidRow(
                        shiny::column(3,
                                      bslib::card(
                                          shiny::h4("Population Profiling"),
                                          shinyWidgets::pickerInput(
                                              inputId = ns("selected_biplot_pop"),
                                              label = "Target Population",
                                              # choices = names(ScIGMA_data$protein_gating_tree$gates_list),
                                              choices = sapply(ScIGMA_data$protein_gating_tree$meta_list[names(ScIGMA_data$protein_gating_tree$gates_list)],
                                                               \(x) x$name)|> as.character(),
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
        # >> Protein UMAP x SNV _
        # [!] Clones
        output$ptnUMAP_DNA_acc_dnaClones_plot <- plotly::renderPlotly({
            watch("umap_computed")
            watch("dna_clones_renamed")
            watch("compass_completed")

            w <- waiter::Waiter$new(id = ns("ptnUMAP_DNA_acc_dnaClones_plot"), html = waiter::spin_3(), color = waiter::transparent(0.5))
            w$show()

            # FIX: Use isTRUE to prevent crashes if input is NULL during UI rebuild
            use_compass <- isTRUE(input$ptnUMAP_DNA_acc_gtMtx_choice)

            if (use_compass) {
                if (is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)) {
                    shiny::showNotification("COMPASS not run. Displaying raw clones.", type = "warning")
                    dna_clones_to_use <- ScIGMA_data$dna.clones_pre_compass
                } else {
                    dna_clones_to_use <- ScIGMA_data$dna.clones
                }
            } else {
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
        # [!] Variants _
        observeEvent({
            watch("dataLoaded")
            watch("dna_clones_renamed")
        }, {
            shiny::req(ScIGMA_data$mae)
            # available_variants <- rownames(ScIGMA_data$mae[["dna_variants"]])
            available_variants <- ScIGMA_data$variants.filtered$label
            shinyWidgets::updatePickerInput(
                session = session,
                inputId = "selected_variant",
                choices = available_variants
            )
        }, ignoreInit = TRUE)

        output$ptnUMAP_DNA_acc_dnaVariants_plot <- plotly::renderPlotly({
            watch("umap_computed")
            watch("dna_clones_renamed")
            watch("compass_completed")

            w <- waiter::Waiter$new(id = ns("ptnUMAP_DNA_acc_dnaVariants_plot"), html = waiter::spin_3(), color = waiter::transparent(0.5))
            w$show()

            shiny::req("DNA Variants" %in% input$ptnUMAP_DNA_acc)
            shiny::req(ScIGMA_data$seurat_object, input$selected_variant)

            # FIX: Ensure targeted variant still exists after filtering
            shiny::req(input$selected_variant %in% ScIGMA_data$variants.filtered$label)

            umap_df <- as.data.frame(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings)
            umap_df$Barcode <- rownames(umap_df)

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

            print("plot_df")
            print(head(plot_df))
            print(str(plot_df))
            print(table(plot_df$Variant_Genotype))

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

        # >> Bi-plot Gates x SNV Distribution _
        output$biplot_dna_distribution_plot <- plotly::renderPlotly({

            watch("dna_clones_renamed")
            watch("compass_completed")
            watch("clusters_computed")

            # Déclencheur : l'onglet doit être actif
            shiny::req(input$selected_biplot_pop, ScIGMA_data$variants.filtered)

            selected_biplot_gate_name <- input$selected_biplot_pop
            selected_biplot_gate_ids_list <- sapply(ScIGMA_data$protein_gating_tree$meta_list, \(x){
                x$name
            })
            selected_biplot_gate_ids <- names(selected_biplot_gate_ids_list)[selected_biplot_gate_ids_list == selected_biplot_gate_name]

            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            assay_colnames <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use) |> colnames()
            cel_barcode <- assay_colnames[ScIGMA_data$protein_gating_tree$gates_list[[selected_biplot_gate_ids]]]

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

        # ---------------------------------------------------------
        # [ NODE_ACCESS : UMAP Clusters x SNV ]
        # ---------------------------------------------------------
        # 1. Rendu du Barplot : Clusters x Clones
        output$plot_clust_clones <- plotly::renderPlotly({
            # NEW: Missing triggers added to force UI reset
            watch("dna_clones_renamed")
            watch("compass_completed")
            watch("clusters_computed")

            shiny::req(ScIGMA_data$seurat_object)

            if (isTRUE(input$clusterPtn_DNA_compass_choice)) {
                if (is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)) {
                    shiny::showNotification("COMPASS not run. Displaying raw clones.", type = "warning")
                    dna_clones_to_use <- ScIGMA_data$dna.clones_pre_compass
                } else {
                    dna_clones_to_use <- ScIGMA_data$dna.clones
                }
            } else {
                dna_clones_to_use <- ScIGMA_data$dna.clones_pre_compass
            }

            shiny::req(dna_clones_to_use)

            cell_barcodes <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], "normalized") |>
                colnames()

            # >> Setup dataframe _
            dna_clones_clusters_df <- data.frame("Barcode" =  cell_barcodes,
                                                 "Clone" = dna_clones_to_use[cell_barcodes],
                                                 "Cluster" = ScIGMA_data$seurat_object$seurat_clusters[cell_barcodes])

            dna_clones_clusters_df$Clone <- forcats::fct_na_value_to_level(dna_clones_clusters_df$Clone, level = "Missing")

            dna_clones_clusters_df_summarized <- dna_clones_clusters_df |>
                dplyr::group_by(Clone, Cluster)|>
                dplyr::summarise(Count = dplyr::n(), .groups = "drop") |>
                dplyr::group_by(Cluster) |>
                dplyr::mutate(
                    Total_In_Variant = sum(Count),
                    Percentage = (Count / Total_In_Variant) * 100
                ) |>
                dplyr::ungroup()


            plotly::plot_ly(
                data = dna_clones_clusters_df_summarized,
                x = ~Cluster,
                y = ~Percentage,
                color = ~Clone,
                # colors = variant_colors,
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


        # 2. Rendu du Barplot : Clusters x Variants (Mutations)
        output$plot_clust_variants <- plotly::renderPlotly({
            watch("dna_clones_renamed")
            watch("compass_completed")
            watch("clusters_computed")

            shiny::req(ScIGMA_data$seurat_object, ScIGMA_data$variants.filtered)

            assay_to_use <- ifelse(
                "normalized" %in% SummarizedExperiment::assayNames(
                    ScIGMA_data$mae[["proteins"]]
                ),
                "normalized",
                "counts"
            )
            assay_colnames <- SummarizedExperiment::assay(
                ScIGMA_data$mae[["proteins"]],
                assay_to_use
            ) |> colnames()

            plot_df <- compute_population_genotype_distribution(
                mae_data = ScIGMA_data$mae,
                variant_ids = rownames(ScIGMA_data$variants.filtered),
                cell_barcodes = assay_colnames,
                use_compass = isTRUE(input$clusterPtn_DNA_compass_choice),
                seurat_cluster = ScIGMA_data$seurat_object$seurat_clusters
            ) |> as.data.frame()

            variant_colors <- c(
                "WT" = "#440154",
                "HET" = "#21918c",
                "HOM" = "#fde725",
                "Missing/ADO" = "#e0e0e0"
            )

            variants_filtered_tmp <- ScIGMA_data$variants.filtered |>
                dplyr::select(-dplyr::any_of("variant_id")) |>
                tibble::rownames_to_column("Variant_ID")

            plot_df_joined <- dplyr::left_join(
                plot_df,
                variants_filtered_tmp,
                by = "Variant_ID"
            )

            id_mapping <- plot_df_joined |>
                dplyr::select(Variant_ID, gene, cdna) |>
                dplyr::distinct() |>
                dplyr::mutate(
                    cdna_clean = ifelse(is.na(cdna) | cdna == "", "Unknown", cdna),
                    Variant_Base = paste(gene, cdna_clean, sep = " - ")
                ) |>
                dplyr::group_by(Variant_Base) |>
                dplyr::mutate(n_variants = dplyr::n()) |>
                dplyr::ungroup() |>
                dplyr::mutate(
                    Variant = ifelse(
                        n_variants > 1,
                        paste0(Variant_Base, " [", Variant_ID, "]"),
                        Variant_Base
                    )
                )

            plot_df_joined <- plot_df_joined |>
                dplyr::left_join(
                    id_mapping |> dplyr::select(Variant_ID, Variant),
                    by = "Variant_ID"
                ) |>
                dplyr::arrange(as.numeric(Cluster), Variant) |>
                dplyr::mutate(
                    X_Axis_Label = paste0(
                        "<b>C", Cluster, "</b> | ", Variant
                    ),
                    X_Axis_Label = factor(
                        X_Axis_Label,
                        levels = unique(X_Axis_Label)
                    )
                ) |>
                tidyr::complete(
                    tidyr::nesting(X_Axis_Label, Cluster, Variant),
                    Variant_Genotype,
                    fill = list(Percentage = 0, Count = 0)
                ) |>
                dplyr::arrange(X_Axis_Label, Variant_Genotype)

            plotly::plot_ly(
                data = plot_df_joined,
                x = ~X_Axis_Label,
                y = ~Percentage,
                color = ~Variant_Genotype,
                colors = variant_colors,
                type = "bar",
                text = ~ifelse(
                    Percentage >= 5,
                    paste0(round(Percentage, 1), "%"),
                    ""
                ),
                textposition = "inside",
                insidetextanchor = "middle",
                textfont = list(size = 11, color = "black", family = "Arial"),
                constraintext = "none",
                hovertext = ~paste0(
                    "<b>Cluster:</b> ", Cluster, "<br>",
                    "<b>Variant:</b> ", Variant, "<br>",
                    "<b>Genotype:</b> ", Variant_Genotype, "<br>",
                    "<b>Frequency:</b> ", round(Percentage, 1),
                    "% (n=", Count, ")"
                ),
                hoverinfo = "text"
            ) |>
                plotly::layout(
                    barmode = "stack",
                    height = 700,
                    xaxis = c(
                        list(title = "", tickangle = -45, automargin = TRUE),
                        prism_axis_style
                    ),
                    yaxis = c(
                        list(title = "<b>Frequency (%)</b>", range = c(0, 105)),
                        prism_axis_style
                    ),
                    legend = list(title = list(text = "<b>Genotype</b>")),
                    margin = list(b = 140, t = 50)
                ) |>
                plotly::config(displaylogo = FALSE)
        })


    })
}
