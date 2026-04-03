#' analysis_right_Protein UI Function
#'
#' @noRd
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS actionButton checkboxInput column fluidRow h3 h4 hr
mod_analysis_Protein_ui <- function(id) {
    ns <- NS(id)
    tagList(
        # --- NEW : Réceptacle 100% R ---
        uiOutput(ns("protein_main_ui"))
    )
}

#' analysis_right_Protein Server Functions
#'
#' @noRd
#' @import ggplot2
#' @importFrom ggprism theme_prism
#' @importFrom plotly renderPlotly plot_ly layout event_data toWebGL config
#' @importFrom shiny moduleServer reactiveValues observe req updateSelectInput updateTextInput renderUI renderText actionLink div observeEvent
#' @importFrom tidyr pivot_longer
mod_analysis_Protein_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # ---------------------------------------------------------
        # NEW : Contrôleur d'affichage 100% R (renderUI)
        # ---------------------------------------------------------
        is_filtered_flag <- shiny::reactiveVal(FALSE)

        observeEvent({
            list(watch("dnaVariant_filtered"),
            watch("dataLoaded"))
        }, {
            if (!is.null(ScIGMA_data$protein.filtered) && isTRUE(ScIGMA_data$protein.filtered)) {
                is_filtered_flag(TRUE)
            } else {
                is_filtered_flag(FALSE)
            }
        }, ignoreNULL = FALSE, ignoreInit = FALSE)

        output$protein_main_ui <- shiny::renderUI({
            if (!isTRUE(is_filtered_flag())) {
                # Cas 1 : Bloqué
                card(
                    br(), br(),
                    h3("Please filter variant and cells first.",
                       style = "text-align: center; color: #7f8c8d;"),
                    br(), br()
                )
            } else {
                # Cas 2 : Autorisé (Injection de l'UI complète)
                navset_card_underline(
                    nav_panel(
                        "Description",
                        accordion(
                            id = ns("acc_ptn"),
                            open = FALSE,
                            accordion_panel(
                                "Protein ridge plot",
                                fluidRow(plotlyOutput(outputId = ns("protein_ridge")))
                            ),
                            accordion_panel(
                                "Protein bar plot",
                                fluidRow(plotlyOutput(outputId = ns("protein_bar")))
                            )
                        )
                    ),
                    nav_panel(
                        "Bi-plot",
                        fluidRow(
                            # ---- 1. Contrôles & Axes ----
                            column(
                                3,
                                grid_card(
                                    area = "sidebar",
                                    h3("Controls"),
                                    selectInput(ns("xvar"), "Axe X", choices = rownames(ScIGMA_data$mae[["proteins"]])),
                                    selectInput(ns("yvar"), "Axe Y", choices = rownames(ScIGMA_data$mae[["proteins"]])),
                                    checkboxInput(ns("logx"), "Log X", FALSE),
                                    checkboxInput(ns("logy"), "Log Y", FALSE),
                                    hr(),
                                    h4("Gating"),
                                    textInput(ns("subset_name"), "Gate name", placeholder = "ex: Clone A"),
                                    actionButton(ns("mk_subset"), "Create sub-sample", class = "btn-primary"),
                                    verbatimTextOutput(ns("selection_info"))
                                )
                            ),
                            # ---- 2. Visualisation (Bi-plot) ----
                            column(
                                6,
                                grid_card(area = "main", uiOutput(ns("biplot_container")))
                            ),
                            # ---- 3. Arborescence (Tree) ----
                            column(
                                3,
                                grid_card(
                                    area = "subsets",
                                    h3("Sub-samples"),
                                    uiOutput(ns("subsets_ui")),
                                    hr(),
                                    actionButton(ns("save_to_r6"), "Save for further analysis",
                                                 class = "btn-success", style = "margin-bottom: 10px; width: 100%;"),
                                    actionButton(ns("reset_root"), "Reset root", class = "btn-warning")
                                )
                            )
                        )
                    ),
                    nav_panel("UMAP",
                              accordion(
                                  accordion_panel(
                                      "Built UMAP model",
                                      fluidRow(
                                          column(3,
                                                 tagList(
                                                     grid_card(
                                                         area = "umap_accordion_first",
                                                         h3("UMAP parameters"),
                                                         numericInput(ns("n_neighbors"), "Neighbors", value = 15, min = 5, max = 100),
                                                         numericInput(ns("min_dist"), "Distance Min", value = 0.2, min = 0.01, max = 1.0, step = 0.1),
                                                         hr(),
                                                         actionBttn(
                                                             inputId = ns("run_umap_btn"), label = "Run UMAP",
                                                             style = "unite", color = "primary"
                                                         ),
                                                         helpText("The calculation may take a few seconds."),
                                                         hr(),
                                                         downloadButton(
                                                             outputId = ns("protein_umap_download_bttn"),
                                                             label = "Download UMAP in high resolution."
                                                         )
                                                     )
                                                 )
                                          ),
                                          column(9,
                                                 grid_card(
                                                     area = "main", h3("2D Projection"),
                                                     plotlyOutput(ns("umap_plot_build"), height = "600px")
                                                 )
                                          )
                                      )
                                  ),
                                  accordion_panel(
                                      "Markers expression",
                                      uiOutput(ns("markers_umap_panel_ui"))
                                  ),
                                  accordion_panel(
                                      "Biplot clones",
                                      uiOutput(ns("biplotClones_umap_panel_ui"))
                                  ),
                                  accordion_panel(
                                      "Unsupervised clustering",
                                      fluidRow(
                                          column(3,
                                                 tagList(
                                                     grid_card(
                                                         area = "sidebar",
                                                         h3("Clustering parameters"),
                                                         numericInput(ns("umap_clust_resolution"), "Resolution", value = 0.15, min = 0.01, max = 5),
                                                         hr(),
                                                         actionButton(ns("find_clusters_btn"), "Find clusters", class = "btn-primary", width = "100%")
                                                     )
                                                 )
                                          ),
                                          column(9,
                                                 grid_card(
                                                     area = "main", h3("2D Projection"),
                                                     plotlyOutput(ns("umap_clustering_plot"), height = "600px")
                                                 )
                                          )
                                      )
                                  )
                              )
                    )
                )
            }
        })

        # --- 1. Initialization & State (Bi-plot Logic) ---

        # State: Stores indices (pointers) only. Zero data duplication.
        r_state <- reactiveValues(
            subsets = list(),
            subset_meta = list(),
            current_view = "root",
            temp_selection = NULL
        )

        # NEW : Déclaré ici pour que la purge puisse le réinitialiser
        bound_listeners <- c()

        # Initialize Inputs & Root Gate from R6 Object
        observeEvent(
            list(watch("dataLoaded"),
                 watch("dnaVariant_filtered"),
                 input$reset_root
                 ),
        {
            req(ScIGMA_data$mae)

            # Sécurité : On s'assure que les protéines sont bien dans le dataset
            if (!("proteins" %in% names(ScIGMA_data$mae))) return()

            assay_to_use_state <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")
            n_cells_current <- ncol(SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use_state))

            # FIX CRITIQUE : Purge inconditionnelle (Déterministe)
            r_state$subsets <- list()
            r_state$subsets[["root"]] <- 1:n_cells_current

            r_state$subset_meta <- list()
            r_state$subset_meta[["root"]] <- list(
                name = "Tout",
                parent = NA,
                depth = 0
            )

            r_state$current_view <- "root"
            r_state$temp_selection <- NULL

            # Atomisation des écouteurs du menu latéral
            bound_listeners <<- c()

            updateSelectInput(
                session = session,
                inputId = "protein_umap_markers",
                choices = c("None", rownames(ScIGMA_data$mae[["proteins"]])),
                selected = "None"
            )
        }, ignoreInit = FALSE)

        # --- 2. Bi-plot Rendering (Optimized) ---

        # Container output to ensure correct sizing
        output$biplot_container <- renderUI({
            plotlyOutput(ns("biplot"), height = "600px")
        })

        output$biplot <- renderPlotly({
            req(input$xvar, input$yvar, r_state$current_view)

            current_indices <- r_state$subsets[[r_state$current_view]]

            # UPDATED : Extraction directe depuis le MAE (Source de vérité unique)
            # On cherche les données normalized en priorité, sinon les counts bruts.
            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            # Matrice native (Protéines x Cellules) -> On extrait la ligne spécifique (ptn) pour les cellules ciblées
            raw_x <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[input$xvar, current_indices]
            raw_y <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[input$yvar, current_indices]

            plot_df <- data.frame(
                x = raw_x,
                y = raw_y,
                custom_id = current_indices
            )

            if (input$logx) plot_df$x <- log1p(plot_df$x)
            if (input$logy) plot_df$y <- log1p(plot_df$y)


            plot_ly(
                data = plot_df,
                x = ~x,
                y = ~y,
                key = ~custom_id,
                type = "scatter",
                mode = "markers",
                source = ns("gating_plot"),
                marker = list(size = 5, opacity = 1, color = "#2c3e50")
            ) %>%
                layout(
                    title = paste("Gate:", r_state$subset_meta[[r_state$current_view]]$name),
                    xaxis = list(title = input$xvar, range = list(0, max(plot_df$x)+1)),
                    yaxis = list(title = input$yvar, range = list(0, max(plot_df$y)+1)),
                    dragmode = "lasso"
                ) %>%
                # toWebGL() %>%
                config(displaylogo = FALSE)
        })

        # --- 3. Selection & Gating Logic ---

        # Listener
        observeEvent(event_data("plotly_selected", source = ns("gating_plot")), {
            sel <- event_data("plotly_selected", source = ns("gating_plot"))

            # 1. Le Bouclier Absolu : Intercepte NULL, list(), et dataframes vides
            if (is.null(sel) || length(sel) == 0 || (is.data.frame(sel) && nrow(sel) == 0)) {
                r_state$temp_selection <- NULL
                return()
            }

            # 2. Cascade d'extraction
            extracted_indices <- NULL
            if ("customdata" %in% colnames(sel)) {
                extracted_indices <- as.numeric(sel$customdata)
            } else if ("key" %in% colnames(sel)) {
                extracted_indices <- as.numeric(sel$key)
            } else if ("pointNumber" %in% colnames(sel)) {
                current_indices <- r_state$subsets[[r_state$current_view]]
                extracted_indices <- current_indices[as.numeric(sel$pointNumber) + 1]
            }

            # 3. Purge des "Fantômes" (NA générés par le lasso sur le fond du plot)
            if (!is.null(extracted_indices)) {
                valid_indices <- extracted_indices[!is.na(extracted_indices)]

                if (length(valid_indices) > 0) {
                    r_state$temp_selection <- valid_indices
                } else {
                    r_state$temp_selection <- NULL
                }
            } else {
                r_state$temp_selection <- NULL
            }
        }, ignoreNULL = FALSE, ignoreInit = TRUE)

        # Feedback Info
        output$selection_info <- renderText({
            req(r_state$temp_selection)
            n <- length(r_state$temp_selection)
            parent_n <- length(r_state$subsets[[r_state$current_view]])
            pct <- round(n / parent_n * 100, 1)
            paste0("Sélection : ", n, " (", pct, "%)")
        })

        # Create Subset
        observeEvent(input$mk_subset, {
            req(r_state$temp_selection)

            new_id <- paste0("sub_", as.numeric(Sys.time()))
            name_val <- input$subset_name
            if (name_val == "") name_val <- paste("Gate", length(r_state$subsets))

            # Save Indices & Meta
            r_state$subsets[[new_id]] <- r_state$temp_selection

            parent_id <- r_state$current_view
            r_state$subset_meta[[new_id]] <- list(
                name = name_val,
                parent = parent_id,
                depth = r_state$subset_meta[[parent_id]]$depth + 1
            )

            # Update View
            r_state$current_view <- new_id
            r_state$temp_selection <- NULL
            updateTextInput(session, "subset_name", value = "")
        })

        # --- 4. Navigation Tree ---

        output$subsets_ui <- renderUI({
            req(length(r_state$subset_meta) > 0)
            meta_list <- r_state$subset_meta

            # FIX CRITIQUE : Moteur de rendu récursif (Depth-First Search)
            # Permet de regrouper visuellement les enfants sous leur VRAI parent
            build_tree_ui <- function(node_id) {
                meta <- meta_list[[node_id]]
                indent <- meta$depth * 15

                # UX : Met en surbrillance bleue le noeud actuellement sélectionné
                style_str <- paste0(
                    "margin-left:", indent, "px; ",
                    "display: block; padding: 2px 0; ",
                    if(node_id == r_state$current_view) "font-weight:bold; color:#007bff;" else "color:#333;"
                )

                # 1. Rendu du noeud parent
                current_ui <- div(
                    actionLink(ns(paste0("go_", node_id)), label = paste("•", meta$name), style = style_str)
                )

                # 2. Identification stricte de ses enfants
                children_ids <- names(meta_list)[sapply(meta_list, function(m) {
                    !is.na(m$parent) && m$parent == node_id
                })]

                # 3. Plongée récursive
                if (length(children_ids) > 0) {
                    children_ui <- lapply(children_ids, build_tree_ui)
                    return(tagList(current_ui, do.call(tagList, children_ui)))
                } else {
                    return(current_ui)
                }
            }

            # Lancement de la construction à partir de la racine
            if ("root" %in% names(meta_list)) {
                return(build_tree_ui("root"))
            } else {
                return(NULL)
            }
        })

        # Tree Navigation Logic
        bound_listeners <- c()

        observe({
            sids <- names(r_state$subsets)

            # On n'isole que les nouveaux clones qui n'ont pas encore d'écouteur
            new_sids <- setdiff(sids, bound_listeners)

            lapply(new_sids, function(sid) {
                # Suppression de once = TRUE. Le clic est permanent.
                observeEvent(input[[paste0("go_", sid)]], {
                    r_state$current_view <- sid
                    r_state$temp_selection <- NULL
                }, ignoreInit = TRUE)
            })

            # On met à jour le registre local à la session (<<- pointe vers l'env du module)
            if (length(new_sids) > 0) {
                bound_listeners <<- c(bound_listeners, new_sids)
            }
        })


        # --- 5. ORIGINAL PLOTS (UNCHANGED) ---
        # These blocks are kept exactly as in the original file provided.

        # Render ridge plot for all proteins
        output$protein_ridge <- renderPlotly({
            watch("dnaVariant_filtered")
            render_protein_ridge_plot(obj = ScIGMA_data)
        })

        # --- 6. PERSISTENCE LOGIC (R6 Bridge) ---

        observeEvent(input$save_to_r6, {
            req(ScIGMA_data)

            # Vérifier qu'il y a quelque chose à sauver
            if (length(r_state$subsets) <= 1) {
                showNotification("Aucune sous-population créée à enregistrer.", type = "warning")
                return()
            }

            # FIX CRITIQUE : Alignement strict sur le nom R6 (protein_gating_tree et non protein.gating_tree)
            # Cela permet à la méthode ScIGMA_data$reset_analysis() de l'effacer proprement
            ScIGMA_data$protein_gating_tree <- list("gates_list" = r_state$subsets,
                                                    "meta_list" = r_state$subset_meta)

            showNotification(
                paste("Success !", length(r_state$subsets), "populations saved."),
                type = "message"
            )
            print("saved  ScIGMA_data$protein_gating_tree")
            print( ScIGMA_data$protein_gating_tree)
        })

        # --- 5. UMAP Calculation ---

        # Trigger global via bouton (recommandé pour éviter les calculs intempestifs)
        observeEvent(input$run_umap_btn, {
            gargoyle::trigger("launch_umap")
        })

        # Bloc de calcul réactif
        observe({
            # On écoute explicitement le trigger ET les paramètres
            watch("launch_umap")

            # show_modal_spinner(text = "Filtering and annotating DNA variants...")
            w <- Waiter$new(
                # id = "umap_clustering_plot", # Peut cibler un plot précis ou toute la page
                html = spin_loaders(1, color = "black"),             # Spinner discret et élégant
                color = "rgba(255, 255, 255, 0.5)" # Fond blanc semi-transparent (pas de bloc blanc opaque)
            )
            w$show()

            # Récupération isolée des paramètres pour ne pas relancer le calcul
            # si l'utilisateur change juste le chiffre sans cliquer
            n_neighbors <- isolate(input$n_neighbors)
            min_dist <- isolate(input$min_dist)

            # req(ScIGMA_data$protein.mtx)
            req(ScIGMA_data$mae[["proteins"]])

            message("Computing UMAP...")

            # Sélection des données
            # data_to_use <- if (!is.null(ScIGMA_data$protein.mtx.filtered.normalized)) {
            #     ScIGMA_data$protein.mtx.filtered.normalized
            # } else {
            #     ScIGMA_data$protein.mtx
            # }

            # Compute UMAP


            ScIGMA_data$seurat_object <- RunUMAP(ScIGMA_data$seurat_object,
                                                 dims = 1:(nrow(ScIGMA_data$seurat_object)-2),
                                                 min.dist = min_dist,
                                                 n.neighbors = n_neighbors)


            # umap_coords <- run_umap_protein(
            #     expression_matrix = data_to_use,
            #     n_neighbors = n_neighbors,
            #     min_dist = min_dist,
            #     n_components = 2,
            #     metric = "cosine"
            # )

            # Stockage dans l'objet R6
            # ScIGMA_data$protein.umap <- umap_coords

            w$hide()
            gargoyle::trigger("umap_computed")
            message("Calcul UMAP terminé.")

        }) %>%
            bindEvent(input$run_umap_btn) # Force l'exécution UNIQUEMENT sur clic



        # render UMAP without clustering
        observeEvent({watch("umap_computed")},{
            umap_cluster <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                as.data.frame()

            ScIGMA_data$umaps$umap_protein_general <- umap_cluster |>
                ggplot(aes(x=umap_1, y=umap_2)) +
                geom_point(size=2) +
                xlab('UMAP 1') +
                ylab('UMAP 2') +
                theme_prism()

            # 2. Rendering
            output$umap_plot_build <- renderPlotly({

                p <- plot_ly(data = umap_cluster,
                             x = ~umap_1,
                             y = ~umap_2,
                             type = 'scatter',
                             mode = 'markers',
                             marker = list(size = 6,  opacity = 1)) %>%
                    layout(xaxis = list(title = "UMAP 1"),
                           yaxis = list(title = "UMAP 2")) %>%
                    toWebGL()

                return(p)
            })
        },ignoreInit = TRUE)



        # Render Barplot for all proteins
        output$protein_bar <- renderPlotly({
            watch("dnaVariant_filtered")
            # ScIGMA_data <- normalizeProtein(ScIGMA_data)
            plot_protein_barplot(obj = ScIGMA_data)
        })

        # render UMAP with clustering
        observeEvent({watch("umap_computed")
            input$find_clusters_btn},{

                clustering_resolution <- input$umap_clust_resolution
                # Compute umap
                if (is.null(ScIGMA_data$seurat_object@reductions$umap)){
                    message("UMAP is missing, computing default UMAP") # TODO: display notification

                    ScIGMA_data$seurat_object <- RunUMAP(ScIGMA_data$seurat_object,
                                                         dims = 1:(nrow(ScIGMA_data$seurat_object)-2),
                                                         min.dist = 0.15,
                                                         n.neighbors = 30,
                                                         future.seed=TRUE)
                }
                # 1. Clustering
                ScIGMA_data$seurat_object <- FindClusters(ScIGMA_data$seurat_object,
                                                          resolution = clustering_resolution)
                umap_cluster <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                    as.data.frame()
                umap_cluster$cluster = ScIGMA_data$seurat_object$seurat_clusters[rownames(umap_cluster)]

                # 2. Rendering
                # output$umap_clustering_plot <- renderPlotly({plot_ly(data = umap_cluster,
                #                                                      x = ~umap_1,
                #                                                      y = ~umap_2,
                #                                                      type = 'scatter',
                #                                                      mode = 'markers',
                #                                                      color = ~cluster,
                #                                                      marker = list(size = 6,  opacity = 1)) %>%
                #         layout(xaxis = list(title = "UMAP 1"),
                #                yaxis = list(title = "UMAP 2")) %>%
                #         toWebGL() # Toujours optimiser pour le single-cell
                # })
                output$umap_clustering_plot <- renderPlotly({

                    w <- Waiter$new(
                        id = "umap_plot_build",
                        html = spin_3(),
                        color = transparent(0.5)
                    )

                    w$show()

                    p <- plot_ly(data = umap_cluster,
                                 x = ~umap_1,
                                 y = ~umap_2,
                                 type = 'scatter',
                                 mode = 'markers',
                                 color = ~cluster,
                                 marker = list(size = 6, opacity = 0.9)) %>%
                        layout(
                            # Configuration de l'axe X
                            xaxis = list(
                                title = list(text = "<b>UMAP 1</b>", font = list(size = 16, color = "black")),
                                tickfont = list(size = 14, color = "black", family = "Arial"),
                                showline = TRUE,
                                linewidth = 2,        # Épaisseur de l'axe (style Prism)
                                linecolor = "black",
                                mirror = FALSE,       # Garde l'axe uniquement en bas
                                ticks = "outside",    # Ticks vers l'extérieur
                                tickwidth = 2,
                                gridcolor = 'transparent',
                                zeroline = FALSE
                            ),
                            # Configuration de l'axe Y
                            yaxis = list(
                                title = list(text = "<b>UMAP 2</b>", font = list(size = 16, color = "black")),
                                tickfont = list(size = 14, color = "black", family = "Arial"),
                                showline = TRUE,
                                linewidth = 2,        # Épaisseur de l'axe
                                linecolor = "black",
                                mirror = FALSE,
                                ticks = "outside",
                                tickwidth = 2,
                                gridcolor = 'transparent',
                                zeroline = FALSE
                            ),
                            # Paramètres généraux
                            plot_bgcolor = "white",  # Fond blanc
                            paper_bgcolor = "white",
                            legend = list(
                                font = list(size = 12, family = "Arial"),
                                title = list(text = "<b>Cluster</b>")
                            ),
                            margin = list(l = 50, r = 50, b = 50, t = 50) # Marges propres
                        ) %>%
                        toWebGL()
                    w$hide()
                    return(p)
                })

            },ignoreInit = TRUE)



        # [ NODE_ACCESS : Markers projection ]
        # ----------------------------------------------------- _
        observeEvent(watch("umap_computed"),
                     {
                         req(ScIGMA_data$seurat_object)

                         if(is.null(ScIGMA_data$seurat_object@reductions$umap)){
                             output$markers_umap_panel_ui <-  renderUI({
                                 tagList(
                                     fluidRow(
                                         h2("Please compute UMAP first")
                                     )
                                 )
                             })
                         } else {
                             # >> Add protein marker information to seurat object _
                             plot_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                                 as.data.frame()


                             assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

                             # Matrice native (Protéines x Cellules) -> On extrait la ligne spécifique (ptn) pour les cellules ciblées
                             protein_markers_df <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use) |> t()

                             print("assay_to_use")
                             print(assay_to_use)
                             print("protein_markers_df")
                             print(dim(protein_markers_df))
                             print("rownames(ScIGMA_data$seurat_object@meta.data)")
                             print(length(rownames(ScIGMA_data$seurat_object@meta.data)))


                             ScIGMA_data$seurat_object@meta.data <- cbind(
                                 ScIGMA_data$seurat_object@meta.data[intersect(rownames(protein_markers_df),
                                                                               rownames(ScIGMA_data$seurat_object@meta.data)),],
                                 # ScIGMA_data$seurat_object@meta.data,
                                 # protein_markers_df[rownames(ScIGMA_data$seurat_object@meta.data),])
                                 protein_markers_df[intersect(rownames(protein_markers_df),
                                                              rownames(ScIGMA_data$seurat_object@meta.data)),])

                             output$markers_umap_panel_ui <-  renderUI({
                                 fluidRow(
                                     tagList(
                                         column(3,
                                                grid_card(
                                                    area = "umap_marker_select",
                                                    virtualSelectInput(
                                                        inputId = ns("protein_umap_markers"),
                                                        label = "Markers :",
                                                        choices = colnames(protein_markers_df),
                                                        multiple = TRUE,
                                                        width = "100%",
                                                        dropboxWrapper = "body",
                                                        selected = "None"
                                                    ),
                                                    hr(),
                                                    actionBttn(
                                                        inputId = ns("protein_umap_markers_bttn"),
                                                        label = "Plot markers projections",
                                                        style = "unite",
                                                        color = "primary"
                                                    ),
                                                    hr(),
                                                    downloadButton(
                                                        outputId = ns("protein_umap_markers_download_bttn"),
                                                        label = "Download UMAPs"
                                                    )

                                                )
                                         ),
                                         column(9,
                                                plotOutput(ns("umap_plot_markers"), height = "600px"))
                                     )
                                 )
                             })
                         }

                     })

        # [ NODE_ACCESS : Biplot clones projection ]
        # ----------------------------------------------------- _
        observeEvent(watch("umap_computed"),
                     {
                         req(ScIGMA_data$seurat_object)

                         if(is.null(ScIGMA_data$seurat_object@reductions$umap)){
                             output$biplotClones_umap_panel_ui <-  renderUI({
                                 tagList(
                                     fluidRow(
                                         h2("Please compute UMAP first")
                                     )
                                 )
                             })
                         } else {
                             # >> Add protein marker information to seurat object _
                             plot_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                                 as.data.frame()


                             assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

                             # Matrice native (Protéines x Cellules) -> On extrait la ligne spécifique (ptn) pour les cellules ciblées
                             protein_markers_df <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use) |> t()


                             ScIGMA_data$seurat_object@meta.data <- cbind(
                                 ScIGMA_data$seurat_object@meta.data[intersect(rownames(protein_markers_df),
                                                                               rownames(ScIGMA_data$seurat_object@meta.data)),],
                                 # ScIGMA_data$seurat_object@meta.data,
                                 # protein_markers_df[rownames(ScIGMA_data$seurat_object@meta.data),])
                                 protein_markers_df[intersect(rownames(protein_markers_df),
                                                              rownames(ScIGMA_data$seurat_object@meta.data)),])

                             output$biplotClones_umap_panel_ui <-  renderUI({
                                 fluidRow(
                                     tagList(
                                         column(3,
                                                grid_card(
                                                    area = "umap_marker_select",
                                                    virtualSelectInput(
                                                        inputId = ns("protein_umap_markers"),
                                                        label = "Markers :",
                                                        choices = colnames(protein_markers_df),
                                                        multiple = TRUE,
                                                        width = "100%",
                                                        dropboxWrapper = "body",
                                                        selected = "None"
                                                    ),
                                                    hr(),
                                                    actionBttn(
                                                        inputId = ns("protein_umap_markers_bttn"),
                                                        label = "Plot markers projections",
                                                        style = "unite",
                                                        color = "primary"
                                                    ),
                                                    hr(),
                                                    downloadButton(
                                                        outputId = ns("protein_umap_markers_download_bttn"),
                                                        label = "Download UMAPs"
                                                    )

                                                )
                                         ),
                                         column(9,
                                                plotOutput(ns("umap_plot_markers"), height = "600px"))
                                     )
                                 )
                             })
                         }

                     })



        observeEvent(input$protein_umap_markers_bttn,{

            # >> handle inputs _

            req(input$protein_umap_markers)
            umap_marker_choice <- input$protein_umap_markers

            # >> Handle markers selection _

            if (is.null(umap_marker_choice)){
                umap_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                    as.data.frame()
                umap_df$ptn_expression <- 0
                umap_df$marker <- "No marker selected"
            } else if (length(umap_marker_choice) == 1) {
                umap_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                    as.data.frame()
                umap_df$ptn_expression <- ScIGMA_data$seurat_object@meta.data[,umap_marker_choice]
                umap_df$marker <- umap_marker_choice
            } else {
                umap_df <- cbind(ScIGMA_data$seurat_object@reductions$umap@cell.embeddings,
                                 ScIGMA_data$seurat_object@meta.data[,umap_marker_choice]) |>
                    as.data.frame()
                umap_df <- umap_df |>
                    pivot_longer(-c(umap_1, umap_2),
                                 values_to = "ptn_expression",
                                 names_to = "marker")
            }

            ScIGMA_data$umaps$umap_protein_markers <- umap_df |>
                ggplot(aes(x=umap_1, y=umap_2, color=ptn_expression)) +
                geom_point(size = 1) +
                facet_wrap(~marker) +
                scale_color_viridis_c("inferno") +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                ggprism::theme_prism()

            output$umap_plot_markers <- renderPlot(ScIGMA_data$umaps$umap_protein_markers )
        },ignoreInit = TRUE)

        # [ NODE_ACCESS : Download primary UMAP ]
        # ----------------------------------------------------- _
        output$protein_umap_download_bttn <- downloadHandler(
            filename = function() {
                "protein_umap.png"
            },
            content = function(file) {
                ggsave(file, ScIGMA_data$umaps$umap_protein_general,
                       units = "in",
                       width = 12, height = 10, dpi = 300,
                       bg = "transparent")
            }
        )

        # [ NODE_ACCESS : Download markers UMAPs ]
        # ----------------------------------------------------- _
        output$protein_umap_markers_download_bttn <- downloadHandler(
            filename = function() {
                "protein_umap_markers.png"
            },
            content = function(file) {
                ggsave(file, ScIGMA_data$umaps$umap_protein_markers ,
                       units = "in",
                       width = 12, height = 10, dpi = 300,
                       bg = "transparent")
            }
        )

    })
}

## To be copied in the UI
# mod_analysis_Protein_ui("analysis_Protein_1")

## To be copied in the server
# mod_analysis_Protein_server("analysis_Protein_1")
