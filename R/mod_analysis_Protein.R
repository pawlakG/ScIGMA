#' analysis_right_Protein UI Function
#'
#' @noRd
#' @importFrom shiny NS tagList navset_card_underline nav_panel accordion accordion_panel fluidRow column grid_card h3 h4 hr selectInput checkboxInput actionButton uiOutput textInput verbatimTextOutput
#' @importFrom plotly plotlyOutput
mod_analysis_Protein_ui <- function(id) {
    ns <- NS(id)
    tagList(
        navset_card_underline(
            nav_panel(
                "Description",
                accordion(
                    id = ns("acc_ptn"),
                    open = FALSE,
                    accordion_panel(
                        "Protein ridge plot",
                        fluidRow(
                            plotlyOutput(outputId = ns("protein_ridge"))
                        )
                    ),
                    accordion_panel(
                        "Protein bar plot",
                        fluidRow(
                            plotlyOutput(outputId = ns("protein_bar"))
                        )
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
                            # Sélection des marqueurs (rempli côté serveur via R6)
                            selectInput(ns("xvar"), "Axe X", choices = NULL),
                            selectInput(ns("yvar"), "Axe Y", choices = NULL),
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
                        grid_card(
                            area = "main",
                            uiOutput(ns("biplot_container"))
                        )
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
                                         class = "btn-success",
                                         style = "margin-bottom: 10px; width: 100%;"),
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
                                  # Colonne de Contrôle
                                  column(3,
                                         tagList(
                                             grid_card(
                                                 area = "umap_accordion_first",
                                                 h3("UMAP parameters"),

                                                 # Paramètres de la fonction standalone
                                                 numericInput(ns("n_neighbors"), "Neighbors", value = 15, min = 5, max = 100),
                                                 numericInput(ns("min_dist"), "Distance Min", value = 0.2, min = 0.01, max = 1.0, step = 0.1),

                                                 hr(),
                                                 # actionButton(ns("run_umap_btn"), "Run UMAP", class = "btn-primary", width = "100%"),

                                                 actionBttn(
                                                     inputId = ns("run_umap_btn"),
                                                     label = "Run UMAP",
                                                     style = "unite",
                                                     color = "primary"
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

                                  # Colonne de Visualisation
                                  column(9,
                                         grid_card(
                                             area = "main",
                                             h3("2D Projection"),

                                             # shinycssloaders::withSpinner(
                                             plotlyOutput(ns("umap_plot_build"), height = "600px")
                                             # type = 6, # Un style de loader (cercle qui tourne)
                                             # color = "#007bff"
                                             # )
                                         )
                                  )
                              )
                          ),
                          accordion_panel(
                              "Markers expression",
                              uiOutput(ns("markers_umap_panel_ui"))
                          ),
                          accordion_panel(
                              "Unsupervised clustering",

                              fluidRow(
                                  # Colonne de Contrôle
                                  column(3,
                                         tagList(
                                             grid_card(
                                                 area = "sidebar",
                                                 h3("Clustering parameters"),

                                                 # Paramètres de la fonction standalone
                                                 numericInput(ns("umap_clust_resolution"), "Resolution", value = 0.15, min = 0.01, max = 5),
                                                 hr(),
                                                 actionButton(ns("find_clusters_btn"), "Find clusters", class = "btn-primary", width = "100%")
                                             )
                                         )
                                  ),
                                  # Colonne de Visualisation
                                  column(9,
                                         grid_card(
                                             area = "main",
                                             h3("2D Projection"),
                                             plotlyOutput(ns("umap_clustering_plot"), height = "600px")
                                         )
                                  )
                              )
                          )
                      )
            )
        )
    )
}

#' analysis_right_Protein Server Functions
#'
#' @noRd
#' @importFrom ggprism theme_prism
#' @importFrom plotly renderPlotly plot_ly layout event_data toWebGL config
#' @importFrom shiny moduleServer reactiveValues observe req updateSelectInput updateTextInput renderUI renderText actionLink div observeEvent
mod_analysis_Protein_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # --- 1. Initialization & State (Bi-plot Logic) ---

        # State: Stores indices (pointers) only. Zero data duplication.
        r_state <- reactiveValues(
            subsets = list(),
            subset_meta = list(),
            current_view = "root",
            temp_selection = NULL
        )

        # Initialize Inputs & Root Gate from R6 Object
        observeEvent({
            watch("dnaVariant_filtered")
        },{

            req(ScIGMA_data$protein.mtx)

            # Populate SelectInputs for Markers
            ptn_names <- colnames(ScIGMA_data$protein.mtx)
            updateSelectInput(session, "xvar", choices = ptn_names, selected = ptn_names[1])
            updateSelectInput(session, "yvar", choices = ptn_names, selected = ptn_names[2])

            # Initialize Root (All Cells)
            if (is.null(r_state$subsets[["root"]])) {
                n_cells <- nrow(ScIGMA_data$protein.mtx)
                r_state$subsets[["root"]] <- 1:n_cells
                r_state$subset_meta[["root"]] <- list(
                    name = "Tout",
                    parent = NA,
                    depth = 0
                )
            }
            updateSelectInput(
                session = session,      # Indispensable dans un module
                inputId = "protein_umap_markers",       # Pas besoin de ns() ici, c'est géré par session
                choices = c("None", colnames(ScIGMA_data$protein.mtx)),
                selected = "None" # Optionnel : définir la sélection par défaut
            )
        })

        # --- 2. Bi-plot Rendering (Optimized) ---

        # Container output to ensure correct sizing
        output$biplot_container <- renderUI({
            plotlyOutput(ns("biplot"), height = "600px")
        })

        output$biplot <- renderPlotly({
            req(input$xvar, input$yvar, r_state$current_view)

            # 1. Get indices (Recursive filtering)
            current_indices <- r_state$subsets[[r_state$current_view]]

            # 2. Fetch data from R6 (Reference access)
            # We assume ScIGMA_data$protein.mtx is the matrix.
            # raw_x <- ScIGMA_data$protein.mtx[current_indices, input$xvar]
            # raw_y <- ScIGMA_data$protein.mtx[current_indices, input$yvar]
            raw_x <- t(ScIGMA_data$seurat_object@assays$RNA$data)[current_indices, input$xvar]
            raw_y <- t(ScIGMA_data$seurat_object@assays$RNA$data)[current_indices, input$yvar]

            plot_df <- data.frame(
                x = raw_x,
                y = raw_y,
                custom_id = current_indices # Key for mapping back to global R6 rows
            )

            # 3. Transform
            if (input$logx) plot_df$x <- log1p(plot_df$x)
            if (input$logy) plot_df$y <- log1p(plot_df$y)

            # 4. Render (WebGL for performance)
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
                    xaxis = list(title = input$xvar,
                                 range = list(0, max(plot_df$x)+1)),
                    yaxis = list(title = input$yvar,
                                 range = list(0, max(plot_df$y)+1)),
                    dragmode = "lasso"
                ) %>%
                toWebGL() %>%
                config(displaylogo = FALSE)
        })

        # --- 3. Selection & Gating Logic ---

        # Listener
        observeEvent(event_data("plotly_selected", source = ns("gating_plot")), {
            sel <- event_data("plotly_selected", source = ns("gating_plot"))
            if (is.null(sel) || nrow(sel) == 0) {
                r_state$temp_selection <- NULL
            } else {
                r_state$temp_selection <- as.numeric(sel$key)
            }
        })

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
            ids <- names(r_state$subsets)
            ui_elems <- lapply(ids, function(sid) {
                meta <- r_state$subset_meta[[sid]]
                indent <- meta$depth * 15
                style_str <- paste0("margin-left:", indent, "px;", if(sid == r_state$current_view) "font-weight:bold;" else "")

                div(
                    actionLink(ns(paste0("go_", sid)), label = paste("•", meta$name), style = style_str)
                )
            })
            do.call(tagList, ui_elems)
        })

        # Tree Navigation Logic
        observe({
            lapply(names(r_state$subsets), function(sid) {
                observeEvent(input[[paste0("go_", sid)]], {
                    r_state$current_view <- sid
                    r_state$temp_selection <- NULL
                }, ignoreInit = TRUE, once = TRUE)
            })
        })

        observeEvent(input$reset_root, {
            r_state$current_view <- "root"
            r_state$temp_selection <- NULL
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
            if (length(r_state$subsets) <= 1) { # <= 1 car "root" existe toujours
                showNotification("Aucune sous-population créée à enregistrer.", type = "warning")
                return()
            }

            # Appel de la méthode R6 définie à l'étape 1
            # # On passe les indices bruts et les métadonnées
            # ScIGMA_data$save_gating_tree(
            #     gates_list = r_state$subsets,
            #     meta_list = r_state$subset_meta
            # )
            ScIGMA_data$protein.gating_tree <- list("gates_list" = r_state$subsets,
                                                    "meta_list" = r_state$subset_meta)

            showNotification(
                paste("Success !", length(r_state$subsets), "populations saved."),
                type = "message"
            )
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

            # Récupération isolée des paramètres pour ne pas relancer le calcul
            # si l'utilisateur change juste le chiffre sans cliquer
            n_neighbors <- isolate(input$n_neighbors)
            min_dist <- isolate(input$min_dist)

            req(ScIGMA_data$protein.mtx)

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

            gargoyle::trigger("umap_computed")
            message("Calcul UMAP terminé.")

        }) %>%
            bindEvent(input$run_umap_btn) # Force l'exécution UNIQUEMENT sur clic



        # render UMAP without clustering
        observeEvent({watch("umap_computed")},{
            umap_cluster <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                as.data.frame()

            print(head(umap_cluster))

            ScIGMA_data$umaps$umap_protein_general <- umap_cluster |>
                ggplot(aes(x=umap_1, y=umap_2)) +
                geom_point(size=2) +
                xlab('UMAP 1') +
                ylab('UMAP 2') +
                theme_prism()

            print("ScIGMA_data$umaps$umap_protein_general")
            print(ScIGMA_data$umaps$umap_protein_general)

            # 2. Rendering
            output$umap_plot_build <- renderPlotly({plot_ly(data = umap_cluster,
                                                            x = ~umap_1,
                                                            y = ~umap_2,
                                                            type = 'scatter',
                                                            mode = 'markers',
                                                            marker = list(size = 6,  opacity = 1)) %>%
                    layout(xaxis = list(title = "UMAP 1"),
                           yaxis = list(title = "UMAP 2")) %>%
                    toWebGL() # Toujours optimiser pour le single-cell
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
                output$umap_clustering_plot <- renderPlotly({plot_ly(data = umap_cluster,
                                                                     x = ~umap_1,
                                                                     y = ~umap_2,
                                                                     type = 'scatter',
                                                                     mode = 'markers',
                                                                     color = ~cluster,
                                                                     marker = list(size = 6,  opacity = 1)) %>%
                        layout(xaxis = list(title = "UMAP 1"),
                               yaxis = list(title = "UMAP 2")) %>%
                        toWebGL() # Toujours optimiser pour le single-cell
                })
            },ignoreInit = TRUE)



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
                             protein_markers_df <- if(!is.null(ScIGMA_data$protein.mtx.filtered.normalized)){
                                 ScIGMA_data$protein.mtx.filtered.normalized
                             } else {
                                 ScIGMA_data$protein.mtx
                             }
                             ScIGMA_data$seurat_object@meta.data <- cbind(ScIGMA_data$seurat_object@meta.data,
                                                                          protein_markers_df[rownames(ScIGMA_data$seurat_object@meta.data),])

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


        observeEvent(input$protein_umap_markers_bttn,{

            # >> handle inputs _

            req(input$protein_umap_markers)
            umap_marker_choice <- input$protein_umap_markers

            message("umap_marker_choice")
            print(umap_marker_choice)

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
