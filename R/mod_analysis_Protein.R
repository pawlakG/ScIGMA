#' analysis_right_Protein UI Function
#'
#' @noRd
#' @importFrom plotly plotlyOutput
#' @importFrom shiny NS actionButton checkboxInput column fluidRow h3 h4 hr
mod_analysis_Protein_ui <- function(id) {
    ns <- NS(id)
    tagList(
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

        # [ NODE_ACCESS : Contrôleur d'affichage 100% R (renderUI)
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
            watch("dnaVariant_filtered")
            watch("compass_completed")
            watch("umap_computed") # NEW : Déclenche l'UI quand l'UMAP est prêt

            if (!isTRUE(is_filtered_flag())) {
                # Cas 1 : Bloqué
                card(
                    br(), br(),
                    h3("Please filter variant and cells first.",
                       style = "text-align: center; color: #7f8c8d;"),
                    br(), br()
                )
            } else {
                # NEW : Vérification de l'existence du modèle UMAP
                umap_ready <- !is.null(ScIGMA_data$seurat_object@reductions$umap)

                # NEW : Persistance des inputs utilisateurs lors des rechargements UI
                curr_features <- shiny::isolate(input$umap_features)
                if (is.null(curr_features) && !is.null(ScIGMA_data$mae[["proteins"]])) {
                    curr_features <- rownames(ScIGMA_data$mae[["proteins"]])
                }

                curr_n_neighbors <- shiny::isolate(input$n_neighbors)
                if (is.null(curr_n_neighbors)) curr_n_neighbors <- 15

                curr_min_dist <- shiny::isolate(input$min_dist)
                if (is.null(curr_min_dist)) curr_min_dist <- 0.2

                # NEW : Capture et persistance de l'onglet principal actif
                curr_tab <- shiny::isolate(input$protein_tabs)
                if (is.null(curr_tab)) curr_tab <- "tab_description"

                # Cas 2 : Autorisé (Injection de l'UI complète)
                navset_card_underline(
                    id = ns("protein_tabs"),     # NEW : Identifiant pour écouter l'onglet actif
                    selected = curr_tab,         # NEW : Restauration immédiate de l'état

                    nav_panel(
                        title = "Description",
                        value = "tab_description", # NEW : Identifiant interne du panneau
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
                        title = "Biplot",
                        value = "tab_biplot",      # NEW : Identifiant interne du panneau
                        fluidRow(
                            # ---- 1. Contrôles & Axes ----
                            column(
                                3,
                                grid_card(
                                    area = "sidebar",
                                    h3("Controls"),
                                    selectInput(ns("xvar"),
                                                "Axe X",
                                                choices = rownames(ScIGMA_data$mae[["proteins"]])),
                                    selectInput(ns("yvar"),
                                                "Axe Y",
                                                choices = rownames(ScIGMA_data$mae[["proteins"]])),
                                    checkboxInput(ns("logx"),
                                                  "Log X",
                                                  FALSE),
                                    checkboxInput(ns("logy"),
                                                  "Log Y",
                                                  FALSE),
                                    selectInput(ns("color_genotype"),
                                                "Color by Genotype",
                                                # choices = c("None", ScIGMA_data$variants.filtered$label)),
                                                choices = c("None", levels(ScIGMA_data$dna.clones))),

                                    # FIX CRITIQUE : Détection propre du nouvel Assay
                                    # if ("compass_imputed" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["dna_variants"]])) {
                                    #     shinyWidgets::materialSwitch(
                                    #         inputId = ns("use_compass_gt"),
                                    #         label = "Use COMPASS imputed data",
                                    #         value = FALSE,
                                    #         status = "success"
                                    #     )
                                    # },
                                    hr(),
                                    h4("Gating"),
                                    textInput(ns("subset_name"), "Gate name", placeholder = "ex: Gate A"),
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
                    nav_panel(
                        title = "UMAP",
                        value = "tab_umap",        # NEW : Identifiant interne du panneau
                        accordion(
                            id = ns("umap_accordion_main"),
                                  accordion_panel(
                                      "Built UMAP model",
                                      fluidRow(
                                          column(3,
                                                 tagList(
                                                     grid_card(
                                                         area = "umap_accordion_first",
                                                         h3("UMAP parameters"),

                                                         shinyWidgets::pickerInput(
                                                             inputId = ns("umap_features"),
                                                             label = "Protein markers",
                                                             choices = rownames(ScIGMA_data$mae[["proteins"]]),
                                                             selected = curr_features,
                                                             multiple = TRUE,
                                                             options = shinyWidgets::pickerOptions(
                                                                 actionsBox = TRUE,
                                                                 liveSearch = TRUE,
                                                                 size = 10,
                                                                 `selected-text-format` = "count > 3"
                                                             )
                                                         ),

                                                         numericInput(ns("n_neighbors"), "Neighbors", value = curr_n_neighbors, min = 5, max = 100),
                                                         numericInput(ns("min_dist"), "Distance Min", value = curr_min_dist, min = 0.01, max = 1.0, step = 0.1),
                                                         hr(),
                                                         actionBttn(
                                                             inputId = ns("run_umap_btn"), label = "Run UMAP",
                                                             style = "unite", color = "primary"
                                                         ),
                                                         helpText("The calculation may take a few seconds."),

                                                         uiOutput(ns("umap_evaluation_panel")),
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

                                  # NEW : Remplacement du contenu si l'UMAP n'est pas calculée
                                  accordion_panel(
                                      "Markers expression",
                                      if (umap_ready) {
                                          uiOutput(ns("markers_umap_panel_ui"))
                                      } else {
                                          fluidRow(h3("Please compute UMAP first", style = "text-align: center; color: #7f8c8d; margin-top: 20px;"))
                                      }
                                  ),
                                  accordion_panel(
                                      "Biplot gates",
                                      if (umap_ready) {
                                          uiOutput(ns("biplotClones_umap_panel_ui"))
                                      } else {
                                          fluidRow(h3("Please compute UMAP first", style = "text-align: center; color: #7f8c8d; margin-top: 20px;"))
                                      }
                                  ),
                                  accordion_panel(
                                      "Unsupervised clustering",
                                      if (umap_ready) {
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
                                      } else {
                                          fluidRow(h3("Please compute UMAP first", style = "text-align: center; color: #7f8c8d; margin-top: 20px;"))
                                      }
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
                 watch("compass_completed"), # <-- FIX CRITIQUE 1 : Reset synchronisé sur la purge de l'entonnoir
                 input$reset_root
            ),
            {
                req(ScIGMA_data$mae)

                # Sécurité : On s'assure que les protéines sont bien dans le dataset
                if (!("proteins" %in% names(ScIGMA_data$mae))) return()

                assay_to_use_state <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")
                n_cells_current <- ncol(SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use_state))

                r_state$subsets <- list()
                r_state$subsets[["root"]] <- 1:n_cells_current

                r_state$subset_meta <- list()
                r_state$subset_meta[["root"]] <- list(
                    name = "All",
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

        shiny::observeEvent(
            list(
                watch("dnaVariant_selected"),  # Quand une heatmap est générée
                watch("dna_clones_renamed"),   # Quand tu renommes manuellement
                watch("compass_completed")     # Quand COMPASS recalcule l'architecture
            ),
            {
                # Sécurité : on ne fait rien si les clones n'existent pas encore
                shiny::req(ScIGMA_data$dna.clones)

                # 1. Extraction brute : on contourne les bugs liés aux "Factors" en R
                current_clones <- unique(as.character(ScIGMA_data$dna.clones))

                # 2. Nettoyage strict : on vire le bruit technique de la liste déroulante
                current_clones <- setdiff(current_clones, c("Missing", "Missing/ADO", "small", "NA", "Unknown"))

                # 3. Tri propre
                current_clones <- sort(current_clones)

                # 4. Injection dans l'interface
                shiny::updateSelectInput(
                    session = session,
                    inputId = "color_genotype",  # CRITIQUE : Jamais de ns() dans le serveur pour un update
                    choices = c("None", current_clones)
                )
            },
            ignoreInit = TRUE
        )

        output$biplot <- renderPlotly({
            req(input$xvar, input$yvar, r_state$current_view)

            current_indices <- r_state$subsets[[r_state$current_view]]
            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            raw_x <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[input$xvar, current_indices]
            raw_y <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[input$yvar, current_indices]

            plot_df <- data.frame(
                x = raw_x, y = raw_y, custom_id = current_indices
            )

            if (input$logx) plot_df$x <- log1p(plot_df$x)
            if (input$logy) plot_df$y <- log1p(plot_df$y)

            # Configuration par défaut (Aucune coloration)
            color_formula <- NULL
            color_palette <- NULL
            marker_config <- list(size = 5, opacity = 0.8, color = "#2c3e50")

            if (!is.null(input$color_genotype) && input$color_genotype != "None") {

                # 1. Mappage strict du clone pour chaque cellule
                plot_df$Genotype <- as.character(ScIGMA_data$dna.clones[rownames(plot_df)])

                # 2. Isolement du clone cible, tout le reste devient "Other"
                plot_df$Genotype <- ifelse(plot_df$Genotype %in% input$color_genotype, plot_df$Genotype, "Other")

                # 3. Verrouillage factoriel (Z-Indexing)
                # En mettant "Other" en premier, Plotly va dessiner le nuage gris au fond,
                # et le clone coloré par-dessus.
                plot_df$Genotype <- factor(plot_df$Genotype, levels = c("Other", input$color_genotype))

                color_formula <- ~Genotype

                # 4. Injection de la transparence via le canal Alpha Hexadécimal (33 = 20%)
                color_palette <- c(ScIGMA_data$dna_clone_colors, "Other" = "#e0e0e033")

                # 5. On retire l'opacité et la couleur forcées pour laisser la palette agir
                marker_config <- list(size = 5)
            }

            plotly::plot_ly(
                data = plot_df,
                x = ~x,
                y = ~y,
                key = ~custom_id,
                color = color_formula,
                colors = color_palette,
                type = "scattergl",
                mode = "markers",
                source = ns("gating_plot"),
                marker = marker_config
            ) %>%
                plotly::layout(
                    title = list(text = paste("<b>Gate:</b>", r_state$subset_meta[[r_state$current_view]]$name),
                                 font = list(family = "Arial", size = 18)),
                    plot_bgcolor = "white",
                    paper_bgcolor = "white",
                    xaxis = c(list(title = paste("<b>", input$xvar, "</b>"), range = list(0, max(plot_df$x)+1)), prism_axis_style),
                    yaxis = c(list(title = paste("<b>", input$yvar, "</b>"), range = list(0, max(plot_df$y)+1)), prism_axis_style),
                    dragmode = "lasso",
                    legend = list(title = list(text = "<b>Genotype</b>")),
                    margin = list(l = 60, r = 30, b = 60, t = 50)
                ) %>%
                plotly::config(displaylogo = FALSE)
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

            safe_name <- sanitize_gate_name(input$subset_name)

            if (safe_name == "") {
                # On génère un nom par défaut qui respecte aussi la convention sans espace
                safe_name <- paste0("Gate_", length(r_state$subsets))
            }

            # name_val <- input$subset_name
            # if (name_val == "") name_val <- paste("Gate", length(r_state$subsets))

            # Save Indices & Meta
            r_state$subsets[[new_id]] <- r_state$temp_selection

            parent_id <- r_state$current_view
            r_state$subset_meta[[new_id]] <- list(
                # name = name_val,
                name = safe_name,
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
                    actionLink(ns(paste0("go_", node_id)), label = meta$name, style = style_str)
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

            trigger("gating_updated")
        })

        # --- 5. UMAP Calculation ---

        # Trigger global via bouton (recommandé pour éviter les calculs intempestifs)
        observeEvent(input$run_umap_btn, {
            gargoyle::trigger("launch_umap")
        })

        # Bloc de calcul réactif
        observe({
            watch("launch_umap")

            w <- waiter::Waiter$new(
                html = waiter::spin_loaders(1, color = "black"),
                color = "rgba(255, 255, 255, 0.5)"
            )
            w$show()

            n_neighbors <- isolate(input$n_neighbors)
            min_dist <- isolate(input$min_dist)

            # NEW : Extraction des features sélectionnées
            umap_features <- isolate(input$umap_features)

            req(ScIGMA_data$mae[["proteins"]])

            # Sécurité : Éviter un crash si l'utilisateur désélectionne tout
            if (is.null(umap_features) || length(umap_features) < 3) {
                shiny::showNotification("Please select at least 3 protein markers to build a UMAP model.", type = "error")
                w$hide()
                return()
            }

            message("Computing UMAP...")

            # UPDATED : Injection stricte des features ciblées dans Seurat
            ScIGMA_data$seurat_object <- Seurat::RunUMAP(ScIGMA_data$seurat_object,
                                                         features = umap_features,
                                                         min.dist = min_dist,
                                                         n.neighbors = n_neighbors)

            # --- FIX CRITIQUE : Vraie Ground Truth ---
            message("Computing rigorous UMAP metrics against raw normalized data...")

            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            # UPDATED : Restriction de la matrice haute-dimension aux features sélectionnées
            true_ground_truth_mat <- t(SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[umap_features, , drop = FALSE])
            umap_mat <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings

            # Calcul du triptyque
            metrics <- compute_umap_metrics(
                high_dim_mat = true_ground_truth_mat,
                low_dim_mat = umap_mat,
                k = n_neighbors
            )

            # Stockage pour l'UI
            ScIGMA_data$umaps$metrics <- list(
                trustworthiness = round(metrics$trustworthiness * 100, 1),
                continuity = round(metrics$continuity * 100, 1),
                spearman = round(metrics$global_spearman, 3)
            )
            # ----------------------------------------------------------------------

            w$hide()
            gargoyle::trigger("umap_computed")
            message("Calcul UMAP et évaluations terminés.")

        }) |>
            shiny::bindEvent(input$run_umap_btn)

        # --- NEW : Rendu dynamique du panneau de score (Radar Plot) ---
        output$umap_evaluation_panel <- renderUI({
            watch("umap_computed")
            req(ScIGMA_data$umaps$metrics)

            div(
                style = "margin-top: 15px; padding: 5px; background-color: #f8f9fa; border-radius: 5px; border: 1px solid #dee2e6;",
                plotlyOutput(ns("umap_radar_plot"), height = "280px")
            )
        })

        # --- NEW : Génération du Spider Plot (Norme Q1/Prism-like) ---
        output$umap_radar_plot <- renderPlotly({
            watch("umap_computed")
            req(ScIGMA_data$umaps$metrics)

            metrics <- ScIGMA_data$umaps$metrics

            # Mise à l'échelle : Spearman [-1, 1] transformé en [0, 100] pour cohérence visuelle
            spearman_scaled <- max(0, metrics$spearman * 100)

            # Fermeture de la boucle pour Plotly (le premier point doit être répété à la fin)
            r_values <- c(metrics$trustworthiness, metrics$continuity, spearman_scaled, metrics$trustworthiness)

            # FIX CRITIQUE 1 : Utilisation des balises HTML <b> et <br> pour le style
            theta_labels <- c('<b>Local</b><br>(Trustworthiness)',
                              '<b>Structure</b><br>(Continuity)',
                              '<b>Global</b><br>(Topology)',
                              '<b>Local</b><br>(Trustworthiness)')

            plot_ly(
                type = 'scatterpolar',
                r = r_values,
                theta = theta_labels,
                fill = 'toself',
                fillcolor = 'rgba(41, 128, 185, 0.3)', # Bleu transparent élégant
                line = list(color = '#2c3e50', width = 2),
                marker = list(color = '#e74c3c', size = 6)
            ) %>%
                layout(
                    polar = list(
                        radialaxis = list(
                            visible = TRUE,
                            range = c(0, 100),
                            tickfont = list(family = "Arial", size = 10, color = "black"),
                            gridcolor = "#e0e0e0",
                            linecolor = "black"
                        ),
                        angularaxis = list(
                            tickfont = list(family = "Arial", size = 11, color = "black"),
                            linecolor = "black"
                        )
                    ),
                    margin = list(t = 40, b = 20, l = 40, r = 40),
                    paper_bgcolor = "rgba(0,0,0,0)", # Fond transparent
                    plot_bgcolor = "rgba(0,0,0,0)",
                    showlegend = FALSE,
                    title = list(
                        text = "<b>Model Fidelity (%)</b>",
                        font = list(family = "Arial", size = 14, color = "black"),
                        y = 0.95
                    )
                ) %>%
                config(displaylogo = FALSE, displayModeBar = FALSE) # Interface épurée
        })



        # render UMAP without clustering
        observeEvent({watch("umap_computed")},{
            req(ScIGMA_data$seurat_object)
            umap_cluster <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                as.data.frame()
            umap_cluster$barcode <- rownames(umap_cluster)

            # --- NEW : Calcul du Top Marqueurs (Haute Performance & Défensif) ---
            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            # Matrice transposée (Cellules x Protéines)
            mat <- t(SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use))
            mat <- mat[umap_cluster$barcode, , drop = FALSE] # Sécurité d'alignement

            ptn_names <- colnames(mat)
            n_top <- min(3, ncol(mat)) # Sécurité si < 3 marqueurs

            # Calcul "one-shot" pour éviter de saturer la RAM ou le rendu
            top_strings <- apply(mat, 1, function(x) {
                idx <- order(x, decreasing = TRUE)[1:n_top]
                paste(paste0("<b>", 1:n_top, ". ", ptn_names[idx], "</b>: ", round(x[idx], 2)), collapse = "<br>")
            })

            umap_cluster$hover_info <- paste0(
                "<b>--- Top Expression ---</b><br>",
                top_strings
            )
            # -------------------------------------------------------------

            # Version Sauvegardée (ggplot2 strict)
            ScIGMA_data$umaps$umap_protein_general <- umap_cluster |>
                ggplot(aes(x=umap_1, y=umap_2)) +
                geom_point(size=1.5, color = "#2c3e50", alpha = 0.8) +
                xlab('UMAP 1') +
                ylab('UMAP 2') +
                ggprism::theme_prism(base_size = 14, base_fontface = "bold")

            # Version Interactive (Plotly)
            output$umap_plot_build <- renderPlotly({
                p <- plot_ly(data = umap_cluster,
                             x = ~umap_1,
                             y = ~umap_2,
                             type = 'scattergl',
                             mode = 'markers',
                             text = ~hover_info,
                             hoverinfo = "text",
                             marker = list(size = 5, opacity = 0.8, color = "#2c3e50")) %>%
                    layout(
                        plot_bgcolor = "white",
                        paper_bgcolor = "white",
                        xaxis = list(visible = FALSE),
                        yaxis = list(visible = FALSE),
                        margin = list(l = 60, r = 30, b = 10, t = 30)
                    ) %>%
                    config(displaylogo = FALSE)
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
                ScIGMA_data$seurat_object <- FindNeighbors(
                    ScIGMA_data$seurat_object,
                    features = rownames(ScIGMA_data$seurat_object),
                    dims = 1:ncol(ScIGMA_data$seurat_object@reductions$pca@cell.embeddings)
                )

                ScIGMA_data$seurat_object <- FindClusters(
                    ScIGMA_data$seurat_object,
                    resolution = clustering_resolution
                )

                umap_cluster <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |>
                    as.data.frame()
                umap_cluster$cluster = ScIGMA_data$seurat_object$seurat_clusters[rownames(umap_cluster)]

                # 2. Rendering-
                output$umap_clustering_plot <- renderPlotly({

                    w <- Waiter$new(
                        id = ns("umap_clustering_plot"),
                        html = spin_3(),
                        color = transparent(0.5)
                    )
                    w$show()

                    # Palette Viridis pour variables catégorielles (Clusters)
                    n_clusters <- length(unique(umap_cluster$cluster))
                    pal <- viridis::viridis(n_clusters)

                    p <- plot_ly(data = umap_cluster,
                                 x = ~umap_1,
                                 y = ~umap_2,
                                 type = 'scattergl',
                                 mode = 'markers',
                                 color = ~cluster,
                                 colors = pal,
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

                trigger("clusters_computed")

            },ignoreInit = TRUE)



        # [ NODE_ACCESS : Markers projection ]
        # ----------------------------------------------------- _
        observeEvent(watch("umap_computed"), {
            req(ScIGMA_data$seurat_object)

            if(is.null(ScIGMA_data$seurat_object@reductions$umap)){
                output$markers_umap_panel_ui <-  renderUI({
                    tagList(fluidRow(h2("Please compute UMAP first")))
                })
            } else {
                # FIX CRITIQUE : Le cbind destructeur a été atomisé.
                # On génère l'UI proprement à partir du MAE.
                protein_names <- rownames(ScIGMA_data$mae[["proteins"]])

                output$markers_umap_panel_ui <-  renderUI({
                    fluidRow(
                        tagList(
                            column(3,
                                   grid_card(
                                       area = "umap_marker_select",
                                       virtualSelectInput(
                                           inputId = ns("protein_umap_markers"),
                                           label = "Markers :",
                                           choices = protein_names,
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
                            column(9, plotOutput(ns("umap_plot_markers"), height = "600px"))
                        )
                    )
                })
            }
        }, ignoreInit = TRUE)

        # [ NODE_ACCESS : Biplot clones projection ]
        # ----------------------------------------------------- _
        observeEvent(
            list(watch("umap_computed"),
                 input$save_to_r6), # Actualise le menu déroulant quand on sauvegarde de nouveaux clones
            {
                req(ScIGMA_data$seurat_object)

                if(is.null(ScIGMA_data$seurat_object@reductions$umap)){
                    output$biplotClones_umap_panel_ui <- renderUI({
                        tagList(fluidRow(h2("Please compute UMAP first")))
                    })
                    return()
                }

                # Récupération de l'arbre sauvegardé
                gating_tree <- ScIGMA_data$protein_gating_tree
                if (is.null(gating_tree) || length(gating_tree$meta_list) <= 1) {
                    output$biplotClones_umap_panel_ui <- renderUI({
                        tagList(fluidRow(h4("No sub-samples saved yet. Please create and save Gates in the Bi-plot tab.", style="color:#7f8c8d; padding:20px;")))
                    })
                    return()
                }

                # Préparation des choix (Exclure la racine "Tout")
                meta_list <- gating_tree$meta_list
                valid_ids <- setdiff(names(meta_list), "root")

                # Formatage : Value = ID (ex: sub_123), Name = Nom utilisateur (ex: TOM1)
                choices_list <- setNames(valid_ids, sapply(valid_ids, function(id) meta_list[[id]]$name))

                output$biplotClones_umap_panel_ui <- renderUI({
                    fluidRow(
                        tagList(
                            column(3,
                                   grid_card(
                                       area = "umap_clone_select",
                                       virtualSelectInput(
                                           inputId = ns("umap_clones_to_project"),
                                           label = "Select sub-gates:",
                                           choices = choices_list,
                                           multiple = TRUE,
                                           width = "100%",
                                           dropboxWrapper = "body"
                                       ),
                                       hr(),
                                       actionBttn(
                                           inputId = ns("project_clones_btn"),
                                           label = "Project gates",
                                           style = "unite",
                                           color = "primary"
                                       ),
                                       helpText("Child gates are automatically drawn on top of their parents.")
                                   )
                            ),
                            column(9, plotlyOutput(ns("umap_plot_biplot_clones"), height = "600px"))
                        )
                    )
                })
            }, ignoreInit = TRUE)

        observeEvent(input$project_clones_btn, {
            req(ScIGMA_data$seurat_object@reductions$umap)
            req(input$umap_clones_to_project)

            selected_ids <- input$umap_clones_to_project
            gates_list <- ScIGMA_data$protein_gating_tree$gates_list
            meta_list <- ScIGMA_data$protein_gating_tree$meta_list
            all_cell_barcodes <- colnames(ScIGMA_data$mae[["proteins"]])

            # Extraction UMAP
            umap_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |> as.data.frame()
            umap_df$barcode <- rownames(umap_df)
            umap_df$Clone <- "Background"

            # 1. Tri par profondeur (Depth-First)
            # Les enfants (depth plus grand) seront traités en dernier pour écraser le label des parents
            depths <- sapply(selected_ids, function(id) meta_list[[id]]$depth)
            ordered_ids <- selected_ids[order(depths)]

            # 2. Assignation des labels
            for (sid in ordered_ids) {
                cell_indices <- gates_list[[sid]]
                cell_barcodes <- all_cell_barcodes[cell_indices]
                valid_barcodes <- intersect(cell_barcodes, umap_df$barcode)

                clone_name <- meta_list[[sid]]$name
                umap_df[valid_barcodes, "Clone"] <- clone_name
            }

            # 3. Factorisation stricte pour garantir l'ordre de dessin dans Plotly
            # On veut que le "Background" soit dessiné en premier (z-index le plus bas)
            clone_levels <- c("Background", sapply(ordered_ids, function(id) meta_list[[id]]$name))
            umap_df$Clone <- factor(umap_df$Clone, levels = unique(clone_levels))
            umap_df <- umap_df[order(umap_df$Clone), ]

            # 4. Palette de couleurs (Gris fixe pour le fond, couleurs discrètes pour les clones)
            n_clones <- length(unique(clone_levels)) - 1
            pal <- if(n_clones > 0) scales::hue_pal()(n_clones) else c()
            color_map <- setNames(c("#e0e0e0", pal), unique(clone_levels))

            # 5. Rendu WebGL
            p <- plot_ly(
                data = umap_df,
                x = ~umap_1,
                y = ~umap_2,
                color = ~Clone,
                colors = color_map,
                type = 'scattergl',
                mode = 'markers',
                marker = list(size = 5, opacity = 0.8),
                text = ~paste("<b>Gate:</b>", Clone),
                hoverinfo = "text"
            ) %>%
                layout(
                    plot_bgcolor = "white",
                    paper_bgcolor = "white",
                    xaxis = list(visible = FALSE),
                    yaxis = list(visible = FALSE),
                    legend = list(
                        title = list(text = "<b>Sub-gates</b>", font = list(family = "Arial", color = "black")),
                        font = list(family = "Arial", size = 12, color = "black")
                    ),
                    margin = list(l = 60, r = 30, b = 10, t = 30)
                ) %>%
                config(displaylogo = FALSE)

            output$umap_plot_biplot_clones <- renderPlotly({ p })
        })



        # [ NODE_ACCESS : PROTEIN MARKERS ]
        # ----------------------------------------------------- _
        observeEvent(input$protein_umap_markers_bttn, {
            req(input$protein_umap_markers)
            umap_marker_choice <- input$protein_umap_markers

            # Extraction de la base UMAP
            umap_df <- ScIGMA_data$seurat_object@reductions$umap@cell.embeddings |> as.data.frame()
            umap_df$barcode <- rownames(umap_df)

            assay_to_use <- ifelse("normalized" %in% SummarizedExperiment::assayNames(ScIGMA_data$mae[["proteins"]]), "normalized", "counts")

            if (length(umap_marker_choice) == 0 || "None" %in% umap_marker_choice) {
                umap_df$ptn_expression <- 0
                umap_df$marker <- "No marker selected"
            } else {
                # 1. On extrait uniquement les protéines sélectionnées (Protéines x Cellules)
                expr_mat <- SummarizedExperiment::assay(ScIGMA_data$mae[["proteins"]], assay_to_use)[umap_marker_choice, , drop = FALSE]

                # 2. Transposition (Cellules x Protéines)
                expr_df <- t(expr_mat) |> as.data.frame()
                expr_df$barcode <- rownames(expr_df)

                # 3. Jointure sécurisée (UMAP + Expression)
                umap_df <- merge(umap_df, expr_df, by = "barcode")

                # 4. Pivotage pour le facet_wrap de ggplot2
                umap_df <- umap_df |>
                    tidyr::pivot_longer(cols = tidyr::all_of(umap_marker_choice),
                                        values_to = "ptn_expression",
                                        names_to = "marker")
            }

            ScIGMA_data$umaps$umap_protein_markers <- umap_df |>
                ggplot(aes(x=umap_1, y=umap_2, color=ptn_expression)) +
                geom_point(size = 1.2, alpha = 0.9) +
                facet_wrap(~marker) +
                scale_color_viridis_c("inferno", name = "Expression") +
                xlab("UMAP 1") +
                ylab("UMAP 2") +
                ggprism::theme_prism(base_size = 14, base_fontface = "bold")

            output$umap_plot_markers <- renderPlot(ScIGMA_data$umaps$umap_protein_markers)
        }, ignoreInit = TRUE)

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
