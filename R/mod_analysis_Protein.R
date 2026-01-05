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
                            h3("Contrôles"),
                            # Sélection des marqueurs (rempli côté serveur via R6)
                            selectInput(ns("xvar"), "Axe X", choices = NULL),
                            selectInput(ns("yvar"), "Axe Y", choices = NULL),
                            checkboxInput(ns("logx"), "Log X", FALSE),
                            checkboxInput(ns("logy"), "Log Y", FALSE),
                            hr(),
                            h4("Gating"),
                            textInput(ns("subset_name"), "Nom du gate", placeholder = "ex: Clone A"),
                            actionButton(ns("mk_subset"), "Créer sous-échantillon ↘︎", class = "btn-primary"),
                            verbatimTextOutput(ns("selection_info"))
                        )
                    ),
                    # ---- 2. Visualisation (Bi-plot) ----
                    column(
                        6,
                        grid_card(
                            area = "main",
                            uiOutput(ns("biplot_container"))
                            # Utilisation d'un UI container pour gérer la hauteur dynamiquement si besoin,
                            # ou directement plotlyOutput ici.
                        )
                    ),
                    # ---- 3. Arborescence (Tree) ----
                    column(
                        3,
                        grid_card(
                            area = "subsets",
                            h3("Sous-échantillons"),
                            uiOutput(ns("subsets_ui")),
                            hr(),
                            actionButton(ns("reset_root"), "Retour Racine", class = "btn-warning")
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
        observe({
            watch("dnaVariant_filtered")
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
            raw_x <- ScIGMA_data$protein.mtx[current_indices, input$xvar]
            raw_y <- ScIGMA_data$protein.mtx[current_indices, input$yvar]

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
                marker = list(size = 3, opacity = 0.6, color = "#2c3e50")
            ) %>%
                layout(
                    title = paste("Gate:", r_state$subset_meta[[r_state$current_view]]$name),
                    xaxis = list(title = input$xvar),
                    yaxis = list(title = input$yvar),
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
            ScIGMA_data <- normalizeProtein(ScIGMA_data)
            render_protein_ridge_plot(obj = ScIGMA_data)
        })

        # Render Barplot for all proteins
        output$protein_bar <- renderPlotly({
            watch("dnaVariant_filtered")
            ScIGMA_data <- normalizeProtein(ScIGMA_data)
            plot_protein_barplot(obj = ScIGMA_data)
        })
    })
}

## To be copied in the UI
# mod_analysis_Protein_ui("analysis_Protein_1")

## To be copied in the server
# mod_analysis_Protein_server("analysis_Protein_1")
