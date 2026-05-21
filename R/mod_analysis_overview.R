#' analysis_right_overview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' @import waiter
mod_analysis_overview_ui <- function(id) {
    ns <- NS(id)
    ## Add a spinner
    add_busy_spinner(spin = "fading-circle", color = "#112446")
    tagList(
        fluidRow(
            column(4,
                   card(
                       fileInput(ns("file_h5file"),
                                 label = "1. Upload you H5 file.",
                                 accept = ".h5")
                   )
            ),
            column(4,
                   card(
                       radioGroupButtons(
                           inputId = ns("file_fileType"),
                           label = "2. DNA or DNA+protein ?",
                           choices = c("DNA only"="DNA",
                                       "DNA & protein"="DNA+protein"),
                           justified = TRUE
                       )
                   )
            ),
            column(4,
                   h6(HTML("3. Process")),
                   actionBttn(
                       inputId = ns("file_process"),
                       label = "Load file",
                       color = "success",
                       style = "unite",
                       icon = icon("magnifying-glass-chart"),
                       block = TRUE
                   )
            )
        ),
        fluidRow(
            br(),
            card(card_header("Summary"),
                 uiOutput(ns("overview"))
            ),
            card(
                card_header("Preprocess"),
                uiOutput(ns("preprocess")),
                uiOutput(ns("dnaFilterResults"))
            )
        )
    )
}

#' analysis_right_overview Server Functions
#'
#' @noRd
#' @import tools
#' @import Seurat
#' @import MultiAssayExperiment
mod_analysis_overview_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
        message(whereami::whereami())

        # --------------------------------------------------------------- #
        # Reactive state for sharing metrics between observeEvent and renderUI
        init_metrics <- reactiveValues(
            init_number_cell = NULL,
            init_number_dna_variant = NULL
        )

        # ---------------------------- #
        # Load data
        observeEvent(input$file_process,{
            filePath <- input$file_h5file$datapath
            fileType <- input$file_fileType
            req(filePath)
            req(fileType)
            sampleName <- file_path_sans_ext(input$file_h5file$name)
            message(whereami::whereami())
            # show_modal_spinner(text = "Loading data ...")
            w <- Waiter$new(
                html = spin_double_bounce(),             # Spinner discret et élégant
                color = "rgba(255, 255, 255, 0.5)" # Fond blanc semi-transparent (pas de bloc blanc opaque)
            )
            w$show()
            if (file.exists(filePath)) {
                if (file.info(filePath)$isdir) {
                    ScIGMA_data <- tryCatch(
                        loadH5_HDF5_biocond(
                            filepath = filePath,
                            sample_name = sampleName,
                            omic_type = fileType
                        ),
                        error = function(e){
                            message("Error during loadH5")
                            stop(e$message)
                        })
                } else {
                    temp_scigma_obj <- tryCatch(
                        loadH5_HDF5_biocond(
                            filepath = filePath,
                            sample_name = sampleName,
                            omic_type = fileType
                        ),
                        error = function(e){
                            message("Error during loadH5")
                            stop(e$message)
                        }
                    )
                }
            } else {
                stop("File or folder path doesn't exists\n")
            }
            # 2. On transfère les données dans l'objet R6 GLOBAL partagé
            # Cela préserve la référence mémoire pour tous tes modules Shiny
            ScIGMA_data$mae <- temp_scigma_obj$mae
            ScIGMA_data$mae_raw <- temp_scigma_obj$mae
            ScIGMA_data$filetype <- temp_scigma_obj$filetype

            ScIGMA_data$reset_analysis()

            # Remise à zéro des métriques UI locales
            init_metrics$init_number_cell <- NULL
            init_metrics$init_number_dna_variant <- NULL

            # Injection des nouvelles données
            ScIGMA_data$mae <- temp_scigma_obj$mae
            ScIGMA_data$mae_raw <- temp_scigma_obj$mae
            ScIGMA_data$filetype <- temp_scigma_obj$filetype

            w$hide()
            message(whereami::whereami())

            trigger("dataLoaded")
        })

        # ---------------------------- #
        # Render Summary UI
        output$overview <- renderUI({
            watch("dataLoaded")
            print("ScIGMA_data render UI")
            print(ScIGMA_data)
            message(whereami::whereami())

            has_mae <- !is.null(ScIGMA_data$mae)

            val_cells <- if (has_mae) nrow(SummarizedExperiment::colData(ScIGMA_data$mae)) else "No data loaded"
            val_dna   <- if (has_mae) nrow(ScIGMA_data$mae[["dna_variants"]]) else "No data loaded"
            val_cnv   <- if (has_mae) nrow(ScIGMA_data$mae[["amplicons"]]) else "No data loaded"

            # Vérification stricte de l'existence de l'assay 'proteins'
            val_prot <- if (!has_mae) {
                "No data loaded"
            } else if (!("proteins" %in% names(ScIGMA_data$mae))) {
                "No Protein data"
            } else {
                ScIGMA_object$protein.filtered <- TRUE
                nrow(ScIGMA_data$mae[["proteins"]])
            }

            fluidRow(
                column(3,
                       card(
                           card_header("Number of cells"),
                           card_body(
                               p(val_cells, style = "text-align:center")
                           )
                       )
                ),
                column(3,
                       card(
                           card_header("DNA variants"),
                           p(val_dna, style = "text-align:center")
                       )
                ),
                column(3,
                       card(
                           card_header("Number of amplicons"),
                           p(val_cnv, style = "text-align:center")
                       )
                ),
                column(3,
                       card(
                           card_header("Number of proteins"),
                           p(val_prot, style = "text-align:center")
                       )
                )
            )
        })

        # ---------------------------- #
        # Render Preprocess UI
        output$preprocess <- renderUI({
            watch("dataLoaded")
            message(whereami::whereami())
            if(!is.null(ScIGMA_data)){
                tagList(
                    fluidRow(
                        h5("Filter Cells and DNA variants:"),

                        div(
                            HTML("DNA Variant filtering step has two parameters : </br>"),
                            align = "justify"
                        )
                    ),
                    fluidRow(
                        column(6,
                               div(
                                   HTML("<b>min.cell.pt</b>	Minimum threshold for cell percentage that has valid
                                   variant call (GT = 0, 1 or 2) after
              applying the filter. When to change: If the variant of interest is in a high GC content
              region, then PCR amplification is hard. In such cases, you may choose to decrease the percent to 30 or 40
              so that your interested variant could come through the filter.</br>"),
                                   align = "justify"
                               ),
                               div(
                                   sliderTextInput(
                                       inputId = ns("overview_preprocess_minCellPt"),
                                       label = "min.cell.pt",
                                       choices = seq(1, 100,1),
                                       grid = TRUE,
                                       selected = 50
                                   ),
                                   align = "center")
                        ),
                        column(6,
                               div(
                                   HTML("<b>min.mut.cell.pt</b>
                                   Minimum threshold for cell percentage that has mutated genotype (GT = 1 or 2) after
              applying the filter.
              When to change: If you know the variant is rare in the data, then you could try lower threshold to try to
              keep the variant in your dataset."),
                                   align = "justify"
                               ),
                               div(
                               sliderTextInput(
                                   inputId = ns("overview_preprocess_minMutCellPt"),
                                   label = "min.mut.cell.pt",
                                   choices = seq(1, 30, 1),
                                   grid = TRUE,
                                   selected = 1
                               ),
                               align = "center")
                        )
                    ),
                    fluidRow(
                        actionBttn(
                            inputId = ns("dna_variant_filtering"),
                            label = "Filter Cells and DNA variants",
                            color = "primary",
                            style = "stretch",
                            icon = icon("magnifying-glass-chart"),
                            block = TRUE
                        )
                    )
                )
            } else {
                fluidRow(
                    h5("No data provided")
                )
            }

        })

        # --------------------------------------------------------------- #
        # Filter DNA variant
        observeEvent(input$dna_variant_filtering, {
            filePath <- input$file_h5file$datapath
            overview_preprocess_minCellPt <- input$overview_preprocess_minCellPt
            overview_preprocess_minMutCellPt <- input$overview_preprocess_minMutCellPt

            req(overview_preprocess_minCellPt)
            req(overview_preprocess_minMutCellPt)
            req(ScIGMA_data$mae) # Sécurité : s'assure que les données sont chargées

            message(whereami::whereami())
            # show_modal_spinner(text = "Filtering and annotating DNA variants...")
            w <- Waiter$new(
                # id = "umap_clustering_plot", # Peut cibler un plot précis ou toute la page
                html = spin_loaders(1, color = "black"),             # Spinner discret et élégant
                color = "rgba(255, 255, 255, 0.5)" # Fond blanc semi-transparent (pas de bloc blanc opaque)
            )
            w$show()

            # ---------------------------- #
            # Store initial cell and DNA variant info (Utilisation du MAE)
            init_metrics$init_number_cell <- ncol(ScIGMA_data$mae)
            init_metrics$init_number_dna_variant <- nrow(ScIGMA_data$mae[["dna_variants"]])

            message(whereami::whereami())

            # ---------------------------- #
            # Execute Pipeline: Filter + Annotate + Proportions
            temp_scigma_obj <- tryCatch({
                filter_and_annotate_variants(
                    obj = ScIGMA_data,
                    paths = cfg$paths, # Assure-toi que cfg$paths est bien accessible ici
                    min_cell_pt = overview_preprocess_minCellPt,
                    min_mut_cell_pt = overview_preprocess_minMutCellPt
                    # Les autres paramètres prendront leurs valeurs par défaut (min_dp=10, etc.)
                )
            }, error = function(e) {
                remove_modal_spinner()
                message("Error during DNA variant filtering and annotation")
                stop(e$message)
            })

            # ---------------------------- #
            # Update Global R6 Reference
            # On transfère le MAE filtré et annoté dans l'objet global
            ScIGMA_data$mae <- temp_scigma_obj$mae

            # >> Sanitize MAE _
            ScIGMA_data$mae <- sanitize_mae_strings(ScIGMA_data$mae)


            # Normalize protein and perform PCA
            if (ScIGMA_data$filetype == "DNA+protein"){
                message("Preprocessing protein data ...")

                #ScIGMA_data <- normalizeProtein(ScIGMA_data)
                ScIGMA_data$seurat_object <- protein_run_pca(ScIGMA_data)
            }

            ScIGMA_data$cnv_dp_filtered <- NULL

            message(whereami::whereami())

            # ---------------------------- #
            # Trigger downstream modules
            # remove_modal_spinner()
            w$hide()
            trigger("dnaVariant_filtered")
            trigger("launch_umap")
        })


        # --------------------------------------------------------------- #
        # Render UI after DNA variant filtering
        output$dnaFilterResults <- renderUI({
            watch("dnaVariant_filtered")
            watch("dataLoaded")
            message(whereami::whereami())

            # Sécurité anti-crash
            if (is.null(ScIGMA_data$mae)) return(NULL)

            filter_status <- S4Vectors::metadata(ScIGMA_data$mae)$variant_filter

            if (is.null(filter_status) || filter_status != "filtered") {
                fluidRow(
                    h5("DNA variant filtering results:"),
                    div("Data not filtered yet", align ="center")
                )
            } else {
                # Extraction des dimensions POST-filtrage
                current_cells <-ncol(ScIGMA_data$mae[["dna_variants"]])
                current_variants <- nrow(ScIGMA_data$mae[["dna_variants"]])

                # UPDATED: Extraction directe des dimensions PRE-filtrage depuis mae_raw
                initial_cells <- ncol(ScIGMA_data$mae_raw[["dna_variants"]])
                initial_variants <- nrow(ScIGMA_data$mae_raw[["dna_variants"]])

                tagList(
                    fluidRow(
                        h5("DNA variant filtering results:"),
                        column(6,
                               HTML(
                                   # UPDATED : Utilisation des variables dynamiques locales
                                   paste0("Number of cells removed: ", initial_cells - current_cells, "</br>"),
                                   paste0("Number of DNA variants removed: ", initial_variants - current_variants)
                               )
                        ),
                        column(6,
                               div(
                                   HTML(
                                       paste0("Actual number of cells: ", current_cells, "</br>"),
                                       paste0("Actual number of DNA variants : ", current_variants)
                                   )
                               ), align = "justify")
                    )
                )
            }
        })

        # [ NODE_ACCESS : RUN COMPASS ]
        # ----------------------------------------------------- _
        observeEvent(input$compass_bttn, {
            req(ScIGMA_data$mae)
        })

    })
}

## To be copied in the UI
# mod_analysis_overview_ui("analysis_right_overview_1")

## To be copied in the server
# mod_analysis_overview_server("analysis_right_overview_1", ScIGMA_data)
