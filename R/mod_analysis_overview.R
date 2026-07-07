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
            column(
                4,
                card(
                    fileInput(ns("file_h5file"),
                        label = "1. Upload you H5 file.",
                        accept = ".h5"
                    )
                )
            ),
            column(
                4,
                card(
                    radioGroupButtons(
                        inputId = ns("file_fileType"),
                        label = "2. DNA or DNA+protein ?",
                        choices = c(
                            "DNA only" = "DNA",
                            "DNA & protein" = "DNA+protein"
                        ),
                        justified = TRUE
                    )
                )
            ),
            column(
                4,
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
            column(
                12,
                br(),
                card(
                    card_header("Summary"),
                    uiOutput(ns("overview"))
                ),
                uiOutput(ns("preprocess_card"))
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
mod_analysis_overview_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
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
        observeEvent(input$file_process, {
            filePath <- input$file_h5file$datapath
            fileType <- input$file_fileType
            req(filePath)
            req(fileType)
            sampleName <- file_path_sans_ext(input$file_h5file$name)
            message(whereami::whereami())
            # show_modal_spinner(text = "Loading data ...")
            w <- Waiter$new(
                html = spin_double_bounce(),
                color = "rgba(255, 255, 255, 0.5)" # Fond blanc semi-transparent (pas de bloc blanc opaque)
            )
            w$show()
            if (file.exists(filePath)) {
                if (file.info(filePath)$isdir) {
                    ScIGMA_data <- tryCatch(
                        ScIGMA_profile("1. Chargement des donnees (Directory)",
                            {
                                loadH5_HDF5_biocond(
                                    filepath = filePath,
                                    sample_name = sampleName,
                                    omic_type = fileType
                                )
                            },
                            filepath = filePath
                        ),
                        error = function(e) {
                            message("Error during loadH5")
                            stop(e$message)
                        }
                    )
                } else {
                    temp_scigma_obj <- tryCatch(
                        ScIGMA_profile("1. Chargement des donnees (File)",
                            {
                                loadH5_HDF5_biocond(
                                    filepath = filePath,
                                    sample_name = sampleName,
                                    omic_type = fileType
                                )
                            },
                            filepath = filePath
                        ),
                        error = function(e) {
                            message("Error during loadH5")
                            stop(e$message)
                        }
                    )
                }
            } else {
                stop("File or folder path doesn't exists\n")
            }
            ScIGMA_data$mae <- temp_scigma_obj$mae
            ScIGMA_data$mae_raw <- temp_scigma_obj$mae
            ScIGMA_data$filetype <- temp_scigma_obj$filetype

            ScIGMA_data$reset_analysis()

            init_metrics$init_number_cell <- NULL
            init_metrics$init_number_dna_variant <- NULL

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
            # print("ScIGMA_data render UI")
            # print(ScIGMA_data)
            message(whereami::whereami())

            has_mae <- !is.null(ScIGMA_data$mae)

            val_cells <- if (has_mae) nrow(SummarizedExperiment::colData(ScIGMA_data$mae)) else "No data loaded"
            val_dna <- if (has_mae) nrow(ScIGMA_data$mae[["dna_variants"]]) else "No data loaded"
            val_cnv <- if (has_mae) nrow(ScIGMA_data$mae[["amplicons"]]) else "No data loaded"

            val_prot <- if (!has_mae) {
                "No data loaded"
            } else if (!("proteins" %in% names(ScIGMA_data$mae))) {
                "No Protein data"
            } else {
                ScIGMA_object$protein.filtered <- TRUE
                nrow(ScIGMA_data$mae[["proteins"]])
            }

            fluidRow(
                column(
                    3,
                    card(
                        card_header("Number of cells"),
                        card_body(
                            p(val_cells, style = "text-align:center")
                        )
                    )
                ),
                column(
                    3,
                    card(
                        card_header("DNA variants"),
                        p(val_dna, style = "text-align:center")
                    )
                ),
                column(
                    3,
                    card(
                        card_header("Number of amplicons"),
                        p(val_cnv, style = "text-align:center")
                    )
                ),
                column(
                    3,
                    card(
                        card_header("Number of proteins"),
                        p(val_prot, style = "text-align:center")
                    )
                )
            )
        })

        # ---------------------------- #
        # Render Preprocess Card Container
        output$preprocess_card <- renderUI({
            watch("dataLoaded")
            if (!is.null(ScIGMA_data$mae)) {
                card(
                    card_header("Preprocess"),
                    uiOutput(ns("preprocess")),
                    uiOutput(ns("dnaFilterResults"))
                )
            }
        })

        # ---------------------------- #
        # Render Preprocess UI
        output$preprocess <- renderUI({
            watch("dataLoaded")
            message(whereami::whereami())
            if (!is.null(ScIGMA_data$mae)) {


                tagList(
                    fluidRow(
                        h5("Filter Cells and DNA variants:"),
                        div(
                            HTML("Adjust genotyping and variant thresholds:</br>"),
                            align = "justify"
                        )
                    ),
                    fluidRow(
                        column(3, numericInput(ns("overview_preprocess_minDp"), "Min DP", value = 10, min = 0)),
                        column(3, numericInput(ns("overview_preprocess_minGq"), "Min GQ", value = 30, min = 0)),
                        column(3, numericInput(ns("overview_preprocess_vafRef"), "Max VAF Ref (%)", value = 5, min = 0, max = 100)),
                        column(3, numericInput(ns("overview_preprocess_vafHom"), "Min VAF Hom (%)", value = 95, min = 0, max = 100))
                    ),
                    fluidRow(
                        column(4, numericInput(ns("overview_preprocess_vafHet"), "Max VAF Het (%)", value = 30, min = 0, max = 100)),
                        column(4, numericInput(ns("overview_preprocess_minCellPt"), "min.cell.pt (%)", value = 50, min = 0, max = 100)),
                        column(4, numericInput(ns("overview_preprocess_minMutCellPt"), "min.mut.cell.pt (%)", value = 1, min = 0, max = 100))
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
            }
        })

        # --------------------------------------------------------------- #
        # Filter DNA variant
        observeEvent(input$dna_variant_filtering, {
            filePath <- input$file_h5file$datapath
            overview_preprocess_minCellPt <- input$overview_preprocess_minCellPt
            overview_preprocess_minMutCellPt <- input$overview_preprocess_minMutCellPt
            overview_preprocess_minDp <- input$overview_preprocess_minDp
            overview_preprocess_minGq <- input$overview_preprocess_minGq
            overview_preprocess_vafRef <- input$overview_preprocess_vafRef
            overview_preprocess_vafHom <- input$overview_preprocess_vafHom
            overview_preprocess_vafHet <- input$overview_preprocess_vafHet
            

            req(overview_preprocess_minCellPt)
            req(overview_preprocess_minMutCellPt)
            req(overview_preprocess_minDp)
            req(overview_preprocess_minGq)
            req(overview_preprocess_vafRef)
            req(overview_preprocess_vafHom)
            req(overview_preprocess_vafHet)
            req(ScIGMA_data$mae)

            message(whereami::whereami())
            w <- Waiter$new(
                html = spin_loaders(1, color = "black"),
                color = "rgba(255, 255, 255, 0.5)" # Fond blanc semi-transparent (pas de bloc blanc opaque)
            )
            w$show()

            # ---------------------------- #
            # Reset pipeline to raw data to prevent recursive subsetting overhead
            ScIGMA_data$reset_analysis()

            # Store initial cell and DNA variant info (Utilisation du MAE_RAW)
            init_metrics$init_number_cell <- ncol(ScIGMA_data$mae)
            init_metrics$init_number_dna_variant <- nrow(ScIGMA_data$mae[["dna_variants"]])

            message(whereami::whereami())

            # ---------------------------- #
            # Execute Pipeline: Filter + Annotate + Proportions
            temp_scigma_obj <- tryCatch(
                {
                    filter_and_annotate_variants(
                        obj = ScIGMA_data,
                        paths = get_config_paths(), # Safely access config paths
                        min_dp = overview_preprocess_minDp,
                        min_gq = overview_preprocess_minGq,
                        vaf_ref = overview_preprocess_vafRef,
                        vaf_hom = overview_preprocess_vafHom,
                        vaf_het = overview_preprocess_vafHet,
                        min_cell_pt = overview_preprocess_minCellPt,
                        min_mut_cell_pt = overview_preprocess_minMutCellPt
                    )
                },
                error = function(e) {
                    remove_modal_spinner()
                    # print(e)
                    message("DNA variant filtering and annotation failed")
                    shiny::showNotification(e$message, type = "error", duration = 10)
                    req(FALSE)
                }
            )

            # ---------------------------- #
            # Update Global R6 Reference
            ScIGMA_data$mae <- temp_scigma_obj$mae

            # >> Sanitize MAE _
            ScIGMA_data$mae <- sanitize_mae_strings(ScIGMA_data$mae)


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

            if (is.null(ScIGMA_data$mae)) {
                return(NULL)
            }

            filter_status <- S4Vectors::metadata(ScIGMA_data$mae)$variant_filter

            if (is.null(filter_status) || filter_status != "filtered") {
                fluidRow(
                    h5("DNA variant filtering results:"),
                    div("Data not filtered yet", align = "center")
                )
            } else {
                # Extraction des dimensions POST-filtrage
                current_cells <- ncol(ScIGMA_data$mae[["dna_variants"]])
                current_variants <- nrow(ScIGMA_data$mae[["dna_variants"]])

                initial_cells <- ncol(ScIGMA_data$mae_raw[["dna_variants"]])
                initial_variants <- nrow(ScIGMA_data$mae_raw[["dna_variants"]])

                tagList(
                    fluidRow(
                        h5("DNA variant filtering results:"),
                        column(
                            6,
                            HTML(
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
                            ),
                            align = "justify"
                        )
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
