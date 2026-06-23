#' analysis_right_CNV UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_analysis_CNV_ui <- function(id) {
    ns <- NS(id)
    uiOutput(ns("cnv_processing"))
}

#' analysis_right_CNV Server Functions
#'
#' @noRd
mod_analysis_CNV_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        prev_plot_type <- shiny::reactiveVal("Heatmap")

        # Dynamic UI rendering
        output$cnv_processing <- renderUI({
            watch("dnaVariant_selected")
            watch("dataLoaded")
            watch("CNV_filtered")

            if (is.null(ScIGMA_data$dna.clones)) {
                tagList(
                    br(),
                    fluidRow(
                        h4(
                            "Please select DNA variant first.",
                            class = "text-muted text-center"
                        )
                    )
                )
            } else {
                cnv_ready <- ScIGMA_data$is_cnv_filtered

                # Dynamic panel focus definition under 80 characters
                active_panel <- if (cnv_ready) {
                    "Plot Configuration"
                } else {
                    "Amplicon & Cell Filters"
                }

                tagList(
                    br(),
                    bslib::accordion(
                        id = ns("cnv_accordion_main"),
                        open = active_panel,

                        # PANEL 1: AMPLICON & CELL FILTERS
                        bslib::accordion_panel(
                            title = "Amplicon & Cell Filters",
                            icon = shiny::icon("filter"),
                            fluidRow(
                                column(
                                    width = 4,
                                    numericInput(
                                        inputId = ns("cnv_ampCompleteness"),
                                        label = "Amplicon completeness (%)",
                                        value = 50
                                    )
                                ),
                                column(
                                    width = 4,
                                    numericInput(
                                        inputId = ns("cnv_ampReadDepth"),
                                        label = "Minimum amplicon read depth",
                                        value = 10
                                    )
                                )
                                # column(
                                #     width = 4,
                                #     numericInput(
                                #         inputId = ns("cnv_meanCellReadDepth"),
                                #         label = "Minimum mean cell read depth",
                                #         value = 10
                                #     )
                                # )
                            ),
                            br(),
                            fluidRow(
                                column(
                                    width = 6,
                                    shinyWidgets::materialSwitch(
                                        inputId = ns("cnv_use_compass_imputed"),
                                        label = "Use COMPASS imputed clones",
                                        value = FALSE,
                                        status = "primary"
                                    )
                                ),
                                column(
                                    width = 6,
                                    div(
                                        actionButton(
                                            inputId = ns("cnv_filter_button"),
                                            label = "Apply Filters",
                                            class = "btn-primary",
                                            icon = icon("play")
                                        ),
                                        align = "right"
                                    )
                                )
                            )
                        ),

                        # PANEL 2: PLOT CONFIGURATION
                        bslib::accordion_panel(
                            title = "Plot Configuration",
                            icon = shiny::icon("sliders"),
                            if (cnv_ready) {
                                tagList(
                                    uiOutput(ns("cnv_plot_parameters")),
                                    uiOutput(ns("cnv_plot_additionalParameters"))
                                )
                            } else {
                                fluidRow(
                                    h4(
                                        "Please filter amplicons first",
                                        class = "text-muted text-center",
                                        style = "margin-top: 20px;"
                                    )
                                )
                            }
                        ),

                        # PANEL 3: VISUALIZATION
                        bslib::accordion_panel(
                            title = "Visualization",
                            icon = shiny::icon("chart-area"),
                            if (cnv_ready) {
                                uiOutput(ns("dynamic_plot_container"))
                            } else {
                                fluidRow(
                                    h4(
                                        "Please filter amplicons first",
                                        class = "text-muted text-center",
                                        style = "margin-top: 20px;"
                                    )
                                )
                            }
                        )
                    )
                )
            }
        })

        output$cnv_placeholder <- renderUI({
            watch("CNV_filtered")
            # print("ScIGMA_data$cnv_dp_filtered")
            # print(ScIGMA_data$cnv_dp_filtered)
            if (is.null(ScIGMA_data$cnv_dp_filtered)) {
                bslib::card(
                    br(),
                    fluidRow(h4("Please filter amplicons first", class = "text-muted text-center")),
                    br()
                )
            } else {
                return(NULL)
            }
        })


        # UI for plot parameters
        observeEvent(
            {
                watch("CNV_filtered")
            },
            {
                output$cnv_plot_parameters <- renderUI({
                    if (is.null(ScIGMA_data$cnv_dp_filtered)) {
                        return(NULL)
                    } else {
                        clones_to_use <- if (!is.null(ScIGMA_data$cnv.active.clones)) ScIGMA_data$cnv.active.clones else ScIGMA_data$dna.clones
                        clone_choices <- levels(clones_to_use)[levels(clones_to_use) != "small"]

                        cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))

                        tagList(
                            fluidRow(
                                column(4, pickerInput(
                                    inputId = ns("cnv_diploidClone"),
                                    label = "Diploid clone in DNA",
                                    choices = clone_choices,
                                    options = pickerOptions(container = "body"),
                                    width = "100%"
                                )),
                                column(4, pickerInput(
                                    inputId = ns("cnv_plotType"),
                                    label = "Plot Type",
                                    choices = c("Heatmap", "Lineplot"),
                                    options = pickerOptions(container = "body"),
                                    width = "100%"
                                )),
                                column(4, pickerInput(
                                    inputId = ns("cnv_xAxis"),
                                    label = "X-axis",
                                    choices = sort_genomic_chromosomes(cnv_id_table$chrom),
                                    multiple = TRUE,
                                    options = pickerOptions(container = "body"),
                                    width = "100%"
                                ))
                            )
                        )
                    }
                })

                session$onFlushed(function() {
                    isolate(trigger("CNV_ui_cnv_plot_parameters_rendered"))
                }, once = TRUE)
            }
        )

        # Render supplementary plot parameters
        observeEvent(
            {
                watch("CNV_ui_cnv_plot_parameters_rendered")
            },
            {
                message("Rendering cnv_plot_additionalParameters")
                output$cnv_plot_additionalParameters <- renderUI({
                    if (is.null(input$cnv_plotType)) {
                        return(NULL)
                    } else {
                        clones_to_use <- if (!is.null(ScIGMA_data$cnv.active.clones)) ScIGMA_data$cnv.active.clones else ScIGMA_data$dna.clones
                        clone_choices <- levels(clones_to_use)[levels(clones_to_use) != "small"]

                        tagList(
                            if (input$cnv_plotType == "Heatmap") {
                                fluidRow(
                                    column(4, pickerInput(
                                        inputId = ns("cnv_xAxis_projection"),
                                        label = "Heatmap projection",
                                        choices = c("Position", "Genes"),
                                        options = pickerOptions(container = "body"),
                                        width = "100%"
                                    ))
                                )
                            } else {
                                fluidRow(
                                    column(4, pickerInput(
                                        inputId = ns("cnv_xAxis_projection"),
                                        label = "Lineplot projection",
                                        choices = c("Position", "Genes+amplicons"),
                                        options = pickerOptions(container = "body"),
                                        width = "100%"
                                    )),
                                    column(4, pickerInput(
                                        inputId = ns("cnv_lineplot_cluster"),
                                        label = "Select Clone",
                                        choices = clone_choices,
                                        options = pickerOptions(container = "body"),
                                        width = "100%",
                                        selected = clone_choices[1]
                                    ))
                                )
                            }
                        )
                    }
                })

                session$onFlushed(function() {
                    isolate(trigger("CNV_ui_cnv_plot_additionalParameters_rendered"))
                }, once = TRUE)
            }
        )


        # When user change a parameter value
        observeEvent(
            {
                input$cnv_filter_button
                watch("dna_clones_renamed")
            },
            ignoreInit = TRUE,
            handlerExpr = {
                message("Filtering cnv ...")
                shiny::req(ScIGMA_data$mae)
                shiny::req(input$cnv_filter_button > 0)
                if (isTRUE(input$cnv_use_compass_imputed)) {
                    if (is.null(S4Vectors::metadata(ScIGMA_data$mae)$compass)) {
                        shiny::showNotification("COMPASS inference missing. Please run COMPASS first.", type = "error")
                        shinyWidgets::updateMaterialSwitch(session, "cnv_use_compass_imputed", value = FALSE)
                        return()
                    }
                    req(ScIGMA_data$variants.filtered)

                    ht_res <- generate_dna_variant_heatmap(
                        obj = ScIGMA_data,
                        selected_variants_df = ScIGMA_data$variants.filtered,
                        heatmap_include_all_samples = FALSE,
                        use_imputed = TRUE
                    )
                    active_clones <- ht_res$clones
                } else {
                    req(ScIGMA_data$dna.clones)
                    active_clones <- ScIGMA_data$dna.clones
                }

                ScIGMA_data$cnv.active.clones <- active_clones
                # ------------------------------------------------------

                # Store values
                cnv_ampCompleteness <- input$cnv_ampCompleteness
                cnv_ampReadDepth <- input$cnv_ampReadDepth
                # cnv_meanCellReadDepth <- input$cnv_meanCellReadDepth
                cnv_meanCellReadDepth <- 10

                # Filters (Remplacement de ScIGMA_data$dna.clones par active_clones)
                filtered_data <- filter_cnv_profile(ScIGMA_data,
                    active_clones,
                    amp_completeness = cnv_ampCompleteness,
                    amp_readDepth = cnv_ampReadDepth,
                    amp_meanCellRead = cnv_meanCellReadDepth
                )

                ScIGMA_data$cnv_dp_filtered <- filtered_data
                ScIGMA_data$is_cnv_filtered <- TRUE
                trigger("CNV_filtered")

                clone_choices <- levels(active_clones)[levels(active_clones) != "small"]
                shinyWidgets::updatePickerInput(
                    session = session,
                    inputId = "cnv_diploidClone",
                    choices = clone_choices
                )
                shinyWidgets::updatePickerInput(
                    session = session,
                    inputId = "cnv_lineplot_cluster",
                    choices = clone_choices
                )
            }
        )

        # Observe event for ploidy computation
        observeEvent(input$cnv_diploidClone, {
            req(input$cnv_diploidClone)
            req(ScIGMA_data$cnv_dp_filtered)
            message("Recomputing ploidy ...")

            clones_to_use <- if (!is.null(ScIGMA_data$cnv.active.clones)) ScIGMA_data$cnv.active.clones else ScIGMA_data$dna.clones

            ploidy_data <- process_cnv_to_clonal_profile(
                ScIGMA_data$cnv_dp_filtered,
                clones_to_use,
                diploid_ref = input$cnv_diploidClone,
                exclude_clone = "small"
            )

            # Change R6 object
            ScIGMA_data$ploidy.mtx <- ploidy_data

            trigger("CNV_ploidy_computed")
        })

        observeEvent(
            {
                watch("CNV_ploidy_computed")
                input$cnv_xAxis_projection
                input$cnv_plotType
            },
            {
                req(ScIGMA_data$ploidy.mtx)

                plot_type_changed <- input$cnv_plotType != prev_plot_type()
                prev_plot_type(input$cnv_plotType)

                if (is.null(input$cnv_xAxis_projection) || input$cnv_xAxis_projection == "Position") {
                    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
                    new_choices <- sort_genomic_chromosomes(cnv_id_table$chrom)
                    label_text <- "Select Chromosome(s)"
                } else {
                    new_choices <- render_annotation_table(obj = ScIGMA_data, ploidy_data = ScIGMA_data$ploidy.mtx)$symbol |>
                        unique() |>
                        sort()
                    label_text <- "Select Gene(s)"
                }

                current_sel <- isolate(input$cnv_xAxis)

                final_selection <- if (plot_type_changed) NULL else current_sel[current_sel %in% new_choices]

                shinyWidgets::updatePickerInput(
                    session = session,
                    inputId = "cnv_xAxis",
                    label = label_text,
                    choices = new_choices,
                    selected = final_selection
                )
            },
            ignoreInit = TRUE
        )

        observeEvent(
            {
                watch("CNV_ploidy_computed")
                watch("CNV_ui_cnv_plot_additionalParameters_rendered")
                list(
                    input$cnv_plotType,
                    input$cnv_xAxis_projection,
                    input$cnv_xAxis,
                    input$cnv_xAxis_projection,
                    input$cnv_lineplot_cluster
                )
            },
            {
                req(input$cnv_plotType, ScIGMA_data$ploidy.mtx)

                cnv_plotType <- input$cnv_plotType
                cnv_xAxis <- input$cnv_xAxis

                if (cnv_plotType == "Heatmap") {
                    output$dynamic_plot_container <- renderUI({
                        plotOutput(ns("static_plot"))
                    })

                    cnv_heatmap_type <- input$cnv_xAxis_projection
                    show_genes <- (!is.null(cnv_heatmap_type) && cnv_heatmap_type != "Position")

                    # --- Filtrage Heatmap ---
                    filtered_ploidy <- ScIGMA_data$ploidy.mtx
                    if (!is.null(cnv_xAxis) && length(cnv_xAxis) > 0) {
                        mat_data_tmp <- t(filtered_ploidy)
                        cnv_id_table_tmp <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
                        genome_v_tmp <- S4Vectors::metadata(ScIGMA_data$mae)$genome_version
                        if (is.null(genome_v_tmp)) genome_v_tmp <- "hg19"

                        tmp_var <- cnv_id_table_tmp |> dplyr::filter(dna_id %in% colnames(mat_data_tmp))

                        is_chr_focus <- any(grepl("^chr([0-9]+|[XYM])$", cnv_xAxis, ignore.case = TRUE))

                        if (!is_chr_focus) {
                            tmp_annot <- annotate_genomic_regions(region_data = tmp_var, build = genome_v_tmp)
                            valid_ids <- tmp_annot$dna_id[tmp_annot$symbol %in% cnv_xAxis]
                            show_genes <- TRUE
                        } else {
                            tmp_var$chr_lit <- paste0("chr", tmp_var$chrom)
                            valid_ids <- tmp_var$dna_id[tmp_var$chr_lit %in% cnv_xAxis]
                            show_genes <- FALSE
                        }
                        filtered_ploidy <- t(mat_data_tmp[, colnames(mat_data_tmp) %in% valid_ids, drop = FALSE])
                    }

                    shiny::validate(
                        shiny::need(
                            ncol(filtered_ploidy) >= 2 && nrow(filtered_ploidy) >= 2,
                            "Selection too narrow. Requires at least 2 cells and 2 amplicons."
                        ),
                        shiny::need(
                            sd(as.vector(filtered_ploidy), na.rm = TRUE) > 0,
                            "Zero variance detected. All cells have the exact same copy number in this region."
                        )
                    )

                    tmp_ht <- plot_cnv_heatmap(
                        obj = ScIGMA_data,
                        ploidy_data = filtered_ploidy,
                        display_gene = show_genes
                    )
                    output$static_plot <- renderPlot(tmp_ht)
                } else {
                    # Lineplot
                    output$dynamic_plot_container <- renderUI({
                        plotlyOutput(ns("interactive_plot"), height = "600px")
                    })

                    req(input$cnv_xAxis_projection, input$cnv_lineplot_cluster)
                    cnv_lineplot_type <- input$cnv_xAxis_projection
                    cnv_lineplot_cluster <- input$cnv_lineplot_cluster

                    mat_data <- t(ScIGMA_data$ploidy.mtx)
                    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
                    genome_v <- S4Vectors::metadata(ScIGMA_data$mae)$genome_version
                    if (is.null(genome_v)) genome_v <- "hg19"

                    tmp_var_table <- cnv_id_table |>
                        dplyr::filter(dna_id %in% colnames(mat_data)) |>
                        dplyr::arrange(as.numeric(chrom), as.numeric(start_pos)) |>
                        dplyr::mutate(chr_lit = paste0("chr", chrom))

                    # --- Filtrage Lineplot ---
                    if (!is.null(cnv_xAxis) && length(cnv_xAxis) > 0) {
                        is_chr <- any(grepl("^chr", cnv_xAxis, ignore.case = TRUE))

                        if (is_chr) {
                            tmp_var_table <- tmp_var_table |> dplyr::filter(chr_lit %in% cnv_xAxis)
                        } else {
                            tmp_annot <- annotate_genomic_regions(region_data = tmp_var_table, build = genome_v)
                            valid_ids <- tmp_annot$dna_id[tmp_annot$symbol %in% cnv_xAxis]
                            tmp_var_table <- tmp_var_table |> dplyr::filter(dna_id %in% valid_ids)
                        }
                    }

                    shiny::validate(
                        shiny::need(
                            nrow(tmp_var_table) >= 2,
                            "The selected region contains fewer than 2 amplicons. Please expand your selection to view the Lineplot."
                        )
                    )

                    mat_data <- mat_data[, tmp_var_table$dna_id, drop = FALSE]

                    tmp_split_table <- tmp_var_table[match(colnames(mat_data), tmp_var_table$dna_id), ]
                    sorted_gen_levels <- sort_genomic_chromosomes(tmp_split_table$chrom)

                    tmp_split_vec <- annotate_genomic_regions(region_data = tmp_split_table, build = genome_v)

                    gene_annotation <- data.frame(
                        "Gene" = tmp_split_vec$symbol,
                        "Chromosome" = tmp_split_vec$chrom,
                        "Probe" = tmp_split_vec$dna_id,
                        "Chrom_pos" = factor(tmp_split_vec$chr_lit,
                            levels = unique(sort_genomic_chromosomes(tmp_split_vec$chrom))
                        ),
                        "Chrom_start" = tmp_split_vec$start_pos
                    )

                    tmp_plot <- plot_cnv_genome(
                        cnv_matrix = mat_data,
                        sub_indices = cnv_lineplot_cluster,
                        gene_annotation = gene_annotation,
                        lineplot_type = cnv_lineplot_type
                    )

                    output$interactive_plot <- renderPlotly(tmp_plot)
                }
            }
        )
    })
}

## To be copied in the UI
# mod_analysis_CNV_ui("analysis_CNV_1")

## To be copied in the server
# mod_analysis_CNV_server("analysis_CNV_1")
