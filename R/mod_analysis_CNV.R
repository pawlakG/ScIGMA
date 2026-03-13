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
mod_analysis_CNV_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
        # Dynamic UI rendering
        output$cnv_processing <- renderUI({
            watch("dnaVariant_selected")
            if(is.null(ScIGMA_data$dna.clones)){
                tagList(
                    br(),
                    fluidRow(h3("Please select DNA variant first."))
                )

            } else {
                tagList(
                    br(),
                    card(fluidRow(div(h2("Amplicon and cell filters"))),
                         fluidRow(
                             column(3,
                                    numericInput(ns("cnv_ampCompleteness"), label = h3("Amplicon completeness"), value = 50)),
                             column(3,
                                    numericInput(ns("cnv_ampReadDepth"), label = h3("Amplicon read depth"), value = 10)),
                             column(3,
                                    numericInput(ns("cnv_meanCellReadDepth"), label = h3("Mean cell read depth"), value = 10)),
                             column(3,
                                    div(actionButton(ns("cnv_filter_button"), "Filter",
                                                     class = "btn-primary"), align = "center"))
                         )
                    ),
                    uiOutput(ns("cnv_plot_parameters")),
                    uiOutput(ns("cnv_plot_additionalParameters")),
                    # plotOutput(ns("cnv_plot"))
                    uiOutput(ns("dynamic_plot_container"))
                )
            }
        })

        # UI for plot parameters
        observeEvent({
            watch("CNV_filtered")
        },{
            output$cnv_plot_parameters <- renderUI({
                if (is.null(ScIGMA_data$cnv_dp_filtered)){
                    card(br(), fluidRow(h3("Please filter CNV data first")))
                } else {
                    # Extraction à la volée pour les choix d'interface
                    cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))

                    card(
                        fluidRow(
                            column(4, pickerInput(
                                inputId = ns("cnv_diploidClone"),
                                label = h2("Diploid clone in DNA"),
                                choices = levels(ScIGMA_data$dna.clones)[levels(ScIGMA_data$dna.clones) != "small"],
                                options = pickerOptions(container = "body"),
                                width = "100%"
                            )),
                            column(4,pickerInput(
                                inputId = ns("cnv_plotType"),
                                label = h2("Plot"),
                                choices = c("Heatmap", "Lineplot"),
                                options = pickerOptions(container = "body"),
                                width = "100%"
                            )),
                            column(4,pickerInput(
                                inputId = ns("cnv_xAxis"),
                                label = "X-axis",
                                choices = sort_genomic_chromosomes(cnv_id_table$chrom), # UPDATED
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
        })

        # Render supplementary plot parameters
        observeEvent({
            watch("CNV_ui_cnv_plot_parameters_rendered")
        },{
            message("Rendering cnv_plot_additionalParameters")
            output$cnv_plot_additionalParameters <- renderUI({
                if(is.null(input$cnv_plotType)){
                    card(
                        br(),
                        fluidRow(h3("Please select a plot type first"))
                    )
                } else {
                    card(
                        if (input$cnv_plotType == "Heatmap"){
                            fluidRow(
                                pickerInput(
                                    inputId = ns("cnv_heatmap_type"),
                                    label = h2("Heatmap type"),
                                    choices = c("Position", "Genes"),
                                    options = pickerOptions(container = "body"),
                                    width = "100%"
                                )
                            )
                        } else {
                            fluidRow(
                                column(6, pickerInput(
                                    inputId = ns("cnv_lineplot_type"),
                                    label = h2("Lineplot projection"),
                                    choices = c("Position", "Genes+amplicons"),
                                    options = pickerOptions(container = "body"),
                                    width = "100%"
                                )),
                                column(6,pickerInput(
                                    inputId = ns("cnv_lineplot_cluster"),
                                    label = h2("Clone"),
                                    choices = levels(ScIGMA_data$dna.clones)[levels(ScIGMA_data$dna.clones) != "small"],
                                    options = pickerOptions(container = "body"),
                                    width = "100%",
                                    selected = levels(ScIGMA_data$dna.clones)[levels(ScIGMA_data$dna.clones) != "small"][1]
                                ))
                            )
                        }
                    )
                }
            })
            session$onFlushed(function() {
                isolate(trigger("CNV_ui_cnv_plot_additionalParameters_rendered"))
            }, once = TRUE)
        })



        # When user change a parameter value
        observeEvent(eventExpr = input$cnv_filter_button,
                     ignoreInit = TRUE,
                     handlerExpr = {
                         message("Filtering cnv ...")
                         # Store values
                         cnv_ampCompleteness <- input$cnv_ampCompleteness
                         cnv_ampReadDepth <- input$cnv_ampReadDepth
                         cnv_meanCellReadDepth <- input$cnv_meanCellReadDepth
                         cnv_diploidClone <- input$cnv_diploidClone

                         # Filters
                         filtered_data <- filter_cnv_profile(ScIGMA_data,
                                                             ScIGMA_data$dna.clones,
                                                             amp_completeness = cnv_ampCompleteness,
                                                             amp_readDepth = cnv_ampReadDepth,
                                                             amp_meanCellRead = cnv_meanCellReadDepth)

                         ScIGMA_data$cnv_dp_filtered <- filtered_data
                         message("Filtering done")
                         trigger("CNV_filtered")
                     })

        # Observe event for ploidy computation
        observeEvent(input$cnv_diploidClone,
                     {
                         req(input$cnv_diploidClone)
                         message("Recomputing ploidy ...")
                         # To put into another observeEvent
                         ploidy_data <- process_cnv_to_clonal_profile(
                             ScIGMA_data$cnv_dp_filtered,
                             ScIGMA_data$dna.clones,
                             diploid_ref = input$cnv_diploidClone,
                             exclude_clone = "small"
                         )

                         # Change R6 object
                         ScIGMA_data$ploidy.mtx <- ploidy_data

                         trigger("CNV_ploidy_computed")

                     })

        # Observe event after ploidy recomputation or any input change
        observeEvent({
            watch("CNV_ploidy_computed")
            watch("CNV_ui_cnv_plot_additionalParameters_rendered")
            list(input$cnv_plotType,
                 input$cnv_heatmap_type,
                 input$cnv_xAxis,
                 input$cnv_lineplot_type,
                 input$cnv_lineplot_cluster)
        },
        {
            req(input$cnv_plotType,
                input$cnv_heatmap_type
                # input$cnv_xAxis
            )

            cnv_plotType <- input$cnv_plotType
            cnv_heatmap_type <- input$cnv_heatmap_type
            cnv_xAxis <- input$cnv_xAxis

            print("cnv_heatmap_type")
            print(cnv_heatmap_type)
            print("cnv_plotType")
            print(cnv_plotType)

            if (cnv_heatmap_type == "Position"){
                show_genes <- FALSE
                cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
                updatePickerInput(
                    session = session,
                    inputId = "cnv_xAxis",
                    choices = sort_genomic_chromosomes(cnv_id_table$chrom), # UPDATED
                    selected = NULL
                )
            } else {
                show_genes <- TRUE
                tmp_choices <- render_annotation_table(obj = ScIGMA_data,
                                                       ploidy_data = ScIGMA_data$ploidy.mtx)$symbol |>
                    unique() |>
                    sort()
                updatePickerInput(
                    session = session,
                    inputId = "cnv_xAxis",
                    choices = tmp_choices,
                    selected = NULL # Reset selection to avoid carry-over from previous state
                )
            }

            if (cnv_plotType == "Heatmap") {
                output$dynamic_plot_container <- renderUI({plotOutput(ns("static_plot"))})

                print("ScIGMA_data$ploidy.mtx after cnv_plotType Heatmap")
                print(head(ScIGMA_data$ploidy.mtx))

                tmp_ht <- plot_cnv_heatmap(obj = ScIGMA_data,
                                           ploidy_data = ScIGMA_data$ploidy.mtx,
                                           display_gene = show_genes)
                output$static_plot <- renderPlot(tmp_ht)
            } else {
                # Lineplot
                output$dynamic_plot_container <- renderUI({plotlyOutput(ns("interactive_plot"), height = "600px")})
                ## Get inputs
                print("cnv_lineplot_type")
                print(input$cnv_lineplot_type)
                print("cnv_lineplot_cluster")
                print(input$cnv_lineplot_cluster)
                req(input$cnv_lineplot_type)
                req(input$cnv_lineplot_cluster)
                cnv_lineplot_type <- input$cnv_lineplot_type
                cnv_lineplot_cluster <- input$cnv_lineplot_cluster
                ## compute

                # --- Lineplot Block ---
                mat_data <- t(ScIGMA_data$ploidy.mtx)
                cnv_id_table <- as.data.frame(SummarizedExperiment::rowData(ScIGMA_data$mae[["amplicons"]]))
                genome_v <- S4Vectors::metadata(ScIGMA_data$mae)$genome_version
                if (is.null(genome_v)) genome_v <- "hg19"

                tmp_var_table <- cnv_id_table |>
                    dplyr::filter(dna_id %in% colnames(mat_data)) |>
                    dplyr::arrange(as.numeric(chrom), as.numeric(start_pos)) |>
                    dplyr::mutate(chr_lit = paste0("chr", chrom))

                mat_data <- mat_data[, tmp_var_table$dna_id, drop = FALSE]

                tmp_split_table <- tmp_var_table[match(colnames(mat_data), tmp_var_table$dna_id), ]
                sorted_gen_levels <- sort_genomic_chromosomes(tmp_split_table$chrom)

                tmp_split_vec <- annotate_genomic_regions(region_data = tmp_split_table, build = genome_v)

                gene_annotation = data.frame('Gene' = tmp_split_vec$symbol,
                                             'Chromosome' = tmp_split_vec$chrom,
                                             'Probe' = tmp_split_vec$dna_id,
                                             'Chrom_pos' = factor(tmp_split_vec$chr_lit,
                                                                  levels = unique(sort_genomic_chromosomes(tmp_split_vec$chrom))),
                                             'Chrom_start' = tmp_split_vec$start_pos)

                tmp_plot <- plot_cnv_genome(cnv_matrix = mat_data,
                                            sub_indices = cnv_lineplot_cluster,
                                            gene_annotation = gene_annotation,
                                            lineplot_type = cnv_lineplot_type)

                output$interactive_plot <- renderPlotly(tmp_plot)
            }
        })
    })
}

## To be copied in the UI
# mod_analysis_CNV_ui("analysis_CNV_1")

## To be copied in the server
# mod_analysis_CNV_server("analysis_CNV_1")
