#' analysis_right_overview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_analysis_right_overview_ui <- function(id) {
    ns <- NS(id)
    ## Add a spinner
    add_busy_spinner(spin = "fading-circle", color = "#112446")
    tagList(
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
}

#' analysis_right_overview Server Functions
#'
#' @noRd
mod_analysis_right_overview_server <- function(id, lensObject){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
        message(whereami::whereami())

        # ---------------------------- #
        # Render Summary UI
        output$overview <- renderUI({
            watch("dataLoaded")
            message(whereami::whereami())
            fluidRow(
                column(3,
                       card(
                           card_header("Number of cells"),
                           card_body(
                               p(ifelse(is.null(lensObject$data),yes = "No data loaded", no = length(lensObject$data@cell.ids)),
                                 style="text-align:center")
                           ))
                       # summaryBox(title = "Number of cells", value = ifelse(is.null(lensObject$data),yes = "No data loaded", no = length(lensObject$data@cell.ids)), icon = icon("credit-card"))
                ),
                column(3,
                       card(
                           card_header("DNA variants"),
                           p(ifelse(is.null(lensObject$data),yes = "No data loaded", no = length(lensObject$data@variants)),
                             style="text-align:center")
                       )
                ),
                column(3,
                       card(
                           card_header("number of CNVs"),
                           p(ifelse(is.null(lensObject$data),yes = "No data loaded", no = length(lensObject$data@amps)),
                             style="text-align:center")
                       )
                ),
                column(3,
                       card(
                           card_header("Number of proteins"),
                           p(ifelse(is.null(lensObject$data),yes = "No data loaded", no = length(lensObject$data@proteins)),
                             style="text-align:center")
                       )
                )
            )
        })

        # ---------------------------- #
        # Render Preprocess UI
        output$preprocess <- renderUI({
            watch("dataLoaded")
            message(whereami::whereami())
            if(!is.null(lensObject$data)){
                tagList(
                    fluidRow(
                        h5("Filter DNA variants:"),

                        div(
                            HTML("DNA Variant filtering step has two parameters : </br>
              <b>min.cell.pt</b>	Minimum threshold for cell percentage that has valid variant call (GT = 0, 1 or 2) after
              applying the filter. Minimum threshold for cell percentage that has valid variant call (GT = 0, 1 or 2)
              after applying the filter. The default value is 50%. This means for one variant; we need at least 50% of
              cells have a valid variant call. When to change: If the variant of interest is in a high GC content
              region, then PCR amplification is hard. In such cases, you may choose to decrease the percent to 30 or 40
              so that your interested variant could come through the filter.</br>
              </br>
              <b>min.mut.cell.pt</b>	Minimum threshold for cell percentage that has mutated genotype (GT = 1 or 2) after
              applying the filter. The default is 1, corresponds to 1%. This filter is used to remove false positives.
              When to change: If you know the variant is rare in the data, then you could try lower threshold to try to
              keep the variant in your dataset."),
                            align = "justify"
                        )
                    ),
                    fluidRow(
                        column(6,
                               sliderTextInput(
                                   inputId = ns("overview_preprocess_minCellPt"),
                                   label = "min.cell.pt",
                                   choices = seq(1, 100,1),
                                   grid = TRUE,
                                   selected = 50
                               )
                        ),
                        column(6,
                               sliderTextInput(
                                   inputId = ns("overview_preprocess_minMutCellPt"),
                                   label = "min.mut.cell.pt",
                                   choices = seq(1, 30, 1),
                                   grid = TRUE,
                                   selected = 1
                               )
                        )
                    ),
                    fluidRow(
                        actionBttn(
                            inputId = ns("overview_process"),
                            label = "Filter DNA variants",
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
        observeEvent(input$overview_process,{
            filePath <- input$file_h5file$datapath
            overview_preprocess_minCellPt <- input$overview_preprocess_minCellPt
            overview_preprocess_minMutCellPt <- input$overview_preprocess_minMutCellPt
            req(overview_preprocess_minCellPt)
            req(overview_preprocess_minMutCellPt)
            message(whereami::whereami())
            show_modal_spinner()
            # ---------------------------- #
            # Filter variants
            lensObject$data <- tryCatch(
                filterVariant(optima.obj = lensObject$data,
                              min.cell.pt = overview_preprocess_minCellPt,
                              min.mut.cell.pt = overview_preprocess_minMutCellPt),
                error = function(e){
                    remove_modal_spinner()

                    message("Error during DNA variant filtering")
                    stop(e$message)
                })

            message(whereami::whereami())
            trigger("dnaVariant_filtered")
            # ---------------------------- #
            # Annotate variants
            print("length(lensObject$data$variants)")
            print(length(lensObject$data$variants))
            if (length(lensObject$data$variants) == 0){
                # CONDITION TO HANDLE
            } else {
                lensObject$variantAnnotation <- tryCatch(
                    annotateVariant(lensObject$data$variants)
                , error = function(e){
                    remove_modal_spinner()
                    message(warning("Error during variant annotation: "),
                            warning(e$message))
                })
                remove_modal_spinner()
            }
            remove_modal_spinner()
        })



        # --------------------------------------------------------------- #
        # Render UI after DNA variant filtering
        output$dnaFilterResults <- renderUI({
            watch("dnaVariant_filtered")
            message(whereami::whereami())
            message(whereami::whereami())
            if (is.null(lensObject$data)){
                fluidRow(
                    h5("DNA variant filtering results:"),
                    div("Data not filtered yet", align ="center")
                )
            } else if (lensObject$data@variant.filter != "filtered") {
                fluidRow(
                    h5("DNA variant filtering results:"),
                    div("Data not filtered yet", align ="center")
                )
            } else {
                tagList(
                    fluidRow(
                        h5("DNA variant filtering results:"),
                        column(6,
                               HTML(
                                   paste0("Number of cells removed: ", lensObject$initNumberCell - length(lensObject$data@cell.ids), "</br>"),
                                   paste0("Number of DNA variants removed: ", lensObject$initNumberDNA_variant - length(lensObject$data@variants))
                               )
                        ),
                        column(6,
                               div(
                                   HTML(
                                       paste0("Actual number of cells: ", length(lensObject$data@cell.ids), "</br>"),
                                       paste0("Actual number of DNA variants : ", length(lensObject$data@variants))
                                   )
                               ), align = "justify")
                    )
                )
            }
        })

    })
}

## To be copied in the UI
# mod_analysis_right_overview_ui("analysis_right_overview_1")

## To be copied in the server
# mod_analysis_right_overview_server("analysis_right_overview_1", lensObject)
