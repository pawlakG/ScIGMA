#' analysis_left UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @import shinybusy
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_analysis_left_ui <- function(id) {
    ns <- NS(id)
    ## Add a spinner
    add_busy_spinner(spin = "fading-circle", color = "#112446")
    tagList(
        fileInput(ns("file_h5file"),
            label = "1. Upload you H5 file.",
            accept = ".h5"
        ),
        textInput(ns("file_name"),
            label = "2. Enter a name for your assay",
            value = ""
        ),
        radioGroupButtons(
            inputId = ns("file_fileType"),
            label = "3. DNA or DNA+protein ?",
            choices = c(
                "DNA only" = "DNA",
                "DNA & protein" = "DNA+protein"
            ),
            justified = TRUE
        ),
        h6(HTML("4. Process")),
        actionBttn(
            inputId = ns("file_process"),
            label = "Load file",
            color = "primary",
            style = "stretch",
            icon = icon("magnifying-glass-chart"),
            block = TRUE
        )
    )
}

#' analysis_left Server Functions
#'
#' @noRd
#' @importFrom gargoyle  init
mod_analysis_left_server <- function(id, ScIGMA_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        # Init watcher
        init("preprocessFile")

        # Uploaded file
        message(whereami::whereami())

        observeEvent(input$file_process, {
            filePath <- input$file_h5file$datapath
            sampleName <- input$file_name
            fileType <- input$file_fileType
            req(filePath)
            req(sampleName)
            req(fileType)
            message(whereami::whereami())
            show_modal_spinner(text = "Loading data ...")
            if (file.exists(filePath)) {
                if (file.info(filePath)$isdir) {
                    ScIGMA_data$data <- tryCatch(
                        loadH5_dir_HDF5(
                            dir = filePath,
                            feature_policy = "intersect",
                            omic.type = fileType
                        ),
                        error = function(e) {
                            message("Error during loadH5")
                            stop(e$message)
                        }
                    )
                } else {
                    ScIGMA_data$data <- tryCatch(
                        loadH5_HDF5(
                            filepath = filePath,
                            sample.name = sampleName,
                            omic.type = fileType
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
            remove_modal_spinner()
            message(whereami::whereami())
            trigger("dataLoaded")
        })
    })
}
