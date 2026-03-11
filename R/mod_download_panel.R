#' download_panel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_download_panel_ui <- function(id) {
  ns <- NS(id)
  tagList(

  )
}

#' download_panel Server Functions
#'
#' @noRd
mod_download_panel_server <- function(id, ScIGMA_data){
  moduleServer(id, function(input, output, session){
    ns <- session$ns

  })
}

## To be copied in the UI
# mod_download_panel_ui("download_panel_1")

## To be copied in the server
# mod_download_panel_server("download_panel_1")
