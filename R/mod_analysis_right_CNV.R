#' analysis_right_CNV UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_analysis_right_CNV_ui <- function(id) {
  ns <- NS(id)
  tagList(
 
  )
}
    
#' analysis_right_CNV Server Functions
#'
#' @noRd 
mod_analysis_right_CNV_server <- function(id){
  moduleServer(id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_analysis_right_CNV_ui("analysis_right_CNV_1")
    
## To be copied in the server
# mod_analysis_right_CNV_server("analysis_right_CNV_1")
