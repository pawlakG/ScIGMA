#' welcome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_welcome_ui <- function(id) {
    ns <- NS(id)
    tagList(
        card(
            full_screen = TRUE,
            card_header(HTML("Welcome to ScIGMA: Single-cell Integrated Genomic & Multi-omics Analyzer")),
            # NEW: Insert graphical abstract with responsive constraints
            card_body(
                div(
                    style = "padding: 4% 8% 0% 8%; width: 100%;",
                    tags$img(
                        src = "www/ScIGMA_graphicalAbstract.png",
                        style = "width: 100%; height: auto; display: block; margin: 0 auto;",
                        alt = "ScIGMA Graphical Abstract"
                    )
                )
            )
        )
    )
}

#' welcome Server Functions
#'
#' @noRd
mod_welcome_server <- function(id){
    moduleServer(id, function(input, output, session){
        ns <- session$ns

    })
}

## To be copied in the UI
# mod_welcome_ui("welcome_1")

## To be copied in the server
# mod_welcome_server("welcome_1")
