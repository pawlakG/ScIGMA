#' analysis_right_Protein UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_analysis_Protein_ui <- function(id) {
    ns <- NS(id)
    tagList(
        navset_card_underline(
            nav_panel("Description",
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
            nav_panel("Bi-plot",
                      fluidRow(
                          column(3,
                                 grid_card(area = "sidebar",
                                           h3("Contrôles"),
                                           selectInput("xvar", "Axe X", choices = NULL),
                                           selectInput("yvar", "Axe Y", choices = NULL),
                                           checkboxInput("logx", "Log X", FALSE),
                                           checkboxInput("logy", "Log Y", FALSE),
                                           actionButton("mk_subset", "Créer sous-échantillon ↘︎"
                                                        ))),
                          column(6,
                                 grid_card(area = "main",
                                           uiOutput("current_panel_ui"))),
                          column(3,
                                 grid_card(area = "subsets",
                                           h3("Sous-échantillons"),
                                           uiOutput("subsets_ui"))),
                          # ---- Cards / zones ----),,
                      )
            )
        )
    )
}

#' analysis_right_Protein Server Functions
#'
#' @noRd
mod_analysis_Protein_server <- function(id, ScIGMA_data){
    moduleServer(id, function(input, output, session){
        ns <- session$ns
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
