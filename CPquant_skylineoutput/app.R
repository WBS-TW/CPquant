# CPquant: Shiny CPs Quantification for Skyline Output
#
#
#
library(shiny)
library(tidyverse)
library(readxl)
library(plotly)
library(DT)
library(nnls)

# UI
ui <- shiny::navbarPage("Quantification by deconvolution from Skyline output",
        shiny::tabPanel("Quantification Inputs", 
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::fileInput("fileInput", "Import excel file from Skyline", 
                                                         accept = c('xlsx')),
                                        shiny::radioButtons("quantstd",
                                                            label = "Quantification based on chain length or mixture standards?",
                                                            choices = c("Chain length",
                                                                        "Mixture",
                                                                        "Both")),
                                        shiny::radioButtons("blanks",
                                                            label = "Blank subtraction",
                                                            choices = c("no blank subtraction",
                                                                        "blank subtraction based on peak area",
                                                                        "blank subtraction based on final concentrations")
                                                            ),
                                        shiny::selectInput("includedCPs", "Include which CPs for quantification?",
                                                           choices = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                           selected = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                           multiple = TRUE),
                                        shiny::uiOutput("defineVariables"),
                                        shiny::actionButton('go', 'Proceed', width = "100%")
                                        ),
                                shiny::mainPanel(
                                        DT::DTOutput("table1")
                                        )
                                )
                                )),
        shiny::tabPanel(
                "Input summary",
                shiny::sidebarLayout(
                        shiny::sidebarPanel(),
                        shiny::mainPanel(
                                shiny::uiOutput("inputSummary")
                                )
                        )
                )
        )
        

# Define server logic required to draw a histogram
server <- function(input, output) {
        
        Skyline_output <- reactive({
                req(input$fileInput) #requires that the input is available
                df <- readxl::read_excel(input$fileInput$datapath)
                
                })
        
        output$defineVariables <- shiny::renderUI({
                shiny::fluidRow(
                        shiny::h4("Define variables"),
                        shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::selectInput(
                                        inputId = "standardType", #select which variable to use to define standards
                                        label = "Variable for annotating standards",
                                        choices = names(Skyline_output()) # select variable. TODO: set default to NOTE
                                        )
                                ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::selectInput(
                                        inputId = "removeSamples", #select if some samples will be removed. allows 
                                        label = 'Samples to remove from quantification?',
                                        choices = unique(Skyline_output()$`Replicate Name`),
                                        multiple = TRUE
                                        )
                                ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::radioButtons("ISRS", label = "Included IS and RS?", #select which Molecule to define IS/RS
                                                    choices = c("None",
                                                                "IS only",
                                                                "RS only",
                                                                "Both IS and RS")) # select variable. TODO: set default to NOTE
                                )      
                        )
                })
        
        # Render table
        output$table1 <- DT::renderDT({
                DT::datatable(Skyline_output(),
                              options = list(
                                      paging = TRUE,
                                      pageLength = 50
                              )
                )
                
        })
        
        # Render summary statistics and plots of raw input BEFORE quantification
        output$inputSummary <- shiny::renderUI({
                shiny::fluidRow(shiny::column(6, plotly::plotlyOutput("plotSummary1"))
                                )
                })
        
        output$plotSummary1 <- plotly::renderPlotly({
                Skyline_output() |>
                        dplyr::filter(`Isotope Label Type` == "Quan") |> 
                                plotly::plot_ly(
                                         x = ~ Molecule,
                                         y = ~ `Normalized Area`,
                                         color = ~ `Sample Type`,
                                         type = "box"
                                )
                })
        
}

# Run the application 
shinyApp(ui = ui, server = server)
