# TODO
# Add textbox to input all adducts and iteratively calculate these
        # (maybe "CP-Cl" to indicate chlorinated paraffins CPs, "CO-Cl" to indicate chlorinated olefins)
# Add also calculations for chlorinated olefins: https://pubs.acs.org/doi/10.1021/acs.analchem.7b00331; https://www.sciencedirect.com/science/article/pii/S002196732030114X?via%3Dihub
# Add filter on C and Cl

# Instructions:
# Add that [M+Cl-HCl]- can be written as [M-H]-

# Reactive log: https://shiny.rstudio.com/articles/debugging.html

library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)
library(readxl)
#library(rcdk)
library(enviPat)

data("isotopes")
source("./R/getAdduct.R")

#-------Functions-------#


# Additional data on adducts, but not needed yet
# data("adducts")
# Cl_adducts <- readxl::read_xlsx("./data/CP_adducts.xlsx") %>%
#         mutate(Formula_add = as.character(Formula_add))
# adducts <- adducts %>% bind_rows(Cl_adducts)


#--------------------------------UI function----------------------------------#

ui <- shiny::navbarPage(
        "Chlorinated Paraffins interference generator",
        theme = shinythemes::shinytheme('spacelab'),
        shiny::tabPanel("Initial settings",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("Cmin", "C atoms min", value = 9, min = 1, max = 30),
                                        shiny::numericInput("Cmax", "C atoms max", value = 30, min = 1, max = 30),
                                        shiny::numericInput("Clmin", "Cl atoms min", value = 3, min = 1, max = 30),
                                        shiny::numericInput("Clmax", "Cl atoms max", value = 15, min = 1, max = 30),
                                        shiny::br(),
                                        selectInput("Adducts", "Add adducts/fragments",
                                                choices = c("[CP-Cl]-", "[CP-HCl]-", "[CO-Cl]-", "[CO-HCl]-"),
                                                selected = "[CP-Cl]-",
                                                multiple = TRUE,
                                                selectize = TRUE,
                                                width = NULL,
                                                size = NULL),
                                        shiny::numericInput("threshold", "Threshold for isotopic calculations", value = 10, min = 1, max = 95),
                                        shiny::actionButton("go", "Submit", width = "100%")
                                        ),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table")
                                        
                                )
                        )
                        )),
        shiny::tabPanel(
                "Instructions",
                shiny::sidebarLayout(
                        shiny::sidebarPanel(shiny::h3("Table of content"),
                                            shiny::h4("File input"),
                                            shiny::h4("Equation"),
                                            width = 3),
                        shiny::mainPanel(shiny::h3("Add instructions here!")
                                                )
                )
        )
)



#---------------------------Shiny Server function------------------------------#



server = function(input, output, session) {
        
        # function to get isotopic patterns for all CPs. Limit the threshold to 10%, -1 charge. data("isotopes") needs to be loaded first
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, chemforms = x, threshold = 10, plotit = FALSE, charge = -1)}
        
        # Set values from user input
        C <- eventReactive(input$go, {as.integer(input$Cmin:input$Cmax)})
        Cl <- eventReactive(input$go, {as.integer(input$Clmin:input$Clmax)})
        threshold <- eventReactive(input$go, {as.integer(input$threshold)})
        selectedAdducts <- eventReactive(input$go, {as.character(input$Adducts)})
        
        # Outputs
  
        shiny::observeEvent(input$go, {
                
                # Create a Progress object
                progress <- shiny::Progress$new()
                # Make sure it closes when we exit this reactive, even if there's an error
                on.exit(progress$close())
                progress$set(message = "Calculating", value = 0)
                
                Adducts <- as.character(selectedAdducts())
                
                # function to get adducts or fragments
                CP_allions_compl <- list()
                for (i in seq_along(Adducts)) {
                        progress$inc(1/length(Adducts), detail = paste0("Adduct: ", Adducts[i], " . Please wait.."))
                        input <- getAdduct(adduct_ions = Adducts[i], C = C(), Cl = Cl(), threshold = threshold())
                        CP_allions_compl <- rbind(CP_allions_compl, input)
                        }
              
                #output$Table <- DT::renderDataTable(CP_allions_compl)
                
                output$Table <- DT::renderDT(server=FALSE,{
                        # Show data
                        datatable(CP_allions_compl, 
                                  rownames = FALSE,
                                  extensions = 'Buttons', 
                                  options = list(
                                          paging = FALSE,
                                          dom = 'Bfrtip',
                                          buttons = list(list(extend = "excel", title = NULL),
                                                         list(extend = "csv", title = NULL))
                                          )
                                  )
                        })
               
                })
        


        
        
        
       
        
# Close the app when the session ends
if(!interactive()) {
        session$onSessionEnded(function() {
                stopApp()
                q("no")
        })
}
        
}

shiny::shinyApp(ui, server)
