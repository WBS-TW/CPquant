# TODO


# Instructions:
# Add that [M+Cl-HCl]- can be written as [M-H]-
# Chlorinated paraffins are written as [CP] and chlorinated olefins as [CO]

# Reactive log: https://shiny.rstudio.com/articles/debugging.html

library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)
library(readxl)
library(plotly)
library(enviPat)

data("isotopes")
source("./R/getAdduct.R")



#--------------------------------UI function----------------------------------#

ui <- shiny::navbarPage(
        "Chlorinated Paraffins interference generator",
        theme = shinythemes::shinytheme('spacelab'),
        shiny::tabPanel("Initial settings",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("Cmin", "C atoms min (3-30)", value = 9, min = 3, max = 30),
                                        shiny::numericInput("Cmax", "C atoms max", value = 30, min = 4, max = 30),
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
                                        shiny::numericInput("threshold", "Isotope rel ab threshold (in %)", value = 10, min = 1, max = 95),
                                        shiny::actionButton("go1", "Submit", width = "100%")
                                        ),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table", width = "100%")
                                        )
                                )
                        )),
        shiny::tabPanel("Interfering ions",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("MSresolution", "MS Resolution", value = 60000, min = 1000, max = 3000000),
                                        shiny::actionButton("go2", "Calculate interfering ions", width = "100%")
                                        ),
                                shiny::mainPanel(
                                        plotly::plotlyOutput("Plotly"),
                                        DT::dataTableOutput("Table2", width = "100%")
                                        
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
        
        # Set values from user input
        C <- eventReactive(input$go1, {as.integer(input$Cmin:input$Cmax)})
        Cl <- eventReactive(input$go1, {as.integer(input$Clmin:input$Clmax)})
        threshold <- eventReactive(input$go1, {as.integer(input$threshold)})
        selectedAdducts <- eventReactive(input$go1, {as.character(input$Adducts)})
        MSresolution <- eventReactive(input$go2, {as.integer(input$MSresolution)})
        

#----Outputs_Start

        CP_allions_compl_glob <- eventReactive(input$go1, {
                
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
                        return(CP_allions_compl)
                }
        })
              
                
                shiny::observeEvent(input$go1, {
                output$Table <- DT::renderDT(server=TRUE,{
                        # Show data
                        DT::datatable(CP_allions_compl_glob(), 
                                  filter = "top", extensions = c("Buttons", "Scroller"),
                                  options = list(scrollY = 650,
                                                 scrollX = 500,
                                                 deferRender = TRUE,
                                                 scroller = TRUE,
                                                 # paging = TRUE,
                                                 # pageLength = 25,
                                                 buttons = list(list(extend = "excel", title = NULL),
                                                                list(extend = "colvis", targets = 0, visible = FALSE)),
                                                 dom = "lBfrtip",
                                                 fixedColumns = TRUE), 
                                  rownames = FALSE)
                        })
                })

                

        
        shiny::observeEvent(input$go2, {
                
                #CP_allions_compl2 <- as.data.frame(CP_allions_compl_glob())

                CP_allions_compl2 <- CP_allions_compl_glob() %>%
                        arrange(`m/z`) %>%
                        mutate(difflag = round(abs(`m/z` - lag(`m/z`, default = first(`m/z`))),6)) %>%
                        mutate(difflead = round(abs(`m/z` - lead(`m/z`, default = last(`m/z`))), 6)) %>%
                        mutate(reslag = round(`m/z`/difflag, 0)) %>%
                        mutate(reslead = round(`m/z`/difflead, 0)) %>%
                        mutate(interference = case_when(
                                difflag == 0 | difflead == 0 ~ FALSE,
                                reslag >= as.integer(MSresolution()) | reslead >= as.integer(MSresolution()) ~ TRUE,
                                reslag < as.integer(MSresolution()) & reslead < as.integer(MSresolution()) ~ FALSE
                                )
                               )

                
                output$Plotly <- plotly::renderPlotly(
                        p <- CP_allions_compl2 %>% plot_ly(
                                x = ~Parent_Formula, 
                                y = ~`m/z`,
                                #size = ~abundance,
                                type = "scatter",
                                mode = "markers",
                                color = ~interference)
                        %>% 
                        plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
                )
                
                output$Table2 <- DT::renderDT(server=TRUE,{
                        # Show data
                        DT::datatable(CP_allions_compl2, 
                                  filter = "top", extensions = c("Buttons", "Scroller"),
                                  options = list(scrollY = 650,
                                                 scrollX = 500,
                                                 deferRender = TRUE,
                                                 scroller = TRUE,
                                                 # paging = TRUE,
                                                 # pageLength = 25,
                                                 buttons = list(list(extend = "excel", title = NULL),
                                                                list(extend = "colvis", targets = 0, visible = FALSE)),
                                                 dom = "lBfrtip",
                                                 fixedColumns = TRUE), 
                                  rownames = FALSE)
                })
                })
        
#----Outputs_End

        
        
        
       
        
# Close the app when the session ends
if(!interactive()) {
        session$onSessionEnded(function() {
                stopApp()
                q("no")
        })
}
        
}

shiny::shinyApp(ui, server)
