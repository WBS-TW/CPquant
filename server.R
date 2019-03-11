library(shiny)
library(shinythemes)
library(dplyr)
library(plotly)
library(DT)
library(crosstalk)

options(shiny.maxRequestSize=2000*1024^2) 

shinyServer(function(input, output, session) {
  
  observeEvent(input$go_1, {
    output$dt_1 <- renderDT({
      req(input$file_1)
      tryCatch(
        {
        df <- read.csv(input$file_1$datapath,
                       header = input$header,
                       sep = input$sep,
                       quote = input$quote)
        },
        error = function(e) {
          stop(safeError(e))
          }
        )
      if(input$disp == "head") {
        return(head(df))
        }
      })
    })
})
  
