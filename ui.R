library(shiny)
library(shinythemes)
library(dplyr)
library(plotly)
library(DT)
library(crosstalk)
library(readxl)

navbarPage("CPquant: Interactive quantification and exploration of CP data",
  theme = shinytheme('spacelab'),
  tabPanel("Data Processing",
           sidebarLayout(
             sidebarPanel(
               fileInput("file_1", "Choose mzXML File",
                         multiple = FALSE,
                         accept = c(".mzXML")),
               tags$hr(),
               selectInput("methods_peak", "Peak picking methods", choices = c("Method1", "Method2", "Method3")),
               tags$hr(),
               radioButtons("disp", "Display",
                            choices = c(Head = "head",
                                        All = "all"),
                            selected = "head"),
               actionButton("go_1", "Process")
             ),
             mainPanel(
               DT::DTOutput("dt_1")
             )
           )),
  tabPanel("Quantification",
           sidebarLayout(
             sidebarPanel(
               selectInput("methods_quant", "Quantification method", choices = c("Quant1", "Quant2", "Quant3"))
             ),
             mainPanel()
           )
           ),
  tabPanel("Visualization"),
  tabPanel("Statistical analysis"),
  tabPanel("Help")
)

