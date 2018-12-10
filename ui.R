library(shiny)
library(shinythemes)
library(dplyr)
library(plotly)
library(DT)
library(crosstalk)

navbarPage("CPquant: Interactive quantification and exploration of CP data",
  theme = shinytheme('spacelab'),
  tabPanel("Data",
           sidebarLayout(
             sidebarPanel(
               fileInput("file_1", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               tags$hr(),
               checkboxInput("header", "Header", TRUE),
               radioButtons("sep", "Separator",
                            choices = c(Comma = ",",
                                        Semicolon = ";",
                                        Tab = "\t"),
                            selected = ","),
               radioButtons("quote", "Quote",
                            choices = c(None = "",
                                        "Double Quote" = '"',
                                        "Single Quote" = "'"),
                            selected = '"'),
               tags$hr(),
               radioButtons("disp", "Display",
                            choices = c(Head = "head",
                                        All = "all"),
                            selected = "head"),
               actionButton("go_1", "Upload")
             ),
             mainPanel(
               DT::DTOutput("dt_1")
             )
           ),
  tabPanel("Quantification"),
  tabPanel("Visualize"),
  tabPanel("Help")
)
)
