# TODO
# New tab Skyline: add checkbox to include monoisotopic formula even if it is not above threshold
# Rename columns according to Fernandes et al. 2023
# Add a `Precursor m/z` column but needs to be calculated using the getAdduct() function since getSkyline gives the monoisotopic and not adduct m/z
# Add [M+Br]-

# New tab: Add a user input (csv) file to see if the chosen quantifier and qualifier ions have interference from other fragment ions... 
# ...user should first export to excel and then filter only those used for quan/qual. Need a new column to indicate this?

# use isowrap instead to get more accurate envelope profile isotopic fine structure
# Add BCP and BCO: https://pubs.acs.org/doi/10.1021/acs.est.2c03576

# Reactive log: https://shiny.rstudio.com/articles/debugging.html

library(shiny)
library(shinythemes)
library(DT)
library(tidyverse)
library(readxl)
library(plotly)
library(crosstalk)
library(enviPat)
library(markdown)

data("isotopes")
source("./R/getAdduct.R")
source("./R/getSkyline.R")



#--------------------------------UI function----------------------------------#

ui <- shiny::navbarPage(
        "Chlorinated paraffins/olefins ion explorer",
        theme = shinythemes::shinytheme('spacelab'),
        shiny::tabPanel("Initial settings",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("Cmin", "C atoms min (allowed 3-40)", value = 9, min = 3, max = 40),
                                        shiny::numericInput("Cmax", "C atoms max (allowed 4-40)", value = 30, min = 4, max = 40),
                                        shiny::numericInput("Clmin", "Cl atoms min (allowed 1-15))", value = 3, min = 1, max = 15),
                                        shiny::numericInput("Clmax", "Cl atoms max (allowed 1-15)", value = 15, min = 1, max = 15),
                                        shiny::br(),
                                        selectInput("Adducts", "Add adducts/fragments",
                                                    choices = c("[CP-Cl]-", 
                                                                "[CP-H]-",
                                                                "[CP-HCl]-", 
                                                                #"[CP-Cl-HCl]-", #Needs to verify that regex extraction in getAdduct can get these ions before adding these
                                                                #"[CP-2Cl-HCl]-", 
                                                                "[CP+Cl]-", 
                                                                "[CO-Cl]-", 
                                                                "[CO-HCl]-", 
                                                                "[CO-H]-",
                                                                "[CO+Cl]-",
                                                                "[CP+Br]-"
                                                                #"[CP-Cl-HCl]+", 
                                                                #"[CP-Cl-2HCl]+", 
                                                                #"[CP-Cl-3HCl]+", 
                                                                #"[CP-Cl-4HCl]+"),
                                                    ),
                                                    selected = "[CP-Cl]-",
                                                    multiple = TRUE,
                                                    selectize = TRUE,
                                                    width = NULL,
                                                    size = NULL),
                                        shiny::numericInput("threshold", "Isotope rel ab threshold (in %)", value = 5, min = 0, max = 99),
                                        shiny::actionButton("go1", "Submit", width = "100%"),
                                        width = 3),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table", width = "100%")
                                )
                        )
                        )),
        shiny::tabPanel("Interfering ions",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        shiny::numericInput("MSresolution", "MS Resolution", value = 60000, min = 1000, max = 3000000),
                                        shiny::actionButton("go2", "Calculate", width = "100%"),
                                        width = 2
                                ),
                                shiny::mainPanel(
                                        plotly::plotlyOutput("Plotly"),
                                        plotly::plotlyOutput("Plotly2"),
                                        DT::dataTableOutput("Table2", width = "100%")
                                        
                                )
                        )
                        )),
        shiny::tabPanel("Skyline",
                        shiny::fluidPage(shiny::sidebarLayout(
                                shiny::sidebarPanel(
                                        #shiny::numericInput("MSresolution2", "MS Resolution", value = 60000, min = 100, max = 3000000),
                                        shiny::actionButton("go3", "Transition List", width = "100%"),
                                        shiny::checkboxInput("MonoisoAsQuant", "Mark monoisotopic peak as quant ion", FALSE),
                                        width = 4
                                ),
                                shiny::mainPanel(
                                        DT::dataTableOutput("Table3", width = "100%")
                                        
                                )
                        )
                        )),
        shiny::tabPanel(
                "Instructions",
                shiny::sidebarLayout(
                        shiny::sidebarPanel(shiny::h3("Manual"),
                                            width = 3),
                        shiny::mainPanel(
                                shiny::includeMarkdown("instructions.md")
                        )
                )
        )
)



#---------------------------Shiny Server function------------------------------#



server = function(input, output, session) {
        
        # Set reactive values from user input
        C <- eventReactive(input$go1, {as.integer(input$Cmin:input$Cmax)})
        Cl <- eventReactive(input$go1, {as.integer(input$Clmin:input$Clmax)})
        threshold <- eventReactive(input$go1, {as.integer(input$threshold)})
        selectedAdducts <- eventReactive(input$go1, {as.character(input$Adducts)})
        MSresolution <- eventReactive(input$go2, {as.integer(input$MSresolution)})
        
        
        #----Outputs_Start
        
        CP_allions_glob <- eventReactive(input$go1, {
                
                # Create a Progress bar object
                progress <- shiny::Progress$new()
                
                # Make sure it closes when we exit this reactive, even if there's an error
                on.exit(progress$close())
                progress$set(message = "Calculating", value = 0)
                
                Adducts <- as.character(selectedAdducts())
                
                # function to get adducts or fragments
                CP_allions <- list()
                for (i in seq_along(Adducts)) {
                        progress$inc(1/length(Adducts), detail = paste0("Adduct: ", Adducts[i], " . Please wait.."))
                        input <- getAdduct(adduct_ions = Adducts[i], C = C(), Cl = Cl(), threshold = threshold())
                        CP_allions <- rbind(CP_allions, input)
                        
                }
                return(CP_allions)
        })
        
        # go1: Calculate the isotopes from initial settings tab
        shiny::observeEvent(input$go1, {
                output$Table <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download the visible rows of the table, but this will also give warning about large tables
                        # Show data
                        DT::datatable(CP_allions_glob(), 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        # go1 end
        
        
        
        # go2: Calculates the interfering ions tab
        shiny::observeEvent(input$go2, {
                
                CP_allions_compl2 <- CP_allions_glob() %>%
                        arrange(`m/z`) %>%
                        mutate(difflag = round(abs(`m/z` - lag(`m/z`, default = first(`m/z`))),6)) %>%
                        mutate(difflead = round(abs(`m/z` - lead(`m/z`, default = last(`m/z`))), 6)) %>%
                        mutate(reslag = round(`m/z`/difflag, 0)) %>%
                        mutate(reslead = round(`m/z`/difflead, 0)) %>%
                        mutate(interference = case_when(
                                difflag == 0 | difflead == 0 ~ TRUE, # need to keep this true to make same mass ions TRUE
                                reslag >= as.integer(MSresolution()) | reslead >= as.integer(MSresolution()) ~ TRUE,
                                reslag < as.integer(MSresolution()) & reslead < as.integer(MSresolution()) ~ FALSE
                        )
                        )
                # change first and last row to false since their lead/lag is zero
                CP_allions_compl2$interference[1] <- FALSE
                CP_allions_compl2$interference[length(CP_allions_compl2$interference)] <- FALSE
                
                # Output scatterplot: #Cl vs #C
                output$Plotly <- plotly::renderPlotly(
                        p <- CP_allions_compl2 %>% plot_ly(
                                x = ~ (`12C`+`13C`), 
                                y = ~(`35Cl`+`37Cl`),
                                type = "scatter",
                                mode = "markers",
                                color = ~interference,
                                hoverinfo = "text",
                                hovertext = paste("Parent Formula:", CP_allions_compl2$Parent_Formula,
                                                  '<br>',
                                                  "Adduct/Fragment ion:", CP_allions_compl2$Fragment,
                                                  '<br>',
                                                  "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                  '<br>',
                                                  "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                             "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`))
                        )
                        %>% 
                                plotly::layout(xaxis = list(title = "Number of carbons (12C+13C)"),
                                               yaxis = list(title = "Number of chlorines (35Cl+37Cl)"),
                                               legend=list(title=list(text='<b> Interference at MS res? </b>')))
                )
                
                # Output the interference bar plot: Rel_ab vs m/z
                output$Plotly2 <- plotly::renderPlotly(
                        p <- CP_allions_compl2 %>% plot_ly(
                                x = ~`m/z`, 
                                y = ~Rel_ab,
                                type = "bar",
                                color = ~interference,
                                #text = ~Fragment,
                                hoverinfo = "text",
                                hovertext = paste("Parent Formula:", CP_allions_compl2$Parent_Formula,
                                                  '<br>',
                                                  "Adduct/Fragment ion:", CP_allions_compl2$Fragment,
                                                  '<br>',
                                                  "Ion Formula:", CP_allions_compl2$Adduct_Formula,
                                                  '<br>',
                                                  "Adduct isotopes:", paste0("[12C]:", CP_allions_compl2$`12C`, "  [13C]:", CP_allions_compl2$`13C`, 
                                                                             "  [35Cl]:", CP_allions_compl2$`35Cl`, "  [37Cl]:", CP_allions_compl2$`37Cl`),
                                                  '<br>',
                                                  "m/z:", CP_allions_compl2$`m/z`,
                                                  '<br>',
                                                  "m/z diff (prev and next):", CP_allions_compl2$difflag, "&", CP_allions_compl2$difflead,
                                                  '<br>',
                                                  "Resolution needed (prev and next):", CP_allions_compl2$reslag, "&", CP_allions_compl2$reslead)
                        )
                        %>% 
                                plotly::layout(legend=list(title=list(text='<b> Interference at MS res? </b>')))
                )
                
                output$Table2 <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download only part of rows
                        # Show data
                        DT::datatable(CP_allions_compl2, 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        # go2 end
        
        # go3: Skyline tab
        
        CP_allions_skyline <- eventReactive(input$go3, {
                
                # Create a Progress bar object
                progress <- shiny::Progress$new()
                
                # Make sure it closes when we exit this reactive, even if there's an error
                on.exit(progress$close())
                progress$set(message = "Calculating", value = 0)
                
                Adducts <- as.character(selectedAdducts())
                
                # function to get Skyline adducts or fragments
                CP_allions <- list()
                for (i in seq_along(Adducts)) {
                        progress$inc(1/length(Adducts), detail = paste0("Adduct: ", Adducts[i], " . Please wait.."))
                        input <- getSkyline(adduct_ions = Adducts[i], C = C(), Cl = Cl(), threshold = threshold())
                        CP_allions <- rbind(CP_allions, input)
                        
                }
                return(CP_allions)
        })
        shiny::observeEvent(input$go3, {
                
                CP_allions_skyline <- CP_allions_skyline() %>%
                        mutate(`Molecule List Name` = case_when(str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ paste0("CP-C", `12C`),
                                                                str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ paste0("CO-C", `12C`))) %>%
                        rename(`Molecule Name` = Parent_Formula) %>%
                        mutate(`Molecular Formula` = case_when(
                                `37Cl` == 0 ~ paste0("C", `12C`, "H", `1H`, "Cl", `35Cl`),
                                `37Cl` > 0 ~ paste0("C", `12C`, "H", `1H`, "Cl", `35Cl`, "Cl'", `37Cl`)
                        )) %>%
                        mutate(`Precursor Adduct` = str_replace(Fragment, "\\].*", "]")) %>% 
                        mutate(`Precursor Adduct` = str_replace(`Precursor Adduct`, "(.+?(?=\\-))|(.+?(?=\\+))", "[M")) %>%
                        rename(`Precursor Charge` = Charge) %>%
                        add_column(`Explicit Retention Time` = NA) %>%
                        add_column(`Explicit Retention Time Window` = NA) %>%
                        group_by(`Molecule Name`) |> 
                        mutate(Note = ifelse(Rel_ab == 100, "Quan", "Qual")) |> # choose the highest rel_ab ion as quan ion and the rest will be qual
                        ungroup() |> 
                        # mutate(Note = case_when(
                        #         `12C` < 10 & str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ "vSCCP",
                        #         `12C` < 10 & str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ "vSCCO",
                        #         `12C` >= 10 & `12C` < 14 & str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ "SCCPs",
                        #         `12C` >= 10 & `12C` < 14 & str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ "SCCOs",
                        #         `12C` >= 14 & `12C` < 18 & str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ "MCCPs",
                        #         `12C` >= 14 & `12C` < 18 & str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ "MCCOs",
                        #         `12C` >= 18 & str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ "LCCPs",
                        #         `12C` >= 18 & str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ "LCCOs")) %>%
                        select(`Molecule List Name`, 
                               `Molecule Name`, 
                               #Fragment, 
                               `Molecular Formula`, 
                               `Precursor Adduct`, 
                               `Precursor Charge`, 
                               `Explicit Retention Time`, 
                               `Explicit Retention Time Window`, 
                               Note)
                
                
                
                output$Table3 <- DT::renderDT(server=FALSE,{ #need to keep server = FALSE otherwise excel download only part of rows
                        # Show data
                        DT::datatable(CP_allions_skyline, 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", filename = "Transition List", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", title = "Transition List",
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                })
        })
        
        
        
        # go3 end
        
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
