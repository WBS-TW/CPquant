# TODO
# Add textbox to input all adducts and iteratively calculate these
        # (maybe "CP-Cl" to indicate chlorinated paraffins CPs, "CO-Cl" to indicate chlorinated olefins)
# Add also calculations for chlorinated olefins: https://pubs.acs.org/doi/10.1021/acs.analchem.7b00331; https://www.sciencedirect.com/science/article/pii/S002196732030114X?via%3Dihub
# Add filter on C and Cl

# Instructions:
# Add that [M+Cl-HCl]- can be written as [M-H]-


library(shiny)
library(shinythemes)
library(tidyverse)
library(readxl)
#library(rcdk)
library(enviPat)

data("isotopes")


#-------Functions-------#

# function to get isotopic patterns for all CPs. Limit the threshold to 1%, neutral form. data("isotopes") needs to be loaded first
getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, chemforms = x, threshold = 1, plotit = FALSE, charge = -1)}

getAdduct <- function(adduct_ions) {
        
        ion_modes <- str_extract(adduct_ions, "(?<=\\]).*") # lookbehind assertion to extract ion mode
        fragment_ions <- str_extract(adduct_ions, "(?<=.{3}).+?(?=\\])") # extract after the 2nd character and before ]
        group <- str_extract(adduct_ions, "[^\\[].{1}") # positive lookbehind for [)
        
        if (group == "CP") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl < 15) %>% # limit chlorination level CHECK WITH LCCPs!!
                        mutate(H = 2*C+2-Cl) %>% # add H atoms
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Formula, C, H, Cl) # move Formula to first column
        } else if (group == "CO") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl < 15) %>% # limit chlorination level CHECK WITH LCCPs!!
                        mutate(H = 2*C-Cl) %>% # add H atoms
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Formula, C, H, Cl) # move Formula to first column
        }  else {
                print("Input not correct, only CP or CO is allowed")
        }
        
        # check chem_forms
        # if (any(check_chemform(isotopes = isotopes, chemforms = data$Formula)$warning == TRUE)) {print("Warning: incorrect formula")} else {"All correct"}
        
        # adding ion modes to the data frame to be inserted to isopattern
        if (ion_modes == "-") {
                data <- data %>%
                        mutate(Charge = as.integer(-1))
        }else if (ion_modes == "+") {
                data <- data %>%
                        mutate(Charge = as.integer(1))
        }
        
        # generate input data for envipat based on fragment_ions
        
        if (fragment_ions == "-Cl") { # Generate fragments M-Cl for each homolog formula  
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Fragment = adduct_ions) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Charge, Fragment, Formula, C, H, Cl)
        } else if (fragment_ions == "-HCl") { # Generate fragments M-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Fragment = adduct_ions) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(H = H-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Charge, Fragment, Formula, C, H, Cl)
        } else if (fragment_ions == "+Cl") { # Generate fragments M+Cl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Fragment = adduct_ions) %>%
                        mutate(Cl = Cl+1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Charge, Fragment, Formula, C, H, Cl)
        }
        
        # All formula for M-Cl
        CP_allions <- list()
        data_ls <- list()
        
        
        # function to get isotopic patterns for all CPs. Limit the threshold to 10%, neutral form. data("isotopes") needs to be loaded first
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, chemforms = x, threshold = 10, plotit = FALSE, charge = Charge)}
        
        for (j in seq_along(data$Formula)) {
                Formula <- data$Formula[j]
                Parent <- data$Parent[j]
                Charge <- data$Charge[j]
                dat <- getisotopes(x = data$Formula[j])
                dat <- as.data.frame(dat[[1]])
                dat <- dat %>% 
                        mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`)) %>%
                        mutate(Parent_Formula = Parent) %>%
                        mutate(Frag_Formula =  Formula) %>%
                        mutate(Charge = Charge) %>%
                        mutate(Isotope = case_when(
                                `13C` + `37Cl`*2 == 0 ~ "",
                                `13C` + `37Cl`*2 == 1 ~ "+1",
                                `13C` + `37Cl`*2 == 2 ~ "+2",
                                `13C` + `37Cl`*2 == 3 ~ "+3",
                                `13C` + `37Cl`*2 == 4 ~ "+4",
                                `13C` + `37Cl`*2 == 5 ~ "+5",
                                `13C` + `37Cl`*2 == 6 ~ "+6",
                                `13C` + `37Cl`*2 == 7 ~ "+7",
                                `13C` + `37Cl`*2 == 8 ~ "+8",
                                `13C` + `37Cl`*2 == 9 ~ "+9",
                                `13C` + `37Cl`*2 == 10 ~ "+10",
                                `13C` + `37Cl`*2 == 11 ~ "+11",
                                `13C` + `37Cl`*2 == 12 ~ "+12",
                                `13C` + `37Cl`*2 == 13 ~ "+13",
                                `13C` + `37Cl`*2 == 14 ~ "+14",
                                `13C` + `37Cl`*2 == 15 ~ "+15",
                                `13C` + `37Cl`*2 == 16 ~ "+16",
                                `13C` + `37Cl`*2 == 17 ~ "+17",
                                `13C` + `37Cl`*2 == 18 ~ "+18",
                                `13C` + `37Cl`*2 == 19 ~ "+19",
                                `13C` + `37Cl`*2 == 20 ~ "+20")) %>%
                        mutate(Fragment = paste0(adduct_ions, Isotope)) %>%
                        select(Parent_Formula, Charge, Fragment, Frag_Formula, Isotope, Isotope_Formula, `m/z`, abundance, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
                data_ls[[j]] <- dat
        }
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}

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
                                        shiny::numericInput("Clmin", "Cl atoms min", value = 1, min = 1, max = 30),
                                        shiny::numericInput("Clmax", "Cl atoms max", value = 30, min = 1, max = 30),
                                        shiny::br(),
                                        #shiny::textInput("Adducts", "Add adducts/fragments", value = "[CP-Cl]-"),
                                        selectInput("Adducts", "Add adducts/fragments",
                                                choices = c("[CP-Cl]-", "[CP-HCl-", "[CO-Cl]-", "[CO-HCl]-"),
                                                selected = "[CP-Cl]-",
                                                multiple = TRUE,
                                                selectize = TRUE,
                                                width = NULL,
                                                size = NULL),
                                        shiny::numericInput("threshold", "Threshold for isotopic calculations", value = 10, min = 1, max = 95),
                                        shiny::actionButton('go', 'Submit', width = "100%")
                                ),
                                shiny::mainPanel(
                                        shiny::textOutput("viewAtoms"),
                                        shiny::textOutput("viewAdducts")
                                        
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
        

        C <- eventReactive(input$go, {as.integer(input$Cmin:input$Cmax)})
        
        Cl <- eventReactive(input$go, {as.integer(input$Clmin:input$Clmax)})
        
        selectedAdducts <- eventReactive(input$go, {input$Adducts})
        output$viewAtoms <- renderText({C()})
        output$viewAdducts <- renderText({selectedAdducts()})
        
        shiny::observeEvent(input$go, {
                # function to get adducts or fragments
                CP_allions_compl <- reactive({list()})
                # for (i in seq_along(selectedAdducts())) {
                #         print(i)
                #         
                # }
                
                for (i in seq_along(selectedAdducts())) {
                        input <- reactive({getAdduct(adduct_ions = selectedAdducts()[i])})
                        CP_allions_compl <- reactive({rbind(CP_allions_compl, input())})
                        }
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
