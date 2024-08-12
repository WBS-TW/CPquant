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
library(openxlsx)


options(shiny.maxRequestSize = 15000 * 1024^2)


# UI
ui <- shiny::navbarPage("Quantification by deconvolution from Skyline output",
                        shiny::tabPanel("Quantification Inputs", 
                                        shiny::fluidPage(shiny::sidebarLayout(
                                                shiny::sidebarPanel(
                                                        shiny::fileInput("fileInput", "Import excel file from Skyline", 
                                                                         accept = c('xlsx')),
                                                        shiny::radioButtons("quantStd",
                                                                            label = "Quantification based on chain length or mixture standards?",
                                                                            choices = c("Chain length",
                                                                                        "Mixture",
                                                                                        "Both")),
                                                        shiny::radioButtons("blankSubtraction",
                                                                            label = "Blank subtraction",
                                                                            choices = c("no blank subtraction",
                                                                                        "blank subtraction based on peak area",
                                                                                        "blank subtraction based on final concentrations")
                                                        ),
                                                        shiny::selectInput("includedCPs", "Include which CPs for quantification?",
                                                                           choices = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                                           selected = c("vSCCPs", "SCCPs", "MCCPs", "LCCPs", "vLCCPs"),
                                                                           multiple = TRUE),
                                                        shiny::actionButton('go', 'Proceed', width = "100%"),
                                                        shiny::uiOutput("defineVariables")
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
                        ),
                        shiny::tabPanel(
                                "Quantification summary",
                                # shiny::sidebarLayout(
                                #         shiny::sidebarPanel(),
                                #         shiny::mainPanel(
                                #                 DT::DTOutput("quantTable")
                                #         )
                                # )
                                shiny::fluidPage(DT::DTOutput("quantTable"))
                        )
)


################################################################################
server <- function(input, output, session) {
        
        # Create a function to process each row
        
        Skyline_output <- reactive({
                req(input$fileInput) #requires that the input is available
                df <- readxl::read_excel(input$fileInput$datapath) |> 
                        mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |> 
                        mutate(Area = as.numeric(Area)) |> 
                        mutate(Area = replace_na(Area, 0)) |> # Replace missing values in the Response_factor column with 0
                        mutate(Area = replace_na(Area, 0)) |>
                        # filter(Area >= 100) |> 
                        mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |> #--< convert to numeric #-->
                        mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) #--< convert to numeric #-->
        })
        ###START: Define input variables       
        output$defineVariables <- shiny::renderUI({
                shiny::fluidRow(
                        shiny::h4("Define variables"),
                        shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::varSelectInput(
                                        inputId = "standardAnnoColumn", #select which variable to use to define standards
                                        label = "Variable for annotating standards",
                                        data = Skyline_output(),
                                        selected = "Note"
                                )
                        ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::selectInput(
                                        inputId = "blanks", #select which variable to use to define standards
                                        label = "Define which samples are blanks",
                                        choices = unique(Skyline_output()$`Replicate Name`),
                                        multiple = TRUE
                                )
                        ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::selectInput(
                                        inputId = "removeSamples", #select if some samples will be removed from quantification 
                                        label = 'Samples to remove from quantification?',
                                        choices = unique(Skyline_output()$`Replicate Name`),
                                        multiple = TRUE
                                )
                        ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::sliderInput(
                                        inputId = "removeAreas", #remove low peak areas
                                        label = "Keep absolute peak areas above this threshold (0 means keep everything)",
                                        min = min(Skyline_output()$Area),
                                        max = max(Skyline_output()$Area),
                                        value = 0,
                                        step = 100
                                )
                        ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::sliderInput(
                                        inputId = "removeRsquared", #keep only Molecule above this rsquared, zero means keep everything 
                                        label = 'Keep only above rsquared values (0 means keep everything)',
                                        min = 0,
                                        max = 1,
                                        value = 0,
                                        step = 0.05
                                )
                        ),
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), 
                        shiny::tags$br(), shiny::tags$br(), shiny::tags$br(), shiny::tags$br(),
                        shiny::column(
                                6,
                                shiny::radioButtons("ISRS", label = "Included IS and RS?", #need to add to select which Molecule to define IS/RS
                                                    choices = c("None",
                                                                "IS only",
                                                                "RS only",
                                                                "Both IS and RS")) # select variable. TODO: set default to NOTE
                        )      
                )
        })
        ### END: Define input variables
        
        # Set reactive values from user input
        removeAreas <- eventReactive(input$go, {as.numeric(input$removeAreas)})
        removeRsquared <- eventReactive(input$go, {as.numeric(input$removeRsquared)})
        
        
        #Render raw table
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
                shiny::fluidRow(shiny::column(6, plotly::plotlyOutput("plotSummary1", width = "70vw"))
                )
        })
        
        output$plotSummary1 <- plotly::renderPlotly({
                Skyline_output() |>
                        dplyr::filter(`Isotope Label Type` == "Quan") |> 
                        plotly::plot_ly(
                                x = ~ Molecule,
                                y = ~ Area,
                                color = ~ `Sample Type`,
                                type = "box"
                        )
        })
        
        ### START: Deconvolution script
        
        shiny::observeEvent(input$go, {
                
                Skyline_output_filt <- Skyline_output() |> 
                        #filter(Area >= removeAreas())
                        mutate(Area = if_else(Area >= removeAreas(), 0, Area))
                
                
                CPs_standards <- Skyline_output_filt |> 
                        filter(`Sample Type` == "Standard",
                               Molecule != "IS",
                               Molecule != "RS",
                               `Isotope Label Type` == "Quan",
                               Note != "NA") |> 
                        group_by(!!!input$standardAnnoColumn, Molecule) |>
                        mutate(rel_int = Area/sum(Area)) |> #why ius it needed? Maybe can be removed
                        nest() |> 
                        mutate(models = map(data, ~lm(Area ~ `Analyte Concentration`, data = .x))) |> 
                        mutate(coef = map(models, coef)) |> 
                        mutate(Response_factor = map(coef, pluck("`Analyte Concentration`"))) |> 
                        mutate(intercept = map(coef, pluck("(Intercept)"))) |> 
                        mutate(rsquared = map(models, summary)) |> 
                        mutate(rsquared = map(rsquared, pluck("r.squared"))) |>
                        select(-coef) |>  # remove coef variable since it has already been plucked
                        unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
                        mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
                        mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |> 
                        #filter(rsquared >= removeRsquared()) |> 
                        mutate(Response_factor = if_else(rsquared < removeRsquared(), 0, Response_factor)) |>       
                        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
                        #filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
                        ungroup() |> 
                        group_by(!!!input$standardAnnoColumn, Chain_length) |> #grouping by the selected Note
                        mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |> 
                        ungroup()
                
                
                CPs_samples <- Skyline_output_filt |> 
                        filter(`Sample Type` == "Unknown",
                               Molecule != "IS",
                               Molecule != "RS",
                               `Isotope Label Type` == "Quan") |> 
                        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
                        #filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
                        group_by(`Replicate Name`) |> 
                        mutate(Relative_distribution = Area / sum(Area, na.rm = TRUE)) |> 
                        select(Molecule, Area, Relative_distribution) |> 
                        mutate(across(Relative_distribution, ~replace(., is.nan(.), 0))) |>    #replace NaN with zero
                        #nest() |> 
                        ungroup()
                
                
                CPs_standards_input <- CPs_standards |> 
                        select(Molecule, !!!input$standardAnnoColumn, Response_factor) |> 
                        pivot_wider(names_from = !!input$standardAnnoColumn, values_from = "Response_factor")
                
                
                CPs_samples_input <- CPs_samples |> 
                        select(Molecule, `Replicate Name`, Relative_distribution) |> 
                        pivot_wider(names_from = "Replicate Name", values_from = "Relative_distribution")
                
                
                # This step ensures that all values are corresponding to the same molecule for std and sample        
                combined <- CPs_samples_input |> 
                        right_join(CPs_standards_input, by = "Molecule")
                
                ############################################################################### DECONVOLUTION #############################################################################
                
                # Ensure combined_matrix is correctly defined as a matrix prior to the deconvolution
                combined_matrix <- CPs_standards_input  |> 
                        select(-Molecule) |> 
                        as.matrix()
                
                
                # Ensure combined_sample is correctly defined with nested data frames prior to the deconvolution
                combined_sample <- CPs_samples  |> 
                        group_by(`Replicate Name`) |> 
                        select(-Molecule, -Area) |> 
                        nest() |> 
                        ungroup()
                
                
                # Function to perform deconvolution on a single data frame
                perform_deconvolution <- function(df, combined_matrix) {
                        df_matrix <- as.matrix(df)
                        
                        print(paste("df_matrix dimensions:", dim(df_matrix)))
                        print(paste("combined_matrix dimensions:", dim(combined_matrix)))
                        
                        if (nrow(combined_matrix) != nrow(df_matrix)) {
                                stop("Dimensions of combined_matrix and df are incompatible.")
                        }
                        
                        # Reshape df_matrix if it has only one column 
                        if (ncol(df_matrix) == 1) { 
                                df_vector <- as.vector(df_matrix)
                        } else {
                                df_vector <- df_matrix
                        }
                        
                        # Perform nnls
                        deconv <- nnls(combined_matrix, df_vector)
                        
                        # Extract deconvolution results
                        deconv_coef <- deconv$x
                        deconv_resolved <- deconv$fitted.values
                        deconv_reconst <- rowSums(combined_matrix %*% deconv_coef)
                        
                        # Ensure that values are positive for chi-square test
                        if (any(deconv_resolved <= 0) || any(df_vector <= 0)) {
                                warning("Non-positive values found, skipping chi-square test")
                                chisq_result <- NULL
                        } else {
                                chisq_result <- chisq.test(deconv_resolved, p = df_vector, rescale.p = TRUE)
                        }
                        
                        return(list(
                                deconv_coef = deconv_coef,
                                deconv_resolved = deconv_resolved,
                                deconv_reconst = deconv_reconst,
                                chisq_result = chisq_result
                        ))
                }
                
                # Apply the perform_deconvolution function to each nested data frame
                Deconvolution <- combined_sample  |> 
                        mutate(result = map(data, ~ perform_deconvolution(.x, combined_matrix)))
                
                
                # Extract deconv_coef from results and create a new data frame
                deconv_coef_df <- Deconvolution  |> 
                        mutate(deconv_coef = map(result, "deconv_coef"))  |> 
                        select(`Replicate Name`, deconv_coef)  |> 
                        unnest_wider(deconv_coef, names_sep = "_")
                
                #Remove the replicate name to generate vectors:
                deconv_coef_df_matrix<- deconv_coef_df |> 
                        select(-`Replicate Name`)        
                
                ########################################################## Calculate the concentration in ng/uL ###############################################################
                
                # Create an empty list to store results
                result_list <- list()
                
                # Define a function to process each row
                process_row <- function(i, deconv_coef_df_matrix, CPs_standards_input) {
                        deconv_coef_vector <- as.numeric(deconv_coef_df_matrix[i, ])
                        combined_matrix <- as.matrix(CPs_standards_input[, -which(names(CPs_standards_input) == "Molecule")])
                        result <- sweep(combined_matrix, 2, deconv_coef_vector, `*`)
                        result_df <- as.data.frame(result)
                        colnames(result_df) <- colnames(CPs_standards_input)[-which(names(CPs_standards_input) == "Molecule")]
                        replicate_name <- deconv_coef_df$`Replicate Name`[i]
                        result_list[[replicate_name]] <- result_df
                }
                
                # Apply process_row to each row of deconv_coef_df_matrix
                result_list <- lapply(1:nrow(deconv_coef_df_matrix), process_row, deconv_coef_df_matrix, CPs_standards_input)
                result_list <- setNames(result_list, deconv_coef_df$`Replicate Name`)
                
                
                # Combine all data frames into a single data frame with 'Replicate Name'
                final_df <- do.call(rbind, Map(function(df, name) {
                        df$`Replicate Name` <- name
                        df <- df[, c("Replicate Name", setdiff(names(df), "Replicate Name"))]
                        df
                }, result_list, names(result_list)))
                
                # Add CPs_standards_input$Molecule column to final_df
                final_df$Molecule <- CPs_standards_input$Molecule
                
                # Print the final combined data frame
                print(final_df)
                
                
                #Organize the data
                final_df_tidy<-final_df|> 
                        group_by(`Replicate Name`) |> 
                        nest() |> 
                        ungroup()
                
                
                #Total sum the values for each replicate
                
                #Remove the molecule
                final_df_matrix<- final_df |> 
                        select(-Molecule) |> 
                        group_by(`Replicate Name`) |> 
                        nest()
                
                # Initialize an empty data frame to store results
                total_sums_df <- data.frame(
                        `Replicate Name` = character(),
                        `Total Sum` = numeric(),
                        stringsAsFactors = FALSE
                )
                
                # Iterate through each row of final_df_grouped
                for (i in 1:nrow(final_df_matrix)) {
                        # Extract nested data frame
                        nested_df <- final_df_matrix$data[[i]]
                        
                        # Calculate total sum for the current `Replicate Name`
                        `Replicate Name` <- final_df_matrix$`Replicate Name`[[i]]
                        total_sum <- sum(colSums(nested_df[, -1]))  # Exclude the grouping column
                        
                        # Append results to total_sums_df
                        total_sums_df <- rbind(total_sums_df, data.frame(
                                `Replicate Name` = `Replicate Name`,
                                `Total Sum` = total_sum
                        ))
                }
                
                # Print the resulting data frame
                print(total_sums_df)
                
                
                ################################################### FINAL RESULTS ####################################################################
                CPs_samples<-CPs_samples |> 
                        rename(`Replicate.Name` = `Replicate Name`)
                
                # Merge total_sums_df into CPs_samples based on Replicate Name
                Concentration <- CPs_samples  |> 
                        left_join(total_sums_df, by = "Replicate.Name")  |> 
                        mutate(Concentration = `Relative_distribution` * `Total.Sum`)
                
                print(Concentration)
                Concentration<-Concentration |> 
                        group_by(Replicate.Name) |> 
                        distinct( `Molecule`, Concentration) |> 
                        nest()
                
                # Perform operations to reorganize data
                reorganized_data <- Concentration  |> 
                        unnest(c(data)) |>  
                        distinct(`Replicate.Name`, `Molecule`, .keep_all = TRUE)  |> 
                        pivot_wider(names_from = `Molecule`, values_from = `Concentration`)
                reorganized_data <- t(reorganized_data) #transpose
                
                #Make the first row (replicate names) the column names
                colnames(reorganized_data) <- reorganized_data[1, ]
                Samples_Concentration <- reorganized_data[-1, ]
                # Convert the result back to a data frame
                Samples_Concentration <- as.data.frame(Samples_Concentration)
                Samples_Concentration <- Samples_Concentration |> 
                        mutate(Molecule = CPs_samples_input$Molecule)|> 
                        relocate(Molecule, .before = everything()) 
                
                ### END: Deconvolution script
                
                
                # Render table
                output$quantTable <- DT::renderDT({
                        DT::datatable(Samples_Concentration, 
                                      filter = "top", extensions = c("Buttons", "Scroller"),
                                      options = list(scrollY = 650,
                                                     scrollX = 500,
                                                     deferRender = TRUE,
                                                     scroller = TRUE,
                                                     buttons = list(list(extend = "excel", filename = "Samples_concentration", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "csv", filename = "Samples_concentration", title = NULL,
                                                                         exportOptions = list(
                                                                                 modifier = list(page = "all")
                                                                         )),
                                                                    list(extend = "colvis", targets = 0, visible = FALSE)),
                                                     dom = "lBfrtip",
                                                     fixedColumns = TRUE), 
                                      rownames = FALSE)
                        
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

# Run the application 
shinyApp(ui = ui, server = server)

