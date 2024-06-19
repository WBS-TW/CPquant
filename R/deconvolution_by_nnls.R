
## Setup automatic nnls deconvolution process from Skyline excel output
library(tidyverse)
library(readxl)
library(nnls)
#### TODO ####
#filter(rsquared > 0.7) #cannot do this as some are still false positive and need to replace formula later with 0 for nnls 

# NEED TO ADD using case_when or elseif..
# if (sum(!is.na(filtered_dataA$`Normalized Area`)) == 0 || sum(!is.na(filtered_dataA$`Analyte Concentration`)) == 0) {
#       cat("No valid cases for fitting the model for molecule:", molecule_name, "\n")
# OR add a filter to remove those rsquared that are very low (perhaps below 0.7?)  

##############

# reading data from excel
Skyline_output <- read_excel("./data/OrbitrapDust.xlsx") |>
        mutate(`Analyte Concentration` = as.numeric(`Analyte Concentration`)) |> 
        mutate(`Normalized Area` = as.numeric(`Normalized Area`)) |> 
        mutate(`Normalized Area` = replace_na(`Normalized Area`, 0)) |> # Replace missing values in the Response_factor column with 0
        mutate(Area = replace_na(Area, 0)) |>
        mutate(RatioQuanToQual = as.numeric(RatioQuanToQual)) |> #--< convert to numeric #-->
        mutate(RatioQualToQuan = as.numeric(RatioQualToQuan)) #--< convert to numeric #-->


#This is currently filtered for C10-C13 only to compare with the Perkons script. Will remove later or add as an argument in function
CPs_standards <- Skyline_output |> 
        filter(`Sample Type` == "Standard",
               Molecule != "IS",
               Molecule != "RS",
               `Isotope Label Type` == "Quan",
               Note != "NA") |> 
        group_by(Note, Molecule) |>
        mutate(rel_int = `Normalized Area`/sum(`Normalized Area`)) |> #why is it needed? Maybe can be removed
        nest() |> 
        mutate(models = map(data, ~lm(`Normalized Area` ~ `Analyte Concentration`, data = .x))) |> 
        mutate(coef = map(models, coef)) |> 
        mutate(Response_factor = map(coef, pluck("`Analyte Concentration`"))) |> 
        mutate(intercept = map(coef, pluck("(Intercept)"))) |> 
        mutate(rsquared = map(models, summary)) |> 
        mutate(rsquared = map(rsquared, pluck("r.squared"))) |>
        select(-coef) |>  # remove coef variable since it has already been plucked
        unnest(c(Response_factor, intercept, rsquared)) |>  #removing the list type for these variables
        mutate(Response_factor = if_else(Response_factor < 0, 0, Response_factor)) |> # replace negative RF with 0
        mutate(rsquared = ifelse(is.nan(rsquared), 0, rsquared)) |> 
        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
        filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
        ungroup() |> 
        group_by(Note, Chain_length) |> 
        mutate(Sum_response_factor_chainlength = sum(Response_factor, na.rm = TRUE)) |> 
        ungroup() |> 
        group_by(Note) |> #group by standard 
        mutate(rel_response_factor = Response_factor/sum(Response_factor)) |> #rel RF by standard group
        ungroup()



CPs_samples <- Skyline_output |> 
        filter(`Sample Type` == "Unknown",
               Molecule != "IS",
               Molecule != "RS",
               `Isotope Label Type` == "Quan") |> 
        mutate(Chain_length = paste0("C", str_extract(Molecule, "(?<=C)[^H]+"))) |> 
        filter(Chain_length == "C10" | Chain_length == "C11" | Chain_length == "C12" | Chain_length == "C13") |> #this will be remove later or added as arg in fn
        group_by(`Replicate Name`) |> 
        mutate(Relative_distribution = `Normalized Area` / sum(`Normalized Area`, na.rm = TRUE)) |> 
        ungroup() |> 
        mutate(across(Relative_distribution, ~replace(., is.nan(.), 0))) #replace NaN with zero

CPs_samples_individual <- CPs_samples |> 
        filter(`Replicate Name` == "NIST_R1") |> 
        select(Molecule, Relative_distribution)


CPs_standards_input <- CPs_standards |>
        select(Molecule, Note, Response_factor) |>
        pivot_wider(names_from = "Note", values_from = "Response_factor")

# This one uses the relative response factor instead of absolute response factor for comparison. Results should be the same?
# CPs_standards_input <- CPs_standards |> 
#         select(Molecule, Note, rel_response_factor) |> 
#         pivot_wider(names_from = "Note", values_from = "rel_response_factor")

# This step ensures that all values are corresponding to the same molecule for std and sample        
combined <- CPs_samples_individual |> 
        right_join(CPs_standards_input, by = "Molecule")

combined_matrix <- combined |> 
        select(-Molecule, -Relative_distribution) |> 
        as.matrix()

combined_sample <- combined$Relative_distribution

# Using the non-negative least squares, nnls, package to deconvolute the distribution
deconv <- nnls(combined_matrix, combined_sample)

deconv

deconv_coef <- deconv$x

deconv_resolved <- deconv[["fitted"]]

deconv_reconst <- rowSums(combined_matrix %*% deconv_coef)

chisq.test(deconv_resolved, p = combined_sample, rescale.p = TRUE)


par(mfrow = c(2,1))
barplot(combined_sample, main = "Sample")
barplot(deconv_reconst, main = "Reconstructed")


