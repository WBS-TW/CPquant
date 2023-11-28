

getSkyline <- function(adduct_ions, C, Cl, threshold) {
        
        ion_modes <- str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
        fragment_ions <- str_extract(adduct_ions, "(?<=.{3}).+?(?=\\])") # extract after the 2nd character and before ]
        group <- str_extract(adduct_ions, "[^\\[].{1}") # Using positive lookbehind for [)
        
        if (group == "CP") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl < 15) %>% # limit chlorine atoms. CHECK WITH LCCPs!!
                        mutate(H = 2*C+2-Cl) %>% # add H atoms
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Formula, C, H, Cl) # move Formula to first column
        } else if (group == "CO") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl < 15) %>% # limit chlorine atoms. CHECK WITH LCCPs!!
                        mutate(H = 2*C-Cl) %>% # add H atoms
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Formula, C, H, Cl) # move Formula to first column
        }  else {
                print("Input not correct, only CP or CO is allowed")
        }
        
        # check chem_forms
        # if (any(check_chemform(isotopes = isotopes, chemforms = data$Formula)$warning == TRUE)) {print("Warning: incorrect formula")} else {"All correct"}
        
        # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
        if (ion_modes == "-") {
                data <- data %>%
                        mutate(Charge = as.integer(-1))
        }else if (ion_modes == "+") {
                data <- data %>%
                        mutate(Charge = as.integer(1))
        }
        
        # generate input data for envipat based on fragment_ions
        
         # Generate fragments M-Cl for each homolog formula  
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Cl_perc = case_when(group == "CP" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 2),
                                                   group == "CO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 2))) %>%
                        mutate(Fragment = adduct_ions) %>%
                        mutate(Cl = Cl) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Cl_perc, Charge, Fragment, Formula, C, H, Cl)
         
        
        # Remove formula without Cl after adduct formations
        data <- data %>%
                filter(Cl > 0)
        
        # All ion formulas
        CP_allions <- list()
        data_ls <- list()
        
        
        # function to get isotopic patterns for all CPs. Limit the threshold to 10%, neutral form. data("isotopes") needs to be loaded first
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, 
                                                        chemforms = x, 
                                                        threshold = threshold, 
                                                        emass = 0.00054857990924,
                                                        plotit = FALSE, 
                                                        charge = Charge)}
        
        for (j in seq_along(data$Formula)) {
                Formula <- data$Formula[j]
                Parent <- data$Parent[j]
                Charge <- data$Charge[j]
                Cl_perc <- data$Cl_perc[j]
                dat <- getisotopes(x = as.character(data$Formula[j]))
                dat <- as.data.frame(dat[[1]])
                dat <- dat %>% 
                        mutate(abundance = round(abundance, 1)) %>%
                        mutate(`m/z` = round(`m/z`, 6)) %>%
                        mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`)) %>%
                        mutate(Parent_Formula = Parent) %>%
                        mutate(Cl_perc = Cl_perc) %>%
                        mutate(Adduct_Formula =  Formula) %>%
                        mutate(Charge = Charge) %>%
                        mutate(Isotopologue = case_when(
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
                        mutate(Fragment = paste0(adduct_ions, " ", Isotopologue)) %>%
                        rename(Rel_ab = abundance) %>%
                        select(Parent_Formula, Cl_perc, Charge, Fragment, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
                data_ls[[j]] <- dat
        }
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}
