

getAdduct <- function(adduct_ions, C, Cl, Clmax, threshold) {
        
        # Regex to extract strings
        ion_modes <- str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
        fragment_ions <- str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
        group <- str_extract(adduct_ions, "[^\\[].{2}") # Using positive lookbehind for [)
        
        if (group == "PCA") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl < Clmax) %>% # limit chlorine atoms. 
                        mutate(H = 2*C+2-Cl) %>% # add H atoms
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Formula, C, H, Cl) # move Formula to first column
        } else if (group == "PCO") {
                data <- crossing(C, Cl) %>% 
                        filter(C >= Cl) %>% 
                        filter(Cl < Clmax) %>% 
                        mutate(H = 2*C-Cl) %>% 
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% 
                        select(Formula, C, H, Cl) 
        }  else {
                print("Input not correct, only PCA or PCO is allowed")
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
        
        if (fragment_ions == "-Cl") { # Generate fragments M-Cl for each homolog formula  
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-H") { # Generate fragments M-H for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(H = H-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-HCl") { # Generate fragments M-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(H = H-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "+Cl") { # Generate fragments M+Cl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl+1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-Cl-HCl") { # Generate fragments M-Cl-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-2) %>%
                        mutate(H = H-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-Cl-2HCl") { # Generate fragments M-Cl2HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-3) %>%
                        mutate(H = H-2) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-Cl-3HCl") { # Generate fragments M-Cl-3HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-4) %>%
                        mutate(H = H-3) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        } else if (fragment_ions == "-Cl-4HCl") { # Generate fragments M-Cl-3HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-5) %>%
                        mutate(H = H-4) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)
        }else if (fragment_ions == "-2Cl-HCl") { # Generate fragments M-2Cl-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-3) %>%
                        mutate(H = H-1) %>%
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl)

        } else if (fragment_ions == "+Br") { # Generate fragments M+Br for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                   group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl) %>%
                        mutate(Br = 1) |> 
                        mutate(Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Formula, C, H, Cl, Br)
        }
        
        # Remove formula without Cl after adduct formations
        data <- data %>%
                filter(Cl > 0)
        
        # Create empty list for all ion formulas
        CP_allions <- list()
        data_ls <- list()
        
        
        # function to get isotopic patterns for all PCAs. Limit the threshold to 10%, neutral form. data("isotopes") needs to be loaded first
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, 
                                                        chemforms = x, 
                                                        threshold = threshold, 
                                                        emass = 0.00054857990924,
                                                        plotit = FALSE, 
                                                        charge = Charge)}
        
        if (fragment_ions == "+Br") { 
                for (j in seq_along(data$Formula)) {
                        Formula <- data$Formula[j]
                        Parent <- data$Parent[j]
                        Charge <- data$Charge[j]
                        Halo_perc <- data$Halo_perc[j]
                        dat <- getisotopes(x = as.character(data$Formula[j]))
                        dat <- as.data.frame(dat[[1]])
                        dat <- dat %>% 
                                mutate(abundance = round(abundance, 1)) %>%
                                mutate(`m/z` = round(`m/z`, 6)) %>%
                                mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`, "[79Br]", `79Br`, "[81Br]", `81Br`)) %>%
                                mutate(Parent_Formula = Parent) %>%
                                mutate(Halo_perc = Halo_perc) %>%
                                mutate(Adduct_Formula =  Formula) %>%
                                mutate(Charge = Charge) %>%
                                mutate(Isotopologue = case_when(
                                        `13C` + (`37Cl`+`81Br`)*2 == 0 ~ "",
                                        `13C` + (`37Cl`+`81Br`)*2 == 1 ~ "+1",
                                        `13C` + (`37Cl`+`81Br`)*2 == 2 ~ "+2",
                                        `13C` + (`37Cl`+`81Br`)*2 == 3 ~ "+3",
                                        `13C` + (`37Cl`+`81Br`)*2 == 4 ~ "+4",
                                        `13C` + (`37Cl`+`81Br`)*2 == 5 ~ "+5",
                                        `13C` + (`37Cl`+`81Br`)*2 == 6 ~ "+6",
                                        `13C` + (`37Cl`+`81Br`)*2 == 7 ~ "+7",
                                        `13C` + (`37Cl`+`81Br`)*2 == 8 ~ "+8",
                                        `13C` + (`37Cl`+`81Br`)*2 == 9 ~ "+9",
                                        `13C` + (`37Cl`+`81Br`)*2 == 10 ~ "+10",
                                        `13C` + (`37Cl`+`81Br`)*2 == 11 ~ "+11",
                                        `13C` + (`37Cl`+`81Br`)*2 == 12 ~ "+12",
                                        `13C` + (`37Cl`+`81Br`)*2 == 13 ~ "+13",
                                        `13C` + (`37Cl`+`81Br`)*2 == 14 ~ "+14",
                                        `13C` + (`37Cl`+`81Br`)*2 == 15 ~ "+15",
                                        `13C` + (`37Cl`+`81Br`)*2 == 16 ~ "+16",
                                        `13C` + (`37Cl`+`81Br`)*2 == 17 ~ "+17",
                                        `13C` + (`37Cl`+`81Br`)*2 == 18 ~ "+18",
                                        `13C` + (`37Cl`+`81Br`)*2 == 19 ~ "+19",
                                        `13C` + (`37Cl`+`81Br`)*2 == 20 ~ "+20")) %>%
                                mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) %>%
                                rename(Rel_ab = abundance) %>%
                                select(Parent_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
                        data_ls[[j]] <- dat
                }
        }else {
                for (j in seq_along(data$Formula)) {
                        Formula <- data$Formula[j]
                        Parent <- data$Parent[j]
                        Charge <- data$Charge[j]
                        Halo_perc <- data$Halo_perc[j]
                        dat <- getisotopes(x = as.character(data$Formula[j]))
                        dat <- as.data.frame(dat[[1]])
                        dat <- dat %>% 
                                mutate(abundance = round(abundance, 1)) %>%
                                mutate(`m/z` = round(`m/z`, 6)) %>%
                                mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`)) %>%
                                mutate(Parent_Formula = Parent) %>%
                                mutate(Halo_perc = Halo_perc) %>%
                                mutate(Adduct_Formula =  Formula) %>%
                                mutate(Charge = Charge) %>%
                                mutate(Isotopologue = case_when(
                                        `13C` + (`37Cl`)*2 == 0 ~ "",
                                        `13C` + (`37Cl`)*2 == 1 ~ "+1",
                                        `13C` + (`37Cl`)*2 == 2 ~ "+2",
                                        `13C` + (`37Cl`)*2 == 3 ~ "+3",
                                        `13C` + (`37Cl`)*2 == 4 ~ "+4",
                                        `13C` + (`37Cl`)*2 == 5 ~ "+5",
                                        `13C` + (`37Cl`)*2 == 6 ~ "+6",
                                        `13C` + (`37Cl`)*2 == 7 ~ "+7",
                                        `13C` + (`37Cl`)*2 == 8 ~ "+8",
                                        `13C` + (`37Cl`)*2 == 9 ~ "+9",
                                        `13C` + (`37Cl`)*2 == 10 ~ "+10",
                                        `13C` + (`37Cl`)*2 == 11 ~ "+11",
                                        `13C` + (`37Cl`)*2 == 12 ~ "+12",
                                        `13C` + (`37Cl`)*2 == 13 ~ "+13",
                                        `13C` + (`37Cl`)*2 == 14 ~ "+14",
                                        `13C` + (`37Cl`)*2 == 15 ~ "+15",
                                        `13C` + (`37Cl`)*2 == 16 ~ "+16",
                                        `13C` + (`37Cl`)*2 == 17 ~ "+17",
                                        `13C` + (`37Cl`)*2 == 18 ~ "+18",
                                        `13C` + (`37Cl`)*2 == 19 ~ "+19",
                                        `13C` + (`37Cl`)*2 == 20 ~ "+20")) %>%
                                mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) %>%
                                rename(Rel_ab = abundance) %>%
                                select(Parent_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
                        data_ls[[j]] <- dat
                }
        }
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}
