
#---START----Test params for Skyline----START---------- 


Cmin <- 10
Cmax  <- 11
Clmin <- 7
Clmax <- 8
Adduct <- "[CP+Cl]-"
threshold <- 50
MSresolution <- 60000

C <- as.integer(Cmin:Cmax)
Cl <- as.integer(Clmin:Clmax)
selectedAdducts <- Adduct
Adducts <- as.character(selectedAdducts)

CP_allions <- list()
for (i in seq_along(Adducts)) {
        input <- getSkyline(adduct_ions = Adducts[i], C = C, Cl = Cl, threshold = threshold)
        CP_allions <- rbind(CP_allions, input)
}

CP_allions_glob <- CP_allions

CP_allions_skyline <- as_tibble(CP_allions_glob) %>%
        mutate(`Molecule List Name` = paste0("C", `12C`)) %>%
        rename(`Molecule Name` = Parent_Formula) %>%
        mutate(`Molecular Formula` = case_when(
                `37Cl` == 0 ~ paste0("C", `12C`, "H", `1H`, "Cl", `35Cl`),
                `37Cl` > 0 ~ paste0("C", `12C`, "H", `1H`, "Cl", `35Cl`, "Cl'", `37Cl`)
                )) %>%
        mutate(`Precursor Adduct` = str_replace(Fragment, "\\].*", "]")) %>%
        mutate(`Precursor Adduct` = str_replace(`Precursor Adduct`, "\\[.*\\+", "[M+")) %>%
        rename(`Precursor Charge` = Charge) %>%
        add_column(`Explicit Retention Time` = NA) %>%
        add_column(`Explicit Retention Time Window` = NA) %>%
        mutate(Note = case_when(
                `12C` < 10 ~ "vSCCP",
                `12C` >= 10 & `12C` < 14 ~ "SCCPs",
                `12C` >= 14 & `12C` < 18 ~ "MCCPs",
                `12C` >= 18 ~ "LCCPs")) %>%
        select(`Molecule List Name`, `Molecule Name`, `Molecular Formula`, `Precursor Adduct`, 
               `Precursor Charge`, `Explicit Retention Time`, `Explicit Retention Time Window`)
   
       
        


#----END---Test params for Skyline-----END--------- 

test <- dat |> 
        mutate(`Molecule List Name` = case_when(str_detect(Fragment, "(?<=.)CP(?=.)") == TRUE ~ paste0("CP-C", `12C`),
                                      str_detect(Fragment, "(?<=.)CO(?=.)") == TRUE ~ paste0("CO-C", `12C`)))

dat2 <- 




