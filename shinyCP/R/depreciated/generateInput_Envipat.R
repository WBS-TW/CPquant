

generateInput_Envipat <- function(data = data, group = group) {
        
       
        
        #######################
        
        if (fragment_ions == "-Cl") { # Generate fragments M-Cl for each homolog formula  
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "-H") { # Generate fragments M-H for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(H = H-1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "-HCl") { # Generate fragments M-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Cl = Cl-1) %>%
                        mutate(H = H-1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "+Cl") { # Generate fragments M+Cl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl+1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl, Br)
                
        } else if (fragment_ions == "-Cl-HCl") { # Generate fragments M-Cl-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-2) %>%
                        mutate(H = H-1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "-Cl-2HCl") { # Generate fragments M-Cl2HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-3) %>%
                        mutate(H = H-2) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "-Cl-3HCl") { # Generate fragments M-Cl-3HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-4) %>%
                        mutate(H = H-3) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "-Cl-4HCl") { # Generate fragments M-Cl-3HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-5) %>%
                        mutate(H = H-4) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        }else if (fragment_ions == "-2Cl-HCl") { # Generate fragments M-2Cl-HCl for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>% 
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl-3) %>%
                        mutate(H = H-1) %>%
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl)
                
        } else if (fragment_ions == "+Br") { # Generate fragments M+Br for each homolog formula
                data <- data %>%
                        mutate(Parent = Formula) %>%
                        mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                                     group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0),
                                                     group == "BCA" ~ round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0))) %>%
                        mutate(Adduct = adduct_ions) %>%
                        mutate(Cl = Cl) %>%
                        mutate(Br = Br + 1) |> # need to add another Br for [BCA+Br] adduct 
                        mutate(Adduct_formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) %>%
                        select(Parent, Halo_perc, Charge, Adduct, Adduct_formula, C, H, Cl, Br)
        }
}
