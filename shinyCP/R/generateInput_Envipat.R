# Two functions in this file
# generateInput_Envipat
# generateInput_Envipat_BCA

generateInput_Envipat <- function(data = data, group = group) {
        
        
        data <- data |> 
                #mutate(Parent = Formula) |>  
                mutate(Halo_perc = case_when(group == "PCA" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C+2-Cl) + 35.45*Cl)*100, 0),
                                             group == "PCO" ~ round(35.45*Cl / (12.01*C + 1.008*(2*C-Cl) + 35.45*Cl)*100, 0))) %>%
                mutate(Adduct = adduct_ions) |> 
                mutate(Cl = case_when(
                        fragment_ions == "-Cl" ~ Cl-1,
                        fragment_ions == "-HCl" ~ Cl-1,
                        fragment_ions == "+Cl" ~ Cl+1,
                        fragment_ions == "-Cl-HCl" ~ Cl-2,
                        fragment_ions == "-Cl-2HCl" ~ Cl-3,
                        fragment_ions == "-Cl-3HCl" ~ Cl-4,
                        fragment_ions == "-Cl-4HCl" ~ Cl-5,
                        fragment_ions == "-2Cl-HCl" ~ Cl-3,
                        .default = Cl)) |> 
                mutate(H = case_when(
                        fragment_ions == "-H" ~ H-1,
                        fragment_ions == "-HCl" ~ H-1,
                        fragment_ions == "-Cl-HCl" ~ H-1,
                        fragment_ions == "-Cl-2HCl" ~ H-2,
                        fragment_ions == "-Cl-3HCl" ~ H-3,
                        fragment_ions == "-Cl-4HCl" ~ H-4,
                        fragment_ions == "-2Cl-HCl" ~ H-1,
                        .default = H)) |> 
                mutate(Br = case_when(
                        fragment_ions == "+Br" ~ 1,
                        TRUE ~0
                )) |> 
                mutate(Adduct_Formula = case_when(
                        fragment_ions != "+Br" ~ paste0("C", C, "H", H, "Cl", Cl),
                        fragment_ions == "+Br" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br))) |> 
                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl)
        
        return(data)
}





# For bromochloroalkanes
generateInput_Envipat_BCA <- function(data = data, group = group) {
        
        
        data <- data |> 
                #mutate(Parent = Formula) |>  
                mutate(Halo_perc = round((35.45*Cl+79.90*Br) / (12.01*C + 1.008*(2*C-Cl-Br) + 35.45*Cl+79.90*Br)*100, 0)) |> 
                mutate(Adduct = adduct_ions) |> 
                mutate(Cl = case_when(
                        fragment_ions == "-Cl" ~ Cl-1,
                        fragment_ions == "-HCl" ~ Cl-1,
                        fragment_ions == "+Cl" ~ Cl+1,
                        fragment_ions == "-Cl-HCl" ~ Cl-2,
                        fragment_ions == "-Cl-2HCl" ~ Cl-3,
                        fragment_ions == "-Cl-3HCl" ~ Cl-4,
                        fragment_ions == "-Cl-4HCl" ~ Cl-5,
                        fragment_ions == "-2Cl-HCl" ~ Cl-3,
                        .default = Cl)) |> 
                mutate(H = case_when(
                        fragment_ions == "-H" ~ H-1,
                        fragment_ions == "-HCl" ~ H-1,
                        fragment_ions == "-Cl-HCl" ~ H-1,
                        fragment_ions == "-Cl-2HCl" ~ H-2,
                        fragment_ions == "-Cl-3HCl" ~ H-3,
                        fragment_ions == "-Cl-4HCl" ~ H-4,
                        fragment_ions == "-2Cl-HCl" ~ H-1,
                        .default = H)) |> 
                mutate(Adduct_Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) |> 
                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, C, H, Cl, Br)
        
        return(data)
}
