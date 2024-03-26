###############################################################################
#
#  Automated quantification of chlorinated paraffins (CPquant), VSCCPs, SCCPs, and MCCPs
#  
#  Prepared by Cleo Davie-Martin
#  The Climate and Environmental Research Institute NILU
#  Last updated: 2024-01-11
#
###############################################################################

# Initiate libraries
library(tidyverse); library(reshape2); library(readxl); library(data.table); library(gridExtra); library(viridis)


###############################################################################
### DOCUMENTS TO PREPARE BEFORE STARTING IN R
###############################################################################

# 1. CSV file with Natural isotope abundances (with Congener, mz, NatAbund, Adduct columns). Use "1_Standards_Natural isotope abundances.csv" as an example/template. Note that the mass-to-charge ratio (mz) column must match exactly with those stated in the SCCP and MCCP standard information files!

# 2. Agilent raw data output (.xlsx file): file should have two Compound Method columns (Name=Congener, User Defined=sequential number), followed by four additional columns for every sample included in the batch (Name=Sample_name, Amt.=Amount(Volume), Area(1)=Peak_area_sample, Area(2)=Peak_area_IS). The user should add row 3 (IS_vol_uL, which is the volume in microliters of internal standard added to each sample), row 4 (RS_vol_uL, which is the volume in microliters of recovery standard added to each sample) and row 5 (Blank_name, which is the name of the blank file (or average of blanks) to be applied to the sample). Use "2_SMCCP_Agilent_rawdata_revised.xlsx" as an example/template. 
# NOTE: Blanks and standards in the raw data output are non-sensical and should be ignored (they are applied to the samples separately)

# 3. CSV file with the standard used for the recovery calculations (with Congener, Standard_name, Std_peak_area, and Std_peak_area_IS). The information is copied directly from the Agilent raw data file. Use "3_CPquant_Standard for recovery calculation.csv" as an example/template.

# 4. CSV file with the blanks(s) used for blank-subtraction calculations (with Congener, Blank_name, Amount, Blank_peak_area, Blank_peak_area_IS, Blank_IS_vol_uL, and Blank_RS_vol_uL). Most of the information is copied directly from the Agilent raw data file, the Blank_IS_vol_uL and Blank_RS_vol_uL are inserted manually. Use "4_CPquant_Blank for SCCP_MCCP.csv" as an example/template.
# NOTE: No TCN rows in the blank file normally (perhaps add a fail-safe in the code instead...)

# 5. CSV file with the SCCP standard information used in the reconstructed patterns (with Congener, Mass, N_Cl, mz, and a further column for each of the standards). This information will need to be updated any time a new calibration set is prepared (not so frequently). Use "5_Standards_SCCP_normalized peak areas.csv" as an example/template. Note that the mass-to-charge ratio (mz) column must match exactly with those stated in the natural isotope abundances file!

# 6. CSV file with the MCCP standard information used in the reconstructed patterns (with Congener, Mass, N_Cl, mz, and a further column for each of the standards). This information will need to be updated any time a new calibration set is prepared (not so frequently). Use "6_Standards_MCCP_normalized peak areas.csv" as an example/template. Note that the mass-to-charge ratio (mz) column must match exactly with those stated in the natural isotope abundances file!

# 7. CSV file with the SCCP and MCCP standard calibration curve coefficients (A, B, and C) from the polynomial regression (with Standard, CP, A, B, and C columns. This information will need to be updated any time a new calibration curves are prepared (not so frequently). Use "7_CPquant_Standard calibration curve fit.csv" as an example/template.


###############################################################################
### LINES OF CODE THAT REQUIRE MANUAL MANIPULATION
###############################################################################

# Set the working directory (i.e., where the raw Agilent data file is saved)
setwd("C:/Users/cldm/OneDrive - NILU/Documents/Projects/CPquant/20231002_Anders_problematic")
#setwd("T:/MassHunter/GC-QTOF_quant/CP/R/Working environment/Import templates_C9"); getwd()

# Set the output directory (i.e., where the processed data and figures are saved)
output_directory <- "T:/org/CP/Resultater ny R metode/2024/240115_SMCCP"

# Define the batch name (i.e., raw data batch)
batch_name <- "240115"

# Define the filename for the Natural isotope abundances
filename_NatAbund <- "1_Standards_Natural isotope abundances.csv"

# Define the filename of the raw Agilent data
filename_raw_data <- "2_SMCCP_Agilent_rawdata_revised.xlsx"

# Set the filename of the standard data used for recovery calculations
filename_recovery_std <- "3_CPquant_Standard for recovery calculation.csv"

# Define the filename of the blank data
filename_blank <- "4_CPquant_Blank for SCCP_MCCP.csv"

# Define the filename of the SCCP standards (normalized peak areas)
filename_SCCP_std <- "5_Standards_SCCP_normalized peak areas.csv"

# Define the filename of the MCCP standards (normalized peak areas)
filename_MCCP_std <- "6_Standards_MCCP_normalized peak areas.csv"

# Define the filename of the calibration curve parameters (polynomial regression coefficients)
filename_cal_curve <- "7_CPquant_Standard calibration curve fit.csv"

# Define the internal standard (Sample_IS_conc_pg_uL/Std_IS_conc_pg_uL/Blank_IS_conc_pg_uL) and recovery standard (Sample_RS_conc_pg_uL/Std_RS_conc_pg_uL/Blank_RS_conc_pg_uL) concentrations in the samples, the standards used for recovery calculations, and the blanks 
# Define the recovery standard volume in the standards used for recovery calculations (Std_RS_vol_uL)
Sample_IS_conc_pg_uL = 330.330512997808 # For SCCP, this is 13C10 Hexachlorodecane (35.23)
Sample_RS_conc_pg_uL = 21.3054146 # This is 1,2,3,4 Tetrachloronaphthalene, TCN
Blank_IS_conc_pg_uL = 330.330512997808 # For SCCP, this is 13C10 Hexachlorodecane (35.23)
Blank_RS_conc_pg_uL = 21.3054146 # This is 1,2,3,4 Tetrachloronaphthalene, TCN
Std_IS_conc_pg_uL = 25.41720269 # For SCCP, this is 13C10 Hexachlorodecane (24.19)
Std_RS_conc_pg_uL = 3.416570395 # This is 1,2,3,4 Tetrachloronaphthalene, TCN
Std_IS_vol_uL = 20 # uL
Std_RS_vol_uL = 20 # uL


###############################################################################
### Load in Agilent raw data export file
###############################################################################

df <- read_excel(filename_raw_data, skip=1) # Test file sent by Anders 2023-06-22 (with matching output)
# Note: R console displays list of "New names" for columns with duplicate headers (we will correct these later)


###############################################################################
### Rearrange and reformat Agilent export file
###############################################################################

# Calculate the number of samples in the raw data file
# Note: this assumes that the first two columns are applicable to all samples (Congener, N) and that every sample contains four additional columns worth of data (Sample_name, Amount, Peak_area (sample), Peak_area_IS)
n_sample <- (ncol(df)-2)/4

# Create a revised header
header <- c("Congener", "N", rep(c("Sample_name", "Amount", "Sample_peak_area", "Sample_peak_area_IS"), times=n_sample))

# Apply new header to dataframe
names(df) <- header

# Create a sample list (data rearrangement for ease of use)
sample_seq <- seq(1, n_sample, 1) # Sequence of sample numbers
dl <- list() # Create an empty list
for (i in 1:length(sample_seq)) { # Populate list with sample data
        start=sample_seq[i]*4-1
        stop=sample_seq[i]*4+2
        dl[[i]] <- df[,c(1,2,start:stop)] # All samples have first two columns + 4 sample columns
        dl[[i]]$Sample_IS_vol_uL <- as.numeric(dl[[i]]$Sample_name[1]) # Apply internal standard volume across all congeners of the same sample
        dl[[i]]$Sample_RS_vol_uL <- as.numeric(dl[[i]]$Sample_name[2]) # Apply recovery standard volume across all congeners of the same sample
        dl[[i]]$Blank_type <- dl[[i]]$Sample_name[3] # Apply blank name across all congeners of the same sample
        dl[[i]] <- dl[[i]][-c(1:3),] # Remove first two rows
        dl[[i]]$Amount <- dl[[i]]$Amount[1] # Apply the amount from the first row across all congeners from the same sample (fixes any manual insertion issues)
        dl[[i]]$Sample_name <- dl[[i]]$Sample_name[1] # Apply the first sample name across all congeners from the same sample (fixes any manual renaming issues)
}

# Convert sample list back into a long-format dataframe
df_melt <- data.frame(do.call("rbind", dl))

# Replace NA Sample_peak_area with 0
df_melt <- df_melt %>% 
        mutate(Sample_peak_area_corr=ifelse(is.na(Sample_peak_area), 0, Sample_peak_area))


###############################################################################
### Calculate sample recovery (using the standard; 1,2,3,4 Tetrachloronaphthalene, TCN)
###############################################################################

# Read in the standard information for calculating the recovery
df_recovery_std <- read.csv(filename_recovery_std, header=TRUE)

# Replace NA Sample_peak_area with 0
df_recovery_std <- df_recovery_std %>% 
        mutate(Std_peak_area_corr=ifelse(is.na(Std_peak_area), 0, Std_peak_area))

# Filter out just the recovery standard congener (TCN) and one other congener (e.g., C10Cl4) to get the internal standard information (all congeners have the same Sample_peak_area_IS)
# Calculate the mass of internal standard (Sample_IS_mass_pg) and recovery standard (Sample_RS_mass_pg) added to samples
df_recovery <- df_melt %>% 
        filter(Congener=="TCN") %>% 
        select(Congener, Sample_name, Sample_peak_area_RS=Sample_peak_area_corr, Sample_RS_vol_uL) %>% 
        left_join(df_melt %>% 
                          filter(Congener=="C10Cl4") %>% # Can filter by any congener (the internal standard is the same for all)
                          select(Sample_name, Sample_peak_area_IS, Sample_IS_vol_uL)) %>% 
        mutate(Sample_IS_conc_pg_uL=Sample_IS_conc_pg_uL,
               Sample_RS_conc_pg_uL=Sample_RS_conc_pg_uL,
               Sample_IS_mass_pg=Sample_IS_vol_uL*Sample_IS_conc_pg_uL,
               Sample_RS_mass_pg=Sample_RS_vol_uL*Sample_RS_conc_pg_uL)
# Joining with `by = join_by(Sample_name)`

# Add the standard information for calculating recovery
df_recovery$Std_peak_area_RS <- df_recovery_std$Std_peak_area_corr[nrow(df_recovery_std)]
df_recovery$Std_peak_area_IS <- df_recovery_std$Std_peak_area_IS[1]

# Calculate recoveries
df_recovery <- df_recovery %>% 
        mutate(Std_IS_vol_uL=Std_IS_vol_uL,
               Std_RS_vol_uL=Std_RS_vol_uL,
               Std_IS_conc_pg_uL=Std_IS_conc_pg_uL,
               Std_RS_conc_pg_uL=Std_RS_conc_pg_uL,
               Std_IS_mass_pg=Std_IS_vol_uL*Std_IS_conc_pg_uL,
               Std_RS_mass_pg=Std_RS_vol_uL*Std_RS_conc_pg_uL,
               RF_sample=(Sample_peak_area_RS/Sample_RS_mass_pg)/(Sample_peak_area_IS/Sample_IS_mass_pg),
               RF_std=(Std_peak_area_RS/Std_RS_mass_pg)/(Std_peak_area_IS/Std_IS_mass_pg),
               Recovery_pcnt=RF_std/RF_sample*100)

# Merge recovery back in with original sample data
df_melt <- df_melt %>% 
        left_join(df_recovery %>% 
                          select(Sample_name, Recovery_pcnt))
# Joining with `by = join_by(Sample_name)`


###############################################################################
### Define congeners
###############################################################################

# Add a chlorinated paraffins designation (i.e., VSCCP, SCCP, MCCP, or LCCP) for merging with standard information
# Remove the recovery standard; 1,2,3,4 Tetrachloronaphthalene, TCN (now that we have calculated the recovery)
df_melt <- df_melt %>% 
        filter(Congener!="TCN") %>% 
        mutate(Carbons=as.numeric(gsub("C", "", sapply(strsplit(Congener, split="Cl", fixed=TRUE), function(x) (x[1])))),
               CP=ifelse(is.na(Carbons), NA, ifelse(Carbons<10, "VSCCP", ifelse(Carbons>=10 & Carbons<14, "SCCP", ifelse(Carbons>17, "LCCP", "MCCP")))))

# Separate out the SCCP and MCCP data (although they will be treated in the same manner, they are calculated separately because they have different numbers of standards applied to each and can't merge directly)
df_SCCP <- df_melt %>% 
        filter(CP=="SCCP" | CP=="VSCCP") # We have merged SCCP and VSCCP, because there is only one VSCCP standard...

df_MCCP <- df_melt %>% 
        filter(CP=="MCCP")

df_VSCCP <- df_melt %>% 
        filter(CP=="VSCCP")


###############################################################################
### Exclude samples for which there is no Sample_peak_area_IS
###############################################################################

# Identify which samples do not have associated internal standard peak areas (and therefore cannot be quantified)
IS_fail1 <- df_SCCP %>% 
        group_by(Sample_name) %>% 
        summarise(IS_NA_IS_area=length(which(is.na(Sample_peak_area_IS))), 
                  IS_NA_IS_vol=length(which(is.na(Sample_IS_vol_uL))),
                  IS_NA_blank=length(which(is.na(Blank_type)))) %>% 
        filter(IS_NA_IS_area>0 | IS_NA_IS_vol>0 | IS_NA_blank>0)
IS_fail_samples1 <- unique(IS_fail1$Sample_name) # Names of samples with no Sample_peak_area_IS

IS_fail2 <- df_MCCP %>% 
        group_by(Sample_name) %>% 
        summarise(IS_NA_IS_area=length(which(is.na(Sample_peak_area_IS))), 
                  IS_NA_IS_vol=length(which(is.na(Sample_IS_vol_uL))),
                  IS_NA_blank=length(which(is.na(Blank_type)))) %>% 
        filter(IS_NA_IS_area>0 | IS_NA_IS_vol>0 | IS_NA_blank>0)
IS_fail_samples2 <- unique(IS_fail2$Sample_name)

IS_fail3 <- df_VSCCP %>% 
        group_by(Sample_name) %>% 
        summarise(IS_NA_IS_area=length(which(is.na(Sample_peak_area_IS))), 
                  IS_NA_IS_vol=length(which(is.na(Sample_IS_vol_uL))),
                  IS_NA_blank=length(which(is.na(Blank_type)))) %>% 
        filter(IS_NA_IS_area>0 | IS_NA_IS_vol>0 | IS_NA_blank>0)
IS_fail_samples3 <- unique(IS_fail3$Sample_name)


# Exclude samples from main dataset
df_SCCP <- df_SCCP %>% 
        filter(!Sample_name %in% IS_fail_samples1)

df_MCCP <- df_MCCP %>% 
        filter(!Sample_name %in% IS_fail_samples2)

df_VSCCP <- df_VSCCP %>% 
        filter(!Sample_name %in% IS_fail_samples3)


###############################################################################
### Read in standards template and natural isotope abundances table
###############################################################################

standards_SCCP <- read.csv(filename_SCCP_std, header=TRUE)
standards_MCCP <- read.csv(filename_MCCP_std, header=TRUE)
nat_iso_abund <- read.csv(filename_NatAbund, header=TRUE)

# Make a note of the length of the standards data frame before merge (to ensure there is no difference after the merge)
N_standards_SCCP <- nrow(standards_SCCP)
N_standards_MCCP <- nrow(standards_MCCP)


###############################################################################
### Merge standards and natural isotope abundances with sample data
###############################################################################

# Merge standards (SCCP and MCCP) and natural isotope abundances
standards_SCCP <- standards_SCCP %>% 
        left_join(nat_iso_abund)
# Joining with `by = join_by(Congener, mz)`

standards_MCCP <- standards_MCCP %>% 
        left_join(nat_iso_abund)
# Joining with `by = join_by(Congener, mz)`

# Confirm that the length of the data frame has not changed after the merge (paired numbers should be the same)
N_standards_SCCP_new <- nrow(standards_SCCP); N_standards_SCCP; N_standards_SCCP_new
N_standards_MCCP_new <- nrow(standards_MCCP); N_standards_MCCP; N_standards_MCCP_new

# Index which standards/congeners are missing natural isotope abundances
indx1=which(is.na(standards_SCCP$NatAbund))
standards_SCCP$Congener[indx1] # character(0) means there are no missing data (if there is a mis-match/missing data, it will print which Congener it applies to in the Console - correct before proceeding)
indx2=which(is.na(standards_MCCP$NatAbund))
standards_MCCP$Congener[indx2]

# Merge standards data frames with the sample data
N_df_SCCP <- nrow(df_SCCP) # Data frame length before the merge
df_SCCP_merge <- df_SCCP %>% 
        left_join(standards_SCCP)
# Joining with `by = join_by(Congener)`

N_df_MCCP <- nrow(df_MCCP) # Data frame length before the merge
df_MCCP_merge <- df_MCCP %>% 
        left_join(standards_MCCP)
# Joining with `by = join_by(Congener)`

# Confirm that the length of the data frame has not changed after the merge (paired numbers should be the same)
N_df_SCCP_merge <- nrow(df_SCCP_merge); N_df_SCCP; N_df_SCCP_merge 
N_df_MCCP_merge <- nrow(df_MCCP_merge); N_df_MCCP; N_df_MCCP_merge 


###############################################################################
### Normalize sample peak area to IS and NatAbund
###############################################################################

# Calculate the mass of internal standard (IS) added to samples (Sample_IS_mass_pg)
# Account for natural isotope abundances and IS in the sample peak area
df_SCCP_merge <- df_SCCP_merge %>% 
        mutate(Sample_IS_added_pg=Sample_IS_vol_uL*Sample_IS_conc_pg_uL, # Calculate mass of IS added to the sample
               Sample_IS_mass_pg=Sample_peak_area_IS/Sample_IS_added_pg, # Normalise IS peak area to the mass of IS added
               Sample_peak_area_norm=Sample_peak_area_corr/(NatAbund*Sample_IS_mass_pg)) # Account for IS and NatAbund in the sample peak areas

df_MCCP_merge <- df_MCCP_merge %>% 
        mutate(Sample_IS_added_pg=Sample_IS_vol_uL*Sample_IS_conc_pg_uL, # Calculate mass of IS added to the sample
               Sample_IS_mass_pg=Sample_peak_area_IS/Sample_IS_added_pg, # Normalise IS peak area to the mass of IS added
               Sample_peak_area_norm=Sample_peak_area_corr/(NatAbund*Sample_IS_mass_pg)) # Account for IS and NatAbund in the sample peak areas


###############################################################################
### Load in the blank information
###############################################################################

# Read in blank file (raw data), which may contain one or many blanks
df_blank <- read.csv(filename_blank, header=TRUE)

# Fill the values of Blank_IS_vol_uL and Blank_RS_vol_uL down the column (in case it was only entered once per blank)
# Replace NA Peak_area with 0
df_blank <- df_blank %>% 
        fill(c("Blank_IS_vol_uL", "Blank_RS_vol_uL"), .direction="down") %>% 
        mutate(Blank_peak_area_corr=ifelse(is.na(Blank_peak_area), 0, Blank_peak_area))


###############################################################################
### Calculate blank recovery (using the standard; 1,2,3,4 Tetrachloronaphthalene, TCN)
###############################################################################

# Filter out just the recovery standard congener (TCN) and one other congener (e.g., C10Cl4) to get the internal standard information (all congeners have the same Sample_peak_area_IS)
# Calculate the mass of internal standard (Blank_IS_mass_pg) and recovery standard (Blank_RS_mass_pg) added to samples
df_recovery_blank <- df_blank %>% 
        filter(Congener=="TCN") %>% 
        select(Congener, Blank_name, Blank_peak_area_RS=Blank_peak_area_corr, Blank_IS_vol_uL, Blank_RS_vol_uL) %>% 
        left_join(df_blank %>% 
                          filter(Congener=="C10Cl4") %>% # Can filter by any congener (the internal standard is the same for all)
                          select(Blank_name, Blank_peak_area_IS)) %>% 
        mutate(Blank_IS_conc_pg_uL=Blank_IS_conc_pg_uL,
               Blank_RS_conc_pg_uL=Blank_RS_conc_pg_uL,
               Blank_IS_mass_pg=Blank_IS_vol_uL*Blank_IS_conc_pg_uL,
               Blank_RS_mass_pg=Blank_RS_vol_uL*Blank_RS_conc_pg_uL)
# Joining with `by = join_by(Blank_name)`

# Add the standard information for calculating recovery
df_recovery_blank$Std_peak_area_RS <- df_recovery_std$Std_peak_area_corr[nrow(df_recovery_std)]
df_recovery_blank$Std_peak_area_IS <- df_recovery_std$Std_peak_area_IS[1]

# Calculate recoveries
df_recovery_blank <- df_recovery_blank %>% 
        mutate(Std_IS_vol_uL=Std_IS_vol_uL,
               Std_RS_vol_uL=Std_RS_vol_uL,
               Std_IS_conc_pg_uL=Std_IS_conc_pg_uL,
               Std_RS_conc_pg_uL=Std_RS_conc_pg_uL,
               Std_IS_mass_pg=Std_IS_vol_uL*Std_IS_conc_pg_uL,
               Std_RS_mass_pg=Std_RS_vol_uL*Std_RS_conc_pg_uL,
               RF_blank=(Blank_peak_area_RS/Blank_RS_mass_pg)/(Blank_peak_area_IS/Blank_IS_mass_pg),
               RF_std=(Std_peak_area_RS/Std_RS_mass_pg)/(Std_peak_area_IS/Std_IS_mass_pg),
               Blank_recovery_pcnt=RF_std/RF_blank*100)

# Merge recovery back in with original sample data
df_blank <- df_blank %>% 
        left_join(df_recovery_blank %>% 
                          select(Blank_name, Blank_recovery_pcnt, Blank_IS_added_pg=Blank_IS_mass_pg))
# Joining with `by = join_by(Blank_name)`


###############################################################################
### Define congeners in the blank
###############################################################################

# Add a chlorinated paraffins designation (i.e., VSCCP, SCCP, MCCP, or LCCP) for merging with standard information
# Remove the recovery standard; 1,2,3,4 Tetrachloronaphthalene, TCN (now that we have calculated the recovery)
df_blank <- df_blank %>% 
        filter(Congener!="TCN") %>% 
        mutate(Carbons=as.numeric(gsub("C", "", sapply(strsplit(Congener, split="Cl", fixed=TRUE), function(x) (x[1])))),
               CP=ifelse(is.na(Carbons), NA, ifelse(Carbons<10, "VSCCP", ifelse(Carbons>=10 & Carbons<14, "SCCP", ifelse(Carbons>17, "LCCP", "MCCP")))))

# Separate out the SCCP and MCCP data (although they will be treated in the same manner, they are calculated separately because they have different numbers of standards applied to each and can't merge directly)
df_SCCP_blank <- df_blank %>% 
        filter(CP=="SCCP" | CP=="VSCCP") #### Added "VSCCP" here (it was missing from the blanks and causing the issues when normalized concentrations were calculated)

df_MCCP_blank <- df_blank %>% 
        filter(CP=="MCCP")


###############################################################################
### Merge standards and natural isotope abundances with blank data
###############################################################################

# Merge standards data frames with the sample data
N_df_SCCP_blank <- nrow(df_SCCP_blank) # Data frame length before the merge
df_SCCP_blank <- df_SCCP_blank %>% 
        left_join(standards_SCCP)
# Joining with `by = join_by(Congener)`

N_df_MCCP_blank <- nrow(df_MCCP_blank) # Data frame length before the merge
df_MCCP_blank <- df_MCCP_blank %>% 
        left_join(standards_MCCP)

# Confirm that the length of the data frame has not changed after the merge (paired numbers should be the same)
N_df_SCCP_blank2 <- nrow(df_SCCP_blank); N_df_SCCP_blank; N_df_SCCP_blank2 
N_df_MCCP_blank2 <- nrow(df_MCCP_blank); N_df_MCCP_blank; N_df_MCCP_blank2 


###############################################################################
### Normalize blank peak area to IS and NatAbund
###############################################################################

# Calculate the mass of internal standard (IS) added to samples (Blank_IS_mass_pg)
# Account for natural isotope abundances and IS in the sample peak area
df_SCCP_blank <- df_SCCP_blank %>% 
        mutate(Blank_IS_mass_pg=Blank_peak_area_IS/Blank_IS_added_pg, # Normalise IS peak area to the mass of IS added
               Blank_peak_area_norm=Blank_peak_area_corr/(NatAbund*Blank_IS_mass_pg)) # Account for IS and NatAbund in the blank peak areas

df_MCCP_blank <- df_MCCP_blank %>% 
        mutate(Blank_IS_mass_pg=Blank_peak_area_IS/Blank_IS_added_pg, # Normalise IS peak area to the mass of IS added
               Blank_peak_area_norm=Blank_peak_area_corr/(NatAbund*Blank_IS_mass_pg)) # Account for IS and NatAbund in the blank peak areas


###############################################################################
### Blank subtraction
###############################################################################

# Calculate the average blank
df_SCCP_blank_ave <- df_SCCP_blank %>% 
        group_by(Congener) %>% 
        summarise(Blank_peak_area_norm=mean(Blank_peak_area_norm, na.rm=TRUE)) %>% 
        ungroup

df_MCCP_blank_ave <- df_MCCP_blank %>% 
        group_by(Congener) %>% 
        summarise(Blank_peak_area_norm=mean(Blank_peak_area_norm, na.rm=TRUE)) %>% 
        ungroup

# Add a new Blank_name and add in additional columns
df_SCCP_blank_ave <- df_SCCP_blank_ave %>% 
        mutate(Blank_name="Blank_ave") %>% 
        left_join(df_SCCP_blank  %>% 
                          filter(Blank_name==Blank_name[1]) %>% 
                          select(-Blank_name, -Blank_peak_area_norm))
# Joining with `by = join_by(Congener)`

df_MCCP_blank_ave <- df_MCCP_blank_ave %>% 
        mutate(Blank_name="Blank_ave") %>% 
        left_join(df_MCCP_blank  %>% 
                          filter(Blank_name==Blank_name[1]) %>% 
                          select(-Blank_name, -Blank_peak_area_norm))

# Bind the average results back with the individual blank results
df_SCCP_blank_all <- df_SCCP_blank %>% 
        full_join(df_SCCP_blank_ave)
# Joining with `by = join_by(Congener, Blank_name, Blank_peak_area, Blank_peak_area_IS, Blank_IS_vol_uL, Blank_RS_vol_uL, Blank_peak_area_corr, Blank_recovery_pcnt, Blank_IS_added_pg, Carbons, CP, Mass, N_Cl, mz, X51pcnt, X55pcnt, X63pcnt, C10_50pcnt, C10_65pcnt, C11_50pcnt, C11_65pcnt, C12_50pcnt, C12_65pcnt, C13_50pcnt, C13_65pcnt, NatAbund, Adduct, Blank_IS_mass_pg, Blank_peak_area_norm)`

df_MCCP_blank_all <- df_MCCP_blank %>% 
        full_join(df_MCCP_blank_ave)

# Calculate the total Blank_peak_area_norm for every sample (used later in the solution to the quadratic equation for standards)
df_SCCP_blank_all <- df_SCCP_blank_all %>% 
        left_join(df_SCCP_blank_all %>% 
                          group_by(Blank_name) %>% 
                          summarise(Blank_total=sum(Blank_peak_area_norm, na.rm=TRUE)) %>% 
                          ungroup())
# Joining with `by = join_by(Blank_name)`

df_MCCP_blank_all <- df_MCCP_blank_all %>% 
        left_join(df_MCCP_blank_all %>% 
                          group_by(Blank_name) %>% 
                          summarise(Blank_total=sum(Blank_peak_area_norm, na.rm=TRUE)) %>% 
                          ungroup())

# Merge blanks with samples
# Blank subtraction
# Replace negative blank-subtracted values with 0
df_SCCP_merge <- df_SCCP_merge %>% 
        left_join(df_SCCP_blank_all %>% 
                          select(Congener, Blank_type=Blank_name, Blank_recovery_pcnt, Blank_peak_area_norm, Blank_total)) %>%
        mutate(Sample_peak_area_blank_sub=Sample_peak_area_norm-Blank_peak_area_norm,
               Sample_peak_area_blank_corr=ifelse(Sample_peak_area_blank_sub<0, 0, Sample_peak_area_blank_sub))
# Joining with `by = join_by(Congener, Blank_type)`

df_MCCP_merge <- df_MCCP_merge %>% 
        left_join(df_MCCP_blank_all %>% 
                          select(Congener, Blank_type=Blank_name, Blank_recovery_pcnt, Blank_peak_area_norm, Blank_total)) %>%
        mutate(Sample_peak_area_blank_sub=Sample_peak_area_norm-Blank_peak_area_norm,
               Sample_peak_area_blank_corr=ifelse(Sample_peak_area_blank_sub<0, 0, Sample_peak_area_blank_sub))


###############################################################################
### Normalization of sample peak areas (to sample total peak area)
###############################################################################

# Calculate summed total peak area per sample
# Merge summed total peak areas back with sample data
# Normalize sample peak areas to the total peak area (i.e., fractional contribution)
df_SCCP_merge <- df_SCCP_merge %>% 
        left_join(df_SCCP_merge %>% 
                          group_by(CP, Sample_name) %>% 
                          summarise(Sample_total=sum(Sample_peak_area_blank_corr, na.rm=TRUE)) %>% 
                          ungroup()) %>% 
        mutate(Normalized_conc=Sample_peak_area_blank_corr/Sample_total,
               Normalized_conc=replace(Normalized_conc, is.nan(Normalized_conc), 0))  # Replace NaN with 0 (occurs when blank-subtraction results in all negative concentrations)
# Joining with `by = join_by(Sample_name, CP)`

df_MCCP_merge <- df_MCCP_merge %>% 
        left_join(df_MCCP_merge %>% 
                          group_by(Sample_name) %>% 
                          summarise(Sample_total=sum(Sample_peak_area_blank_corr, na.rm=TRUE)) %>% 
                          ungroup()) %>% 
        mutate(Normalized_conc=Sample_peak_area_blank_corr/Sample_total,
               Normalized_conc=replace(Normalized_conc, is.nan(Normalized_conc), 0))  # Replace NaN with 0 (occurs when blank-subtraction results in all negative concentrations)


###############################################################################
### Separate out final sample and standard information for calculations
###############################################################################

df_SCCP_final <- df_SCCP_merge %>% 
        select(CP, Congener, Carbons, N_Cl, Sample_name, Recovery_pcnt, Amount, Normalized_conc, Sample_total, Blank_type, Blank_recovery_pcnt, Blank_total, C9_48pcnt:C13_60pcnt)

df_MCCP_final <- df_MCCP_merge %>% 
        select(CP, Congener, Carbons, N_Cl, Sample_name, Recovery_pcnt, Amount, Normalized_conc, Sample_total, Blank_type, Blank_recovery_pcnt, Blank_total, C14_49pcnt:C17_56pcnt)


###############################################################################
### The Lawson-Hanson NNLS implemention of non-negative least squares
###############################################################################

# Read in the library packages necessary
library(nnls)
# Note: an algorithm for non-negative linear least squares that solves min||Ax-b||2 with the constraint x>=0

# Set parameters for model fitting
options(digits=3)

# Create a for loop to run through all of the samples in the file and save the nnls model parameters (x estimates)
# First, split the data by Sample_name (each sample will have a separate model)
x1 <- split(df_SCCP_final, f=df_SCCP_final$Sample_name)
x2 <- split(df_MCCP_final, f=df_MCCP_final$Sample_name)

# Create an empty list for the model results to populate
mod_list1=list()
mod_list2=list()
sample_name1=vector()
sample_name2=vector()

# Create a loop that runs the nnls model for each sample and saves the x estimate output to mod_list
# Note: check column numbers for the standards are correct (in creation of M matrix)
for (i in 1:length(x1)) {
        M <- as.matrix(x1[[i]][,13:21]) # Perhaps there is a way to separate out "SCCP standards" (based on name). These column numbers will change if the number of standards changes.
        V <- as.vector(x1[[i]]$Normalized_conc)
        mod <- nnls(M, V)
        mod_list1[[i]] <- mod$x
        sample_name1[i] <- unique(x1[[i]]$Sample_name)
}

for (i in 1:length(x2)) {
        M <- as.matrix(x2[[i]][,13:19]) # Perhaps there is a way to separate out "MCCP standards" (based on name). These column numbers will change if the number of standards changes.
        V <- as.vector(x2[[i]]$Normalized_conc)
        mod <- nnls(M, V)
        mod_list2[[i]] <- mod$x
        sample_name2[i] <- unique(x2[[i]]$Sample_name)
}

# For each model output, calculate the sum and then mod_list_rev is mod_list divided by the sum 
# Note: this step is required because mod_list sums to just below 1 for each sample
mod_list_rev1=list()
mod_list_rev2=list()

for (i in 1:length(mod_list1)) {
        mod_sum <- sum(mod_list1[[i]])
        mod_list_rev1[[i]] <- mod_list1[[i]]/mod_sum
}

for (i in 1:length(mod_list2)) {
        mod_sum <- sum(mod_list2[[i]])
        mod_list_rev2[[i]] <- mod_list2[[i]]/mod_sum
}
# Note: mod_list_rev should now sum to 1 for each sample

# Convert model results list to a dataframe
mod_list_rev1 <- data.frame(do.call("rbind", mod_list_rev1))
mod_list_rev2 <- data.frame(do.call("rbind", mod_list_rev2))

# Rename headers
names(mod_list_rev1) <- paste0(names(df_SCCP_final)[13:21], "_mod")
names(mod_list_rev2) <- paste0(names(df_MCCP_final)[13:19], "_mod")

# Add in sample name
mod_list_rev1$Sample_name <- sample_name1
mod_list_rev2$Sample_name <- sample_name2

# Reverse parameter settings so that it doesn't affect the accuracy 
options(digits=7) # default is 7 digits printed


###############################################################################
### Calculate the reconstructed concentrations
###############################################################################

# Merge model results with df_final
df_SCCP_final <- df_SCCP_final %>% 
        left_join(mod_list_rev1)
# Joining with `by = join_by(Sample_name)`

df_MCCP_final <- df_MCCP_final %>% 
        left_join(mod_list_rev2)
# Joining with `by = join_by(Sample_name)`

# Calculate reconstructed concentrations (normalized standard concentrations x model outputs for the corresponding chlorination level)
# Note: these new columns could change in name and length if the standards change
df_SCCP_final <- df_SCCP_final %>% 
        mutate(C9_48pcnt_recon=C9_48pcnt*C9_48pcnt_mod,
               C10_52pcnt_recon=C10_52pcnt*C10_52pcnt_mod,
               C10_58pcnt_recon=C10_58pcnt*C10_58pcnt_mod,
               C11_52pcnt_recon=C11_52pcnt*C11_52pcnt_mod,
               C11_57pcnt_recon=C11_57pcnt*C11_57pcnt_mod,
               C12_53pcnt_recon=C12_53pcnt*C12_53pcnt_mod,
               C12_57pcnt_recon=C12_57pcnt*C12_57pcnt_mod,
               C13_45pcnt_recon=C13_45pcnt*C13_45pcnt_mod,
               C13_60pcnt_recon=C13_60pcnt*C13_60pcnt_mod)

df_MCCP_final <- df_MCCP_final %>% 
        mutate(C14_49pcnt_recon=C14_49pcnt*C14_49pcnt_mod,
               C14_58pcnt_recon=C14_58pcnt*C14_58pcnt_mod,
               C15_47pcnt_recon=C15_47pcnt*C15_47pcnt_mod,
               C15_59pcnt_recon=C15_59pcnt*C15_59pcnt_mod,
               C16_51pcnt_recon=C16_51pcnt*C16_51pcnt_mod,
               C16_58pcnt_recon=C16_58pcnt*C16_58pcnt_mod,
               C17_56pcnt_recon=C17_56pcnt*C17_56pcnt_mod)


###############################################################################
### Calculate the Normalized sample sum (Normalized_conc_recon) from the reconstructed patterns
###############################################################################

# Normalized_conc_recon is the sum of the reconstructed patterns. It is important that the correct columns are identified to sum (first and last _recon standard)
df_SCCP_final <- df_SCCP_final %>% 
        mutate(Normalized_conc_recon=rowSums(select(., C9_48pcnt_recon:C13_60pcnt_recon), na.rm=TRUE))

df_MCCP_final <- df_MCCP_final %>% 
        mutate(Normalized_conc_recon=rowSums(select(., C14_49pcnt_recon:C17_56pcnt_recon), na.rm=TRUE))


###############################################################################
### Calculate r2 value from linear regression between Normalized_conc and Normalized_conc_recon
###############################################################################

df_SCCP_final <- df_SCCP_final %>% 
        left_join(df_SCCP_final %>%
                          group_by(CP, Sample_name) %>%
                          summarise(Normalized_conc_r2=summary(lm(Normalized_conc ~ Normalized_conc_recon))$r.squared))
# Joining with `by = join_by(CP, Sample_name)`

df_MCCP_final <- df_MCCP_final %>% 
        left_join(df_MCCP_final %>%
                          group_by(Sample_name) %>%
                          summarise(Normalized_conc_r2=summary(lm(Normalized_conc ~ Normalized_conc_recon))$r.squared))


###############################################################################
### Calculate the contribution of each standard
###############################################################################

# Load in standard calibration curve information
df_cal_curve <- read.csv(filename_cal_curve, header=TRUE)

# Separate out summed sample information for each sample (need to calculate each standards contribution for each sample)
df_SCCP_sample <- df_SCCP_final %>% 
        select(CP, Sample_name, Sample_total, Amount, C9_48pcnt_mod:C13_60pcnt_mod) %>% 
        distinct()

df_MCCP_sample <- df_MCCP_final %>% 
        select(Sample_name, Sample_total, Amount, C14_49pcnt_mod:C17_56pcnt_mod) %>% 
        distinct()

# Solve the quadratic equation for x (for each standard and sample combination individually)
# Note: this is a messy and long-winded section of code. I couldn't think how to do it differently. The number in [brackets] represents the row number of the df_cal_curve dataframe (one row per standard). Thus, when the standards change, the row numbers and column names will also change.
for (i in 1:nrow(df_SCCP_sample)) {
        df_SCCP_sample$C9_48pcnt_solve_x[i]=(-df_cal_curve$B[1] + sqrt(df_cal_curve$B[1]^2-(4*df_cal_curve$A[1]*(df_cal_curve$C[1]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[1])
        df_SCCP_sample$C10_52pcnt_solve_x[i]=(-df_cal_curve$B[2] + sqrt(df_cal_curve$B[2]^2-(4*df_cal_curve$A[2]*(df_cal_curve$C[2]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[2])
        df_SCCP_sample$C10_58pcnt_solve_x[i]=(-df_cal_curve$B[3] + sqrt(df_cal_curve$B[3]^2-(4*df_cal_curve$A[3]*(df_cal_curve$C[3]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[3])
        df_SCCP_sample$C11_52pcnt_solve_x[i]=(-df_cal_curve$B[4] + sqrt(df_cal_curve$B[4]^2-(4*df_cal_curve$A[4]*(df_cal_curve$C[4]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[4])
        df_SCCP_sample$C11_57pcnt_solve_x[i]=(-df_cal_curve$B[5] + sqrt(df_cal_curve$B[5]^2-(4*df_cal_curve$A[5]*(df_cal_curve$C[5]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[5])
        df_SCCP_sample$C12_53pcnt_solve_x[i]=(-df_cal_curve$B[6] + sqrt(df_cal_curve$B[6]^2-(4*df_cal_curve$A[6]*(df_cal_curve$C[6]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[6])
        df_SCCP_sample$C12_57pcnt_solve_x[i]=(-df_cal_curve$B[7] + sqrt(df_cal_curve$B[7]^2-(4*df_cal_curve$A[7]*(df_cal_curve$C[7]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[7])
        df_SCCP_sample$C13_45pcnt_solve_x[i]=(-df_cal_curve$B[8] + sqrt(df_cal_curve$B[8]^2-(4*df_cal_curve$A[8]*(df_cal_curve$C[8]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[8])
        df_SCCP_sample$C13_60pcnt_solve_x[i]=(-df_cal_curve$B[9] + sqrt(df_cal_curve$B[9]^2-(4*df_cal_curve$A[9]*(df_cal_curve$C[9]-df_SCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[9])
}
# Note: it is possible to get warning messages if the solutions to the quadratic equation produce negatives (these are corrected for in the next steps)
# Warning messages:
#1: In sqrt(df_cal_curve$B[9]^2 - (4 * df_cal_curve$A[9] * (df_cal_curve$C[9] -  :
#   NaNs produced

for (i in 1:nrow(df_MCCP_sample)) {
        df_MCCP_sample$C14_49pcnt_solve_x[i]=(-df_cal_curve$B[12] + sqrt(df_cal_curve$B[12]^2-(4*df_cal_curve$A[12]*(df_cal_curve$C[12]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[12])
        df_MCCP_sample$C14_58pcnt_solve_x[i]=(-df_cal_curve$B[13] + sqrt(df_cal_curve$B[13]^2-(4*df_cal_curve$A[13]*(df_cal_curve$C[13]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[13])
        df_MCCP_sample$C15_47pcnt_solve_x[i]=(-df_cal_curve$B[14] + sqrt(df_cal_curve$B[14]^2-(4*df_cal_curve$A[14]*(df_cal_curve$C[14]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[14])
        df_MCCP_sample$C15_59pcnt_solve_x[i]=(-df_cal_curve$B[15] + sqrt(df_cal_curve$B[15]^2-(4*df_cal_curve$A[15]*(df_cal_curve$C[15]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[15])
        df_MCCP_sample$C16_51pcnt_solve_x[i]=(-df_cal_curve$B[16] + sqrt(df_cal_curve$B[16]^2-(4*df_cal_curve$A[16]*(df_cal_curve$C[16]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[16])
        df_MCCP_sample$C16_58pcnt_solve_x[i]=(-df_cal_curve$B[17] + sqrt(df_cal_curve$B[17]^2-(4*df_cal_curve$A[17]*(df_cal_curve$C[17]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[17])
        df_MCCP_sample$C17_56pcnt_solve_x[i]=(-df_cal_curve$B[18] + sqrt(df_cal_curve$B[18]^2-(4*df_cal_curve$A[18]*(df_cal_curve$C[18]-df_MCCP_sample$Sample_total[i]))))/(2*df_cal_curve$A[18])
}

# Replace NaN with NA in the solutions to the quadratic equation
df_SCCP_sample[,c(5:22)] <- df_SCCP_sample[,c(5:22)] %>% 
        mutate_all(~replace(., is.nan(.), NA))

df_MCCP_sample[,c(4:17)] <- df_MCCP_sample[,c(4:17)] %>% 
        mutate_all(~replace(., is.nan(.), NA))


# Calculate the contribution of each standard (x * R_mod algorithm)
df_SCCP_sample <- df_SCCP_sample %>% 
        mutate(C9_48pcnt_cont=C9_48pcnt_solve_x*C9_48pcnt_mod,
               C10_52pcnt_cont=C10_52pcnt_solve_x*C10_52pcnt_mod,
               C10_58pcnt_cont=C10_58pcnt_solve_x*C10_58pcnt_mod,
               C11_52pcnt_cont=C11_52pcnt_solve_x*C11_52pcnt_mod,
               C11_57pcnt_cont=C11_57pcnt_solve_x*C11_57pcnt_mod,
               C12_53pcnt_cont=C12_53pcnt_solve_x*C12_53pcnt_mod,
               C12_57pcnt_cont=C12_57pcnt_solve_x*C12_57pcnt_mod,
               C13_45pcnt_cont=C13_45pcnt_solve_x*C13_45pcnt_mod,
               C13_60pcnt_cont=C13_60pcnt_solve_x*C13_60pcnt_mod)

df_MCCP_sample <- df_MCCP_sample %>% 
        mutate(C14_49pcnt_cont=C14_49pcnt_solve_x*C14_49pcnt_mod,
               C14_58pcnt_cont=C14_58pcnt_solve_x*C14_58pcnt_mod,
               C15_47pcnt_cont=C15_47pcnt_solve_x*C15_47pcnt_mod,
               C15_59pcnt_cont=C15_59pcnt_solve_x*C15_59pcnt_mod,
               C16_51pcnt_cont=C16_51pcnt_solve_x*C16_51pcnt_mod,
               C16_58pcnt_cont=C16_58pcnt_solve_x*C16_58pcnt_mod,
               C17_56pcnt_cont=C17_56pcnt_solve_x*C17_56pcnt_mod)

# Calculate the amount in nanograms from the samples (i.e., the sum of the contributions from each standard)
df_SCCP_sample <- df_SCCP_sample %>% 
        mutate(SCCP_mass_ng=rowSums(select(., C10_52pcnt_cont:C13_60pcnt_cont), na.rm=TRUE),
               SCCP_conc_pg_m3=SCCP_mass_ng*1000/Amount,
               SCCP_conc_ng_g=SCCP_conc_pg_m3/1000,
               VSCCP_mass_ng=C9_48pcnt_cont,
               VSCCP_conc_pg_m3=VSCCP_mass_ng*1000/Amount,
               VSCCP_conc_ng_g=VSCCP_conc_pg_m3/1000) # Not sure how we go from pg/m3 to ng/g here (must be some other assumptions required, not shown in the spreadsheet...)

df_MCCP_sample <- df_MCCP_sample %>% 
        mutate(MCCP_mass_ng=rowSums(select(., C14_49pcnt_cont:C17_56pcnt_cont), na.rm=TRUE),
               MCCP_conc_pg_m3=MCCP_mass_ng*1000/Amount,
               MCCP_conc_ng_g=MCCP_conc_pg_m3/1000)

# Merge back with original sample information (at this stage, we only need the final parameters)
df_SCCP_final <- df_SCCP_final %>% 
        left_join(df_SCCP_sample %>% 
                          select(CP, Sample_name, SCCP_mass_ng, SCCP_conc_pg_m3, SCCP_conc_ng_g, VSCCP_mass_ng, VSCCP_conc_pg_m3, VSCCP_conc_ng_g))
# Joining with `by = join_by(Sample_name)`

df_MCCP_final <- df_MCCP_final %>% 
        left_join(df_MCCP_sample %>% 
                          select(Sample_name, MCCP_mass_ng, MCCP_conc_pg_m3, MCCP_conc_ng_g))

# Correct VSCCP and SCCP columns so they are mutually exclusive
indx1=which(df_SCCP_final$CP!="VSCCP")
indx2=which(df_SCCP_final$CP=="VSCCP")
df_SCCP_final[indx1,45:47] <- NA
df_SCCP_final[indx2,42:44] <- NA


###############################################################################
### Calculate the fractional contributions of all of the congeners to the sample
###############################################################################

df_SCCP_final <- df_SCCP_final %>% 
        mutate(Sample_conc_pg_m3=ifelse(CP=="SCCP", Normalized_conc*SCCP_conc_pg_m3, Normalized_conc*VSCCP_conc_pg_m3),
               Sample_conc_ng_g=ifelse(CP=="SCCP", Normalized_conc*SCCP_conc_ng_g, Normalized_conc*VSCCP_conc_ng_g),
               Fractional_cont=ifelse(CP=="SCCP", Normalized_conc_recon*SCCP_conc_ng_g, Normalized_conc_recon*VSCCP_conc_ng_g))

df_MCCP_final <- df_MCCP_final %>% 
        mutate(Sample_conc_pg_m3=Normalized_conc*MCCP_conc_pg_m3,
               Sample_conc_ng_g=Normalized_conc*MCCP_conc_ng_g,
               Fractional_cont=Normalized_conc_recon*MCCP_conc_ng_g)


###############################################################################
### Join SCCP and MCCP results together
###############################################################################

df_final <- rbind(df_SCCP_final %>%
                          filter(CP=="SCCP") %>% 
                          select(Sample_name, Amount, Recovery_pcnt, Blank_type, Blank_recovery_pcnt, CP, Congener, Carbons, N_Cl, Normalized_conc, Normalized_conc_recon, Normalized_conc_r2, Sample_conc_pg_m3, Sample_conc_ng_g, Fractional_cont, VSMCCP_mass_ng=SCCP_mass_ng, VSMCCP_conc_pg_m3=SCCP_conc_pg_m3, VSMCCP_conc_ng_g=SCCP_conc_ng_g),
                  df_SCCP_final %>% 
                          filter(CP=="VSCCP") %>% 
                          select(Sample_name, Amount, Recovery_pcnt, Blank_type, Blank_recovery_pcnt, CP, Congener, Carbons, N_Cl, Normalized_conc, Normalized_conc_recon, Normalized_conc_r2, Sample_conc_pg_m3, Sample_conc_ng_g, Fractional_cont, VSMCCP_mass_ng=VSCCP_mass_ng, VSMCCP_conc_pg_m3=VSCCP_conc_pg_m3, VSMCCP_conc_ng_g=VSCCP_conc_ng_g),
                  df_MCCP_final %>% select(Sample_name, Amount, Recovery_pcnt, Blank_type, Blank_recovery_pcnt, CP, Congener, Carbons, N_Cl, Normalized_conc, Normalized_conc_recon, Normalized_conc_r2, Sample_conc_pg_m3, Sample_conc_ng_g, Fractional_cont, VSMCCP_mass_ng=MCCP_mass_ng, VSMCCP_conc_pg_m3=MCCP_conc_pg_m3, VSMCCP_conc_ng_g=MCCP_conc_ng_g))

# Set working directory where you want the outputs saved
setwd(output_directory); getwd()

# Export a .CSV file of results per sample
export <- split(df_final, f=df_final$Sample_name)

# Export sample file(s)
for(i in 1:length(export)){
        write.csv(export[[i]], paste0(unique(export[[i]]$Sample_name), "_Quantification_VSMCCP.csv"), row.names=FALSE)
}

# Batch-based sample overview
batch_final <- df_final %>%
        group_by(Sample_name, CP, Amount, Recovery_pcnt, Blank_type, Blank_recovery_pcnt, Normalized_conc_r2, VSMCCP_mass_ng, VSMCCP_conc_pg_m3, VSMCCP_conc_ng_g) %>% 
        summarise(N=n()) %>% 
        ungroup() %>% 
        select(-N)

# Export batch_based sample overview
batch_filename <- paste0(batch_name, "_VSMCCP_Batch sample overview.csv")
write.csv(batch_final, batch_filename, row.names=FALSE)


###############################################################################
### Prepare figures of the results
###############################################################################

# Define order of Congeners (by Carbons, then N_Cl)
df_SCCP_final$Congener <- factor(df_SCCP_final$Congener, levels=unique(df_SCCP_final$Congener))
df_MCCP_final$Congener <- factor(df_MCCP_final$Congener, levels=unique(df_MCCP_final$Congener))

# Create VSCCP and SCCP data frames
df_VSCCP <- df_SCCP_final %>% filter(CP=="VSCCP")
df_SCCP <- df_SCCP_final %>% filter(CP=="SCCP")

# Split df_final by Sample_name
y1 <- split(df_SCCP, f=df_SCCP$Sample_name)
y2 <- split(df_MCCP_final, f=df_MCCP_final$Sample_name)
y3 <- split(df_VSCCP, f=df_VSCCP$Sample_name)

# Prepare a figure showing original normalized concentrations (i.e., Normalized_conc)
# Create an empty list to store the figure information in
plot_list1a=list()
plot_list2a=list()
plot_list3a=list()

# For each sample, create a plot
for (i in 1:length(y1)) {
        F <- ggplot(data=y1[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Normalized_conc), fill="darkblue") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ggtitle(paste0(unique(y1[[i]]$Sample_name), "_Original normalized SCCP concentrations")) + annotate(geom="text", x=27, y=max(y1[[i]]$Normalized_conc)*0.95, label=paste0("r2 = ", signif(unique(y1[[i]]$Normalized_conc_r2), digits=3)), col="darkred") + annotate(geom="text", x=27, y=max(y1[[i]]$Normalized_conc)*0.8, label=paste0("total conc (ng/g) = ", ifelse(is.na(unique(y1[[i]]$SCCP_conc_ng_g)), NA, round(unique(y1[[i]]$SCCP_conc_ng_g), 1))), col="darkred")
        plot_list1a[[i]]=F
}

for (i in 1:length(y2)) {
        F <- ggplot(data=y2[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Normalized_conc), fill="darkblue") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ggtitle(paste0(unique(y2[[i]]$Sample_name), "_Original normalized MCCP concentrations")) + annotate(geom="text", x=22, y=max(y2[[i]]$Normalized_conc)*0.95, label=paste0("r2 = ", signif(unique(y2[[i]]$Normalized_conc_r2), digits=3)), col="darkred") + annotate(geom="text", x=22, y=max(y2[[i]]$Normalized_conc)*0.8, label=paste0("total conc (ng/g) = ", ifelse(is.na(unique(y2[[i]]$MCCP_conc_ng_g)), NA, round(unique(y2[[i]]$MCCP_conc_ng_g), 1))), col="darkred")
        plot_list2a[[i]]=F
}

for (i in 1:length(y3)) {
        F <- ggplot(data=y3[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Normalized_conc), fill="darkblue") + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + ggtitle(paste0(unique(y3[[i]]$Sample_name), "_Original normalized VSCCP concentrations")) + annotate(geom="text", x=5, y=(0.1+max(y3[[i]]$Normalized_conc))*0.95, label=paste0("r2 = ", signif(unique(y3[[i]]$Normalized_conc_r2), digits=3)), col="darkred") + annotate(geom="text", x=5, y=(0.1+max(y3[[i]]$Normalized_conc))*0.8, label=paste0("total conc (ng/g) = ", ifelse(is.na(unique(y3[[i]]$VSCCP_conc_ng_g)), NA, round(unique(y3[[i]]$VSCCP_conc_ng_g), 1))), col="darkred")
        plot_list3a[[i]]=F
}

# Prepare a figure showing reconstructed concentrations (i.e., stacked bar plot)
# Select the data we need
df_SCCP_red <- df_SCCP %>% 
        select(Congener, Carbons, N_Cl, Sample_name, C10_52pcnt_recon:C13_60pcnt_recon) %>% 
        arrange(Sample_name, Carbons, N_Cl)

df_MCCP_red <- df_MCCP_final %>% 
        select(Congener, Carbons, N_Cl, Sample_name, C14_49pcnt_recon:C17_56pcnt_recon) %>% 
        arrange(Sample_name, Carbons, N_Cl)

df_VSCCP_red <- df_VSCCP %>% 
        select(Congener, Carbons, N_Cl, Sample_name, C9_48pcnt_recon) %>% 
        arrange(Sample_name, Carbons, N_Cl)

# Rearrange into long format (needed to make a stacked plot)
df_SCCP_melt <- reshape2::melt(df_SCCP_red, id.vars=c("Sample_name", "Congener", "Carbons", "N_Cl"), variable.name="Chlorination", value.name="Reconstructed_conc")
df_MCCP_melt <- reshape2::melt(df_MCCP_red, id.vars=c("Sample_name", "Congener", "Carbons", "N_Cl"), variable.name="Chlorination", value.name="Reconstructed_conc")
df_VSCCP_melt <- reshape2::melt(df_VSCCP_red, id.vars=c("Sample_name", "Congener", "Carbons", "N_Cl"), variable.name="Chlorination", value.name="Reconstructed_conc")

# Number of standards and colour palette
n_SCCP_std <- length(unique(df_SCCP_melt$Chlorination))
n_MCCP_std <- length(unique(df_MCCP_melt$Chlorination))
n_VSCCP_std <- length(unique(df_VSCCP_melt$Chlorination))

# Split df_final_melt by Sample_name
z1 <- split(df_SCCP_melt, f=df_SCCP_melt$Sample_name)
z2 <- split(df_MCCP_melt, f=df_MCCP_melt$Sample_name)
z3 <- split(df_VSCCP_melt, f=df_VSCCP_melt$Sample_name)

# Create an empty list to store the figure information in
plot_list1b=list()
plot_list2b=list()
plot_list3b=list()

# For each sample, create a plot
for (i in 1:length(z1)) {
        F <- ggplot(data=z1[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Reconstructed_conc, fill=Chlorination), position="stack") + scale_fill_viridis(discrete=TRUE)  + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="bottom") + ggtitle(paste0(unique(z1[[i]]$Sample_name), "_Reconstructed normalized SCCP concentrations"))
        plot_list1b[[i]]=F
}

for (i in 1:length(z2)) {
        F <- ggplot(data=z2[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Reconstructed_conc, fill=Chlorination), position="stack") + scale_fill_viridis(discrete=TRUE)  + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="bottom") + ggtitle(paste0(unique(z2[[i]]$Sample_name, "_Reconstructed normalized MCCP concentrations")))
        plot_list2b[[i]]=F
}

for (i in 1:length(z3)) {
        F <- ggplot(data=z3[[i]], na.rm=TRUE) + geom_col(aes(x=Congener, y=Reconstructed_conc, fill=Chlorination), position="stack") + scale_fill_viridis(discrete=TRUE)  + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position="bottom") + ggtitle(paste0(unique(z3[[i]]$Sample_name), "_Reconstructed normalized VSCCP concentrations"))
        plot_list3b[[i]]=F
}

# Set working directory where you want the outputs saved
setwd(output_directory); getwd()

# Prepare a .PDF output for each Sample_name with both original and reconstructed normalized concentrations displayed in a single .PDF
for(i in 1:length(plot_list1a)){
        ggsave(file=paste0(unique(z1[[i]]$Sample_name), ".pdf"), plot=grid.arrange(plot_list3a[[i]], plot_list1a[[i]], plot_list2a[[i]], plot_list3b[[i]], plot_list1b[[i]], plot_list2b[[i]], ncol=3), width=20, height=8, unit="in")
}

