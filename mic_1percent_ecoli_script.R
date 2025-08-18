# Author: Noel Waters
# Purpose: Calculate 1% MICs for E. coli based on EUCAST data.
# Algorithm inspired by :Bengtsson-Palme et al 2016 (Concentrations of antibiotics predicted to select for resistant bacteria: Proposed limits for environmental regulation)


# Essentially, it calculates the 1% mic of the wilt-type strains, but also considers that there needs to be at least 10 observations with a lower-or equal that MIC than the MIC selected to represent the 1% mic.

library(tidyverse)


# 1. Load example MIC data from EUCAST - modify the data path as needed but note the required format of the input file from the example.
# Read MIC histogram data

data_path<- "data/example_eucast_ecoli_mic_data.xlsx"
eucast_mics <- readxl::read_xlsx(path = data_path, sheet = "Sheet1", range = "A1:T86")




# Transpose and clean
t_eucast_mics <- data.frame(t(eucast_mics[, -1]))
colnames(t_eucast_mics) <- eucast_mics$Antibiotic_Conc
eucast_mics <- t_eucast_mics |>
  tibble::rownames_to_column(var = "Antibiotic_Conc") |>
  mutate(Antibiotic_Conc = as.numeric(Antibiotic_Conc))

rm(t_eucast_mics)

# Read ECOFF values
  eucast_eccofs <- readxl::read_xlsx(path = data_path, sheet = "Sheet1", range = "A1:X86")[, c(1, 23)]

# Ensure concentration column is numeric
  eucast_mics$Antibiotic_Conc <- as.numeric(eucast_mics$Antibiotic_Conc)

# 2. Filter MICs above ECOFF (keep wild-type only)
{
  eucast_mic_wt <- eucast_mics
  
  for (colname in colnames(eucast_mic_wt)[-1]) {
    print(colname)
    eucast_eccof <- eucast_eccofs$`(T)ECOFF`[eucast_eccofs$Antibiotic_Conc == colname]
    print(eucast_eccof)
    
    # Use ECOFF if available, otherwise set to 512
    if (!eucast_eccof %in% c("ID", "-")) {
      as_numeric_euec <- as.numeric(eucast_eccof)
    } else {
      as_numeric_euec <- 512
    }
    print(as_numeric_euec)
    
    # Set MIC counts above ECOFF to 0
    eucast_mic_wt[eucast_mics$Antibiotic_Conc > as_numeric_euec, colname] <- 0
  }
}
rm(as_numeric_euec, eucast_eccof, colname)


# 3. Determine MIC level with at least 10 observations lower than or equal to this level 
{
  quantlst <- list()
  tmpidx <- 0
  
  for (colname in colnames(eucast_mic_wt)[-1]) {
    tmpidx <- tmpidx + 1
    tmpsum <- 0
    min10 <- NA
    
    for (row in 1:nrow(eucast_mic_wt)) {
      tmpsum <- tmpsum + eucast_mic_wt[row, colname]
      if (tmpsum >= 10) {
        min10 <- eucast_mic_wt$Antibiotic_Conc[row]
        break
      }
    }
    
    quantlst[[tmpidx]] <- data.frame(Antibiotic = colname, min10 = min10)
  }
  
  tenthobsdf <- do.call(rbind, quantlst)
}
rm(quantlst, row, tmpidx, tmpsum, colname, min10)

# 4. Calculate MIC quantiles (1%, 10%, 50%), only 1% is needed for the output.
{
mic100lst <- list()
  
  for (colname in colnames(eucast_mic_wt)[-1]) {
    count_version <- eucast_mic_wt |>
      dplyr::select(Antibiotic_Conc, !!rlang::sym(colname)) |>
      uncount(!!rlang::sym(colname))
    
    mic100lst[[colname]] <- data.frame(
      Antibiotic = colname,
      mic100 = quantile(count_version$Antibiotic_Conc, probs = 0.01),
      mic10 = quantile(count_version$Antibiotic_Conc, probs = 0.10),
      mic50 = quantile(count_version$Antibiotic_Conc, probs = 0.50)
    )
  }
  
  
  mic100df <- do.call(rbind, mic100lst)
  rm(mic100lst,count_version)
  }

# 5. Merge quantile and observation threshold data, and select the 10th lowest value as mic1percent if this value is higher than the 1% mic ( mic100).

merged_min10_and_mic <- merge(tenthobsdf, mic100df, by = "Antibiotic") |>
  mutate(mic1perc = ifelse(min10 > mic100, min10, mic100))



# Final output, a data frame of antibiotic and estimated 1% MIC
mic_ecoli<-merged_min10_and_mic |> dplyr::select(Antibiotic,mic1perc)

writexl::write_xlsx(mic_ecoli, "output/mic_1percent_ecoli.xlsx")


