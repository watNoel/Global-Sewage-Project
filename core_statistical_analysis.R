# Author: Noel Waters

# Overview:
# - Load and clean data from synthetic and natural wastewater community assays.
# - Fit Poisson and Negative Binomial models; compare fits and select NB models if they fit better.
# - Upon selecting models, Apply Benjamini-Hochberg multiple comparison correction.
# - Calculate confidence intervals for resistance rates and rate ratios using emmeans.
# - Three regression setups:
#   1. Synthetic community with 0h as reference.
#   2. Synthetic community with 72h as reference.
#   3. Natural wastewater community with 0h as reference.

{
  library(tidyverse)    # Data wrangling
  library(MASS)         # Negative Binomial modeling
  library(broom)        # Tidying model outputs
  library(magrittr)     # Pipe operators (%<>%, %>%)
  library(multcomp)     # Multiple comparison adjustments
  library(emmeans)      # Estimated marginal means
}

# Load synthetic community data

indat_all <- readxl::read_xlsx(path = "./data/CFU_data.xlsx",sheet = "SEC", guess_max = 2000)

# For the script to run, the column names of the files have to be renamed to these: 
indat_all <- indat_all |> rename (Number=Number,
                                  Name=Name,
                                  Run_Experiment=Run_Experiment,
                                  Target_Antibiotic=Antibiotic,
                                  Antibiotic_Plate_1=Antibiotic_Plates,
                                  Antibiotic_Plate_2="...6",
                                  Antibiotic_Plate_3="...7",
                                  MH_Plate_1=No_Antibiotic_Plates,
                                  MH_Plate_2="...9",
                                  MH_Plate_3="...10"
)


{
  # skip first 4 rows since they have no data
  indat <- indat_all[4:nrow(indat_all), ]
  indat$Name <- as.factor(indat$Name)
  indat$Batch <- as.factor(substr(indat$Run_Experiment, 1, 1))

  
  # Assign water names and numbers
  tmp_name <- "_NAME_"
  water_num <- 0
  for (i in 1:nrow(indat)) {
    if (!is.na(indat$Number[i])) {
      tmp_name <- indat$Name[i]
      water_num <- water_num + 1
    } else {
      indat$Name[i] <- tmp_name
      indat$Number[i] <- water_num
    }
  }
  
  # Convert plate counts to numeric
  indat %<>% mutate_at(vars(starts_with("Antibiotic_Plate"), starts_with("MH_Plate")), as.numeric)
  
  # Calculate sum and median of replicates
  indat$sum_count <- apply(indat[, c("Antibiotic_Plate_1", "Antibiotic_Plate_2", "Antibiotic_Plate_3")], 1, sum)
  indat$sum_mh <- apply(indat[, c("MH_Plate_1", "MH_Plate_2", "MH_Plate_3")], 1, sum)
  indat$sum_rate <- indat$sum_count / indat$sum_mh
  
  indat$med_count <- apply(indat[, c("Antibiotic_Plate_1", "Antibiotic_Plate_2", "Antibiotic_Plate_3")], 1, median)
  indat$med_mh <- apply(indat[, c("MH_Plate_1", "MH_Plate_2", "MH_Plate_3")], 1, median)
  indat$med_rate <- indat$med_count / indat$med_mh
  
  # Clear temporary variables
  rm(list = c("tmp_name", "i", "water_num"))
}

# Set reference level for modeling
indat <- within(indat, Name <- relevel(Name, ref = "Saline_0h_control"))
indat_sec <- indat
rm(indat, indat_all)

# Define antibiotics used
antibiotics <- unique(indat_sec$Target_Antibiotic)
#("AMC", "CIP", "CTX", "TOB", "SXT")





# Load natural wastewater community data
indat_all_nwc <- readxl::read_xlsx(path = "./CFU_data.xlsx",sheet = "NWC", guess_max = 2000) |> 
  rename (Number=Number,
          Name=Name,
          Run_Experiment=Run_Experiment,
          Target_Antibiotic=Antibiotic,
          ECC_Atb_Plate_1=Antibiotic_Plates,
          ECC_Atb_Plate_2="...6",
          ECC_Atb_Plate_3="...7",
          ECC_no_Atb_Plate_1=No_Antibiotic_Plates,
          ECC_no_Atb_Plate_2="...9",
          ECC_no_Atb_Plate_3="...10"
  )
  

{
  # skip first 4 rows since they have no data
  indat_nwc <- indat_all_nwc[4:nrow(indat_all_nwc), ]
  indat_nwc$Name <- as.factor(indat_nwc$Name)
  indat_nwc$Batch <- as.factor(substr(indat_nwc$Run_Experiment, 1, 1))
  
  # Assign water names and numbers
  tmp_name <- "_NAME_"
  water_num <- 0
  for (i in 1:nrow(indat_nwc)) {
    if (!is.na(indat_nwc$Number[i])) {
      tmp_name <- indat_nwc$Name[i]
      water_num <- water_num + 1
    } else {
      indat_nwc$Name[i] <- tmp_name
      indat_nwc$Number[i] <- water_num
    }
  }
  
  # Convert plate counts to numeric
  indat_nwc %<>% mutate_at(vars(starts_with("ECC_Atb_Plate"), starts_with("ECC_no_Atb_Plate")), as.numeric)
  
  # Calculate sum and median of replicates
  indat_nwc$sum_count <- apply(indat_nwc[, c("ECC_Atb_Plate_1", "ECC_Atb_Plate_2", "ECC_Atb_Plate_3")], 1, sum)
  indat_nwc$sum_mh <- apply(indat_nwc[, c("ECC_no_Atb_Plate_1", "ECC_no_Atb_Plate_2", "ECC_no_Atb_Plate_3")], 1, sum)
  indat_nwc$sum_rate <- indat_nwc$sum_count / indat_nwc$sum_mh
  
  indat_nwc$med_count <- apply(indat_nwc[, c("ECC_Atb_Plate_1", "ECC_Atb_Plate_2", "ECC_Atb_Plate_3")], 1, median)
  indat_nwc$med_mh <- apply(indat_nwc[, c("ECC_no_Atb_Plate_1", "ECC_no_Atb_Plate_2", "ECC_no_Atb_Plate_3")], 1, median)
  indat_nwc$med_rate <- indat_nwc$med_count / indat_nwc$med_mh
  
  # Clean temporary variables
  rm(list = c("tmp_name", "i", "water_num"))
}

rm(indat_all_nwc)

# Set reference level for modeling
indat_nwc <- within(indat_nwc, Name <- relevel(Name, ref = "Saline_0h_control"))
antibiotics_nwc <- unique(indat_nwc$Target_Antibiotic)






# Calculate mean, SD, and intervals- NOT based on any models but for an overview:
{
  # Synthetic community
  manual_calc_sec <- indat_sec |>
    group_by(Target_Antibiotic, Name) |>
    summarise(mean_rate = mean(med_rate), nobs = n(), sd_rate = sd(med_rate)) |>
    mutate(
      lower_ci_rate = mean_rate - 1.96 * sd_rate / sqrt(nobs),
      upper_ci_rate = mean_rate + 1.96 * sd_rate / sqrt(nobs),
      rate_ratio = mean_rate / mean_rate[Name == "Saline_0h_control"],
      rate_diff = mean_rate - mean_rate[Name == "Saline_0h_control"]
    ) |>
    ungroup() |>
    group_by(Name) |>
    mutate(
      avg_rate_diff = mean(rate_diff),
      avg_rate_ratio = mean(rate_ratio)
    ) |>
    ungroup() |>
    arrange(avg_rate_diff, Name) |>
    mutate(name_ord = factor(Name, levels = unique(Name)))
  
  # Natural wastewater community
  manual_calc_nwc <- indat_nwc |>
    group_by(Target_Antibiotic, Name) |>
    summarise(mean_rate = mean(med_rate), nobs = n(), sd_rate = sd(med_rate)) |>
    mutate(
      lower_ci_rate = mean_rate - 1.96 * sd_rate / sqrt(nobs),
      upper_ci_rate = mean_rate + 1.96 * sd_rate / sqrt(nobs),
      rate_ratio = mean_rate / mean_rate[Name == "Saline_0h_control"],
      rate_diff = mean_rate - mean_rate[Name == "Saline_0h_control"]
    ) |>
    ungroup() |>
    group_by(Name) |>
    mutate(
      avg_rate_diff = mean(rate_diff),
      avg_rate_ratio = mean(rate_ratio)
    ) |>
    ungroup() |>
    arrange(avg_rate_diff, Name) |>
    mutate(name_ord = factor(Name, levels = unique(Name)))
}




# Next step:
# Fit Poisson and Negative Binomial models to assess overdispersion and model fit.
# Note: Although functionality for using 'sum' as a response variable is included,
# it was never actually used in modeling. 'median' remains the default.

counts_modeller <- function(indat, antibiotic, response = "median") {
  # Purpose:
  # Fit Poisson and Negative Binomial models for a given antibiotic.
  # Uses water source (Name) as a covariate.
  # Supports two response types: 'median' (default) and 'sum'.
  
  # Filter and select relevant columns for modeling
  for_model <- indat |>
    filter(Target_Antibiotic == antibiotic) |>
    dplyr::select(Name, Batch, med_count, med_mh, sum_count, sum_mh)
  
  # Define model formula based on selected response type
  if (response == "median") {
    model.formula <- as.formula(med_count ~ Name + offset(log(med_mh)))
  } else if (response == "sum") {
    model.formula <- as.formula(sum_count ~ Name + offset(log(sum_mh)))
  } else {
    print("ERROR: 'response' must be either 'median' or 'sum'")
  }
  
  # Fit Poisson model
  pois.fit <- glm(
    formula = model.formula,
    data = for_model,
    family = "poisson"
  )
  
  # Extract Poisson model coefficients and fit statistics
  pois.dat <- broom::tidy(pois.fit, conf.int = TRUE, conf.level = 0.95) |>
    mutate(Family = "poisson", Antibiotic = antibiotic, response = response)
  
  pois.glance <- broom::glance(pois.fit) |>
    mutate(Family = "poisson", Antibiotic = antibiotic)
  
  # Fit Negative Binomial model
  
  nb.fit <- glm.nb(
    formula = model.formula,
    data = for_model,
    control=glm.control(maxit = 50) # increasing iterations to aid convergence
  )
  
  # Extract NB model coefficients and fit statistics
  nb.dat <- broom::tidy(nb.fit, conf.int = TRUE, conf.level = 0.95) |>
    mutate(Family = "negbin", Antibiotic = antibiotic, response = response)
  
  nb.glance <- broom::glance(nb.fit) |>
    mutate(Family = "negbin", Antibiotic = antibiotic)
  
  # Package results into lists
  pois_list <- list(
    model = pois.fit,
    coeffs = pois.dat,
    fit = pois.glance
  )
  
  nb_list <- list(
    model = nb.fit,
    coeffs = nb.dat,
    fit = nb.glance
  )
  
  # Return both model results in a named list
  return(antibiotic = list(
    pois_res = pois_list,
    nb_res = nb_list
  ))
}

# Example usage (commented out):
# Fit models for each antibiotic using median response
# median_models_0 <- purrr::map(as.list(antibiotics), .f = \(x) {
#   counts_modeller(indat = indat_sec, response = "median", antibiotic = x)
# })

# Combine model coefficients across antibiotics
# median_coeffs_nb_and_poi_0 <- purrr::map_df(median_models_0, .f = \(x) {
#   rbind(x[["pois_res"]][["coeffs"]], x[["nb_res"]][["coeffs"]])
# })

# Combine model fit statistics across antibiotics
# median_fits_0 <- purrr::map_df(median_models_0, .f = \(x) {
#   rbind(x[["pois_res"]][["fit"]], x[["nb_res"]][["fit"]])
# })




corrected_pval_and_ci_2sided <- function(x, calc_fwer = FALSE) {
  # Purpose:
  # Takes output from counts_modeller and calculates p-values and confidence intervals
  # for estimated coefficients using multiple comparison corrections.
  # Supports both FDR (Benjamini-Hochberg) and FWER (single-step) adjustments.
  
  # Perform Dunnett-style multiple comparisons (two-sided)
  rht_eq <- glht(x[["nb_res"]][["model"]], linfct = mcp(Name = "Dunnett"), alternative = "two.sided")
  
  # FDR adjustment (Benjamini-Hochberg)
  summary_fdr_eq <- summary(rht_eq, test = adjusted("BH"))
  
  if (calc_fwer == TRUE) {
    # FWER adjustment using single-step method (not used in paper)
    set.seed(1234)
    summary_ss_eq <- summary(rht_eq, test = adjusted("single-step", maxpts = 2500 * 100))$test$pvalues
    
    set.seed(1234)
    ci_eq <- confint(rht_eq, level = 0.95, calpha = adjusted_calpha()) |>
      broom::tidy() |>
      dplyr::select(conf.low, conf.high)
    
    # Combine results into a data frame (NA for intercept)
    mc_df <- data.frame(
      term = names(summary_fdr_eq$coef),
      fdr_pval_eq = c(NA, summary_fdr_eq$test$pvalues),
      singlestep_pval_eq = c(NA, summary_ss_eq),
      simul_conf_low_eq = c(NA, ci_eq$conf.low),
      simul_conf_high_eq = c(NA, ci_eq$conf.high)
    )
  } else {
    # Only FDR adjustment (NA for intercept)
    mc_df <- data.frame(
      term = names(summary_fdr_eq$coef),
      fdr_pval_eq = c(NA, summary_fdr_eq$test$pvalues)
    )
  }
  
  # Merge with model coefficients and format output
  out <- merge(x[["nb_res"]][["coeffs"]], mc_df, by = "term") |>
    mutate(Name = gsub("Name", "", term))
  
  # Rename intercept to reference level
  out$Name[out$Name == "(Intercept)"] <- x[["nb_res"]][["model"]]$xlevels$Name[1]
  
  return(out)
}


correcting_with_emmeans <- function(x, .adjust = "none") {
  # Purpose:
  # Uses emmeans to calculate marginal means and contrasts for the negative binomial model.
  # Returns resistance rates and rate ratios with confidence intervals and p-values.
  # Adjustment methods include "none", "mvt" (FWER), or any supported by emmeans.
  
  set.seed(1234)  # Required for reproducibility when using stochastic methods like mvt
  
  # Estimate marginal means (rates)
  marginal_mean_rates <- emmeans(x[["nb_res"]][["model"]], ~Name, type = "response", offset = 0)
  
  # Confidence intervals for marginal means
  marginal_mean_rate_cis <- confint(marginal_mean_rates, level = 0.95, side = "both", adjust = .adjust)
  
  # Dunnett contrasts for rate ratios
  dunnett_contrasts <- contrast(marginal_mean_rates, method = "dunnett")
  
  set.seed(1234)
  dunnett_contrasts_summary <- summary(
    dunnett_contrasts,
    infer = c(TRUE, TRUE),
    side = "both",
    adjust = .adjust,
    level = 0.95
  )
  
  # Format marginal means
  marginal_mean_rate_cis_for_join <- marginal_mean_rate_cis |>
    broom::tidy() |>
    rename(
      rate = response,
      SE_rate = std.error,
      lower_ci_rate = conf.low,
      upper_ci_rate = conf.high
    ) |>
    dplyr::select(-df)
  
  # Format contrasts
  for_left_join <- dunnett_contrasts_summary |>
    mutate(Name = gsub(" /.*", "", contrast)) |>
    rename(
      rate_ratio = ratio,
      SE_rate_ratio = SE,
      lower_ci_rate_ratio = asymp.LCL,
      upper_ci_rate_ratio = asymp.UCL,
      p_value_rate_ratio = p.value
    ) |>
    dplyr::select(Name, everything()) |>
    dplyr::select(-contrast, -df) |>
    mutate(ci_adjustment = .adjust)
  
  # Combine marginal means and contrasts
  result <- marginal_mean_rate_cis_for_join |>
    left_join(for_left_join, by = "Name") |>
    mutate(
      Antibiotic = x[["nb_res"]][["coeffs"]]$Antibiotic[1],
      response = x[["nb_res"]][["coeffs"]]$response[1]
    )
  
  return(result)
}



## Synthetic Community – 0h as Reference

# Fit models for each antibiotic using median response
median_models_sec <- purrr::map(as.list(antibiotics), .f = \(x) {
  counts_modeller(indat = indat_sec, response = "median", antibiotic = x)
})

# Optional: Model diagnostics
# par(mfrow = c(2, 2))
# lapply(1:5, \(i) plot(median_models_sec[[i]][["nb_res"]]$model))

# Extract coefficients and emmeans results
median_coeffs_sec <- purrr::map_df(median_models_sec, .f = corrected_pval_and_ci_2sided)
median_emmeans_sec <- purrr::map_df(median_models_sec, .f = correcting_with_emmeans)

# Merge and clean
merge_df_sec <- merge(median_coeffs_sec, median_emmeans_sec, by = c("Antibiotic", "Name", "response")) |>
  dplyr::select(-term, -response, -Family)

# Select relevant columns for paper
sec_0_control <- merge_df_sec |>
  dplyr::select(
    Antibiotic, Name, estimate, std.error, conf.low, conf.high, p.value, fdr_pval_eq,
    rate, SE_rate, lower_ci_rate, upper_ci_rate,
    rate_ratio, SE_rate_ratio, lower_ci_rate_ratio, upper_ci_rate_ratio,
    ci_adjustment
  )


## Synthetic Community – 72h as Reference

# Relevel reference to 72h and remove 0h
indat_sec_72 <- within(indat_sec, Name <- relevel(Name, ref = "Saline_72h_control")) |>
  filter(Name != "Saline_0h_control")

# Fit models
median_models_sec_72 <- purrr::map(as.list(antibiotics), .f = \(x) {
  counts_modeller(indat = indat_sec_72, response = "median", antibiotic = x)
})

# Extract results
median_coeffs_sec_72 <- purrr::map_df(median_models_sec_72, .f = corrected_pval_and_ci_2sided)
median_emmeans_sec_72 <- purrr::map_df(median_models_sec_72, .f = correcting_with_emmeans)

# Merge and clean
merge_df_sec_72 <- merge(median_coeffs_sec_72, median_emmeans_sec_72, by = c("Antibiotic", "Name", "response")) |>
  dplyr::select(-term, -response, -Family)

# Select relevant columns
sec_72_control <- merge_df_sec_72 |>
  dplyr::select(
    Antibiotic, Name, estimate, std.error, conf.low, conf.high, p.value, fdr_pval_eq,
    rate, SE_rate, lower_ci_rate, upper_ci_rate,
    rate_ratio, SE_rate_ratio, lower_ci_rate_ratio, upper_ci_rate_ratio,
    ci_adjustment
  )


### Natural Wastewater Community – 0h as Reference
# Fit models
median_models_nwc <- purrr::map(as.list(antibiotics_nwc), .f = \(x) {
  counts_modeller(indat = indat_nwc, response = "median", antibiotic = x)
})

# Extract results
median_coeffs_nwc <- purrr::map_df(median_models_nwc, .f = corrected_pval_and_ci_2sided)
median_emmeans_nwc <- purrr::map_df(median_models_nwc, .f = correcting_with_emmeans)

# Merge and clean
merge_df_nwc <- merge(median_coeffs_nwc, median_emmeans_nwc, by = c("Antibiotic", "Name", "response")) |>
  dplyr::select(-term, -response, -Family)

# Select relevant columns
nwc_0_control <- merge_df_nwc |>
  dplyr::select(
    Antibiotic, Name, estimate, std.error, conf.low, conf.high, p.value, fdr_pval_eq,
    rate, SE_rate, lower_ci_rate, upper_ci_rate,
    rate_ratio, SE_rate_ratio, lower_ci_rate_ratio, upper_ci_rate_ratio,
    ci_adjustment
  )



# Export cleaned results to Excel files

  # Modify the output path as applicable:
  
  writexl::write_xlsx(sec_0_control, "./output/sec_fys0h_as_ref.xlsx")
  writexl::write_xlsx(nwc_0_control, "./output/nwc_fys0h_as_ref.xlsx")
  writexl::write_xlsx(sec_72_control, "./output/sec_fys72h_as_ref.xlsx")


