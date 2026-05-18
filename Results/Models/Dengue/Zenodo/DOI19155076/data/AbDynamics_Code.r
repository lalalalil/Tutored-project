
# PROJECT:  Longitudinal antibody profiling after dengue reveals distinct dynamics by antibody specificity over 18 months
# FILE:     Supporting code for manuscript
# AUTHORS (Roles):   Jose Victor Zambrana (Data analysis, Visualization), Sandra Bos (Visualization) 
# CREATED:  07/01/2024
# UPDATED:  23/03/2024



######################################################################3
# Libraries ----------------------------------------------------------
######################################################################3

library(tidyverse) # tidyverse 2.0.0
library(lme4) # version 1.1-37, linear mixed models
library(lspline) # version 1.0-0, construct linear splines
library(ggeffects) # version 2.3.1, for calculate marginal effects
library(broom.mixed) # broom.mixed version 0.2.9.6, for model output summaries 
library(lmeresampler) # version 0.2.4, for linear mixed model case bootstrapping – Note: Make sure the package is well installed.
library(emmeans) # version 2.3.1, for calculate marginal effects
library(tictoc) # version 1.2.1, for measuring running time
library(foreach) # version 1.5.2, for parallel processing
library(doParallel) # version 1.0.17, for parallel processing
library(ggsignif) # for plotting statistical significance in ggplot
library(brms) # version 0.6.4, for bayesian modeling
library(cmdstanr) # version 0.8.0, for bayesian modeling – if not working use : install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
library(bayesplot) # version 1.14.0, for plotting bayesian modeles
library(rstatix) # version 0.7.2, for calculating statistical tests by stratification
library(tidyverse)


######################################################################3
# Read dataset ----------------------------------------------------
######################################################################3

D_AP4.00125 = readRDS("AbDynamics_SimulatedDataset.RData")


######################################################################3
# Data management ------------------------------------------------- 
######################################################################3

D_AP4.00125 = D_AP4.00125 %>%
  mutate(
    # Identify if the virus matches the infecting serotype (DV)
    Virus.resp = case_when(DV == 1 & Virus == "DENV1"~ "HOMOTYPIC",
                           DV == 3 & Virus == "DENV3"~ "HOMOTYPIC",
                           Virus == "DENV2"~ "DV2_XR",
                           Virus == "DENV4"~ "DV4_XR",
                           Virus == "ZIKV"~ "ZV_XR"),
    # Simplify Secondary Antibody labels (remove "-PE" suffix)
    `Secondary Ab` = case_when(`Secondary Ab` == "Total IgG-PE" ~ "IgG",
                               `Secondary Ab` == "IgG1-PE" ~ "IgG1",
                               `Secondary Ab` == "IgG2-PE" ~ "IgG2",
                               `Secondary Ab` == "IgG3-PE" ~ "IgG3",
                               `Secondary Ab` == "IgG4-PE" ~ "IgG4",
                               `Secondary Ab` == "IgA-PE" ~ "IgA",
                               `Secondary Ab` == "IgM-PE" ~ "IgM", T ~ `Secondary Ab`),
    Timepoint = factor(Timepoint, levels = c("Cv","3","6","18")),
    # Simplify Antigen levels  
    Antigen.bead = case_when(Antigen.bead == "recE" ~ "E", T ~ Antigen.bead)) 


aim1.dataset2 = D_AP4.00125 %>%
  # Subset the data to focus specifically on Homotypic vs. Cross-Reactive dynamics
  filter((Virus != "ZIKV" & Virus.resp == "HOMOTYPIC")| # Homotypic
           (Virus == "DENV2"& Virus.resp == "DV2_XR") |  # Heterotypic DENV2
           (Virus == "DENV4"& Virus.resp == "DV4_XR") |  # Heterotypic DENV4
           (Virus == "ZIKV" & Virus.resp == "ZV_XR")) %>% # Heterotypic ZIKV
  # Log2 transformation of MFI to normalize variance
  mutate(MFI = log2(MFI),
         cutoff = log2(cutoff),
         # categorization of variables
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR", "DV4_XR","ZV_XR")),
         Antigen.bead = as.factor(Antigen.bead),
         `Secondary Ab` = as.factor(`Secondary Ab`),
         Code = factor(Code),
         IR = as.factor(IR),
         DV = as.factor(DV))



######################################################################3
# Bootstrapping by question  ------------------------------------------ 
######################################################################3

#  modeling of antibody decay using a stratified, parallelized bootstrapping approach. It is designed to evaluate three core biological variables: Antigen Type, Viral Serotype, and Infection History (IR).
# The pipeline consists of three identical structural blocks that slice the data to answer specific questions
# Each block follows this internal logic:
#   
#   Stratification: Subsets the data by Isotype, Infection History, and Virus Response to isolate specific immune profiles.
# 
#   Spline Modeling: Fits a Linear Mixed-Effects Model (lme4) using a linear spline with knots at 100 and 200 days. This captures non-linear decay phases.
# 
#   Case Bootstrapping: Performs Level-2 resampling (resampling participants/Code with replacement) while keeping all longitudinal observations for a participant intact. This ensures the variance estimates account for individual-level heterogeneity.


# CI Overlap Test 
#  This function determines if a reference group's estimate falls within the 
# confidence intervals of all other groups in the dataset.


ci_overlap_test <- function(data, filter_var, reference_label) {
  
  # Convert the filter_var into a symbol for tidy evaluation
  filter_var <- enquo(filter_var)
  
  # Identify the reference row based on the specified variable and label
  ref_row <- data %>% filter(!!filter_var == reference_label)
  
  
  # Extract reference bound
  ref <- ref_row$norm.estimate
  
  
  # Add a column indicating whether each row overlaps with the reference
  data %>%
    mutate(
      overlap_with_ref = norm.lower <= ref & norm.upper >= ref
    )
}




# Define the variables
secondary_ab <- c("IgA", "IgG1", "IgG2", "IgG3", "IgG4", "IgM", "IgG")
antigen_bead <- c("EDIII", "NS1", "E")
IR <- c("P", "S")
virus_resp <- c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR")



######################################################################3
## Antigen question  ------------------------------------------ 
######################################################################3

### PART 1: Parallel bootstrapping pipeline  ##########


# Estimate MFI kinetics and slopes  using nested case-bootstrapping
# by Antigen =  MFI ~ Antigen.bead*lspline(DPSO, knots = c(100,200)) + Age + Sex  + (1|Code)

tic()
n.simulations = 3

# 1. Parallel Configuration

# Set up parallel backend to use multiple processors
cl <- makeCluster(detectCores() - 1) # Leave one core free
registerDoParallel(cl)

# Loop through each combination in parallel
# 2. Nested Parallel Loop 
# Iterates through strata: Secondary Ab -> IR Status -> Virus Response

ant.results_list_strat_paral<- foreach(i = secondary_ab, .combine=list, .multicombine=TRUE, .packages=c('lme4', 'boot', 'ggeffects', 'dplyr', "lmeresampler", "lspline",  "emmeans"), .errorhandling = 'pass') %:%
  foreach(k = IR, .combine='c', .errorhandling = 'pass') %:%
  foreach(l = virus_resp, .combine='c', .errorhandling = 'pass') %dopar% {
    
    
    # A. Data Subsetting
    subset_data <- subset(aim1.dataset2, `Secondary Ab` == i  & IR == k & Virus.resp == l) 
    
    # B. Reference Model Fitting
    model.ref = lmer(MFI ~ Antigen.bead*lspline(DPSO, knots = c(100,200)) + Age + Sex  + (1|Code), 
                     data = subset_data) 
    
    # C. Case Bootstrapping - Fixed Effects
    # Resamples 'Code' (Level 2) but not through timepoints per participant
    bootstraps = case_bootstrap(model.ref, fixef, n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    # D. Case Bootstrapping - Slopes (Trends)
    # Extracts the rate of change (trends) at specific DPSO milestones using emtrends
    bootstraps.at18 = case_bootstrap(model.ref, 
                                     
                                     function(m) {
                                       out <- emtrends(
                                         m,               # pass the refitted model here
                                         specs = ~ Antigen.bead | DPSO, 
                                         var = "DPSO", 
                                         at = list(DPSO = c(20,100,200)),
                                         # (other emtrends arguments go here)
                                         lmerTest.limit = 19530, 
                                         pbkrtest.limit = 19530,
                                         data = subset_data
                                       )
                                       as_tibble(summary(out, infer = T, confint = T))
                                     }
                                     
                                     , n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    
    # E. Case Bootstrapping - Marginal Predictions
    bootstraps.predictions = case_bootstrap(
      model.ref, 
      function(model) {
        tibble(ggpredict(model, terms = c("DPSO [20, 100, 200, 550]", "Antigen.bead"), ci.lvl = NA))
      }, 
      n.simulations, .refit = TRUE, resample = c(TRUE, FALSE), orig_data = subset_data
    )
    
    # F. Summarization
    bootstraps.confint = bootstraps %>% confint()
    
    # G. Result Packaging
    combo_name <- paste(i, k, l, sep = "_")
    named_result <- list(combo_name = list(boot.lmer = bootstraps, boot.confint = bootstraps.confint, predictions = bootstraps.predictions, ratesat18 = bootstraps.at18))
    names(named_result) <- combo_name
    named_result
    
  }


stopCluster(cl) # Stop the cluster
toc()

### --- PART 2: Result Aggregation &  Formatting --- ######################

# 1. Flatten the nested list structure by one level
# do not print this, it is very expensive
ant.results_list_strat_paral = unlist(ant.results_list_strat_paral, recursive = FALSE)



# 2. Initialize empty data frames (containers) for each result type
ant.combined_predictions_strat          <- data.frame() # Fixed effect replicates
ant.combined_predictions_strat_observed <- data.frame() # Observed point estimates
ant.combined_predictions_strat_summary  <- data.frame() # Confidence Intervals (summary)
ant.combined_predictions_strat_timepoints <- data.frame() # Marginal predictions (for plots)
ant.combined_predictions_strat_rates    <- data.frame() # Slopes/Trends from emtrends

# 3. Iterate through each unique Antibody_IR_Virus combination

# Loop through the results_list to extract and combine predicted.lmer data frames
for (name in names(ant.results_list_strat_paral)) {
  
  # # --- Extraction ---
  predicted_df <- ant.results_list_strat_paral[[name]]$boot.lmer$replicates 
  summary_df <- ant.results_list_strat_paral[[name]]$boot.lmer %>% confint()
  predicted.mg <- ant.results_list_strat_paral[[name]]$predictions$replicates 
  observed_df <- ant.results_list_strat_paral[[name]]$boot.lmer$observed 
  extracted_rates <- ant.results_list_strat_paral[[name]]$ratesat18$replicates
  
  # --- Transformation & Labeling ---
  predicted_df <- predicted_df %>%
    mutate(Combination = name)
  
  summary_df <- summary_df %>%
    mutate(Combination = name)
  
  predicted.mg <- predicted.mg %>%
    mutate(Combination = name)
  
  #  handling for Observed Estimates
  observed_df <- observed_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "coef") %>%
    rename(obs = 2) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    mutate(Combination = name) 
  
  
  extracted_rates <- extracted_rates %>%
    mutate(Combination = name)
  
  
  # --- Re-integration ---
  # Stack (rbind) the current results into the master data frames
  
  ant.combined_predictions_strat <- bind_rows(ant.combined_predictions_strat, predicted_df)
  
  ant.combined_predictions_strat_observed <- bind_rows(ant.combined_predictions_strat_observed, observed_df)
  
  ant.combined_predictions_strat_summary <- bind_rows(ant.combined_predictions_strat_summary, summary_df)
  
  ant.combined_predictions_strat_timepoints <- bind_rows(ant.combined_predictions_strat_timepoints, predicted.mg)
  
  ant.combined_predictions_strat_rates <- bind_rows(ant.combined_predictions_strat_rates, extracted_rates)
}


### --- PART 3: Statistical Inference & Slope Categorization ###########################

ant.combined_predictions_strat_summary2 = ant.combined_predictions_strat_summary %>%
  # 1. Map spline terms to actual Timepoints (DPSO)
  mutate(DPSO = case_when(grepl("))1", term) ~ 20, 
                          grepl("))2", term) ~ 100,
                          grepl("))3", term) ~ 200),
         # 2. Identify which Antigen is being compared
         Antigen = case_when(grepl("EDIII", term) ~ "EDIII",
                             grepl("NS1", term) ~ "NS1"),
         # 3. Flag Interaction Terms
         interaction = grepl(":", term)) %>%
  
  # 4. Filter for specific results
  filter(interaction == T,
         type == "norm") %>% # Use the 'normal' bootstrap calculation method
  
  # 5. Deconstruct the 'Combination' string into Tidy variables
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         # Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                Virus.resp == "ZV" ~ "ZV_XR",
                                T ~ Virus.resp)) %>%
  # 6. Enforce Factor levels for proper plot ordering
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR", "DV4_XR","ZV_XR"))) %>%
  select(-term) %>%
  # 7. SIGNIFICANCE TESTING: Categorize the relative slope
  # If the CI (lower to upper) does not cross 0, the difference is significant.
  mutate(rel.slope.sig = case_when(  (lower < 0 & upper < 0) ~ "Sig. Waning",
                                     (lower > 0 & upper > 0) ~ "Sig. Rising",
                                     T ~ "Not Sig.")) %>%
  # Rename columns to indicate these are RELATIVE estimates (compared to E)
  rename(rel.estimate = estimate, rel.lower = lower, rel.upper = upper)


### --- PART 4: Trend Summarization & Half-Life Calculation #############################

ant.combined_predictions_strat_rates.summary = ant.combined_predictions_strat_rates %>%
  # 1. Clean up labels and extract metadata from the 'Combination' string
  mutate(Antigen = Antigen.bead) %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         # Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                T ~ Virus.resp)) %>%
  # 2. Aggregate Bootstrap Replicates
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(DPSO.trend),
            sd = sd(DPSO.trend),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd, # Lower bound of 95% CI
            norm.upper = mean + qnorm(0.975) * sd # Upper bound of 95% CI
  ) %>%
  select(-mean, -sd, - n) %>%
  # 3. Enforce Factor Ordering
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR"))) %>%
  group_by(Combination, DPSO ,Isotype, IR, Virus.resp) %>%
  # 4. Statistical Comparisons (Point-in-Interval Test)
  # Check if the 'E' antigen's trend falls within the CI of other antigens
  group_modify(~ ci_overlap_test(.x, filter_var = Antigen, reference_label = "E")) %>%
  mutate(overlap_with_ref = case_when(Antigen == "E" ~ "Ref.", T ~ as.character(overlap_with_ref ))) %>%
  # 5. Significance Categorization (Absolute Slope)
  mutate(slope.sig = case_when(  (norm.lower < 0 & norm.upper < 0) ~ "Sig. Waning",
                                 (norm.lower > 0 & norm.upper > 0) ~ "Sig. Rising",
                                 T ~ "Not Sig.")) %>%
  # 6. Join with Relative Slopes & Final Labeling
  # Pull in the 'rel.slope.sig' calculated in the previous step
  left_join(ant.combined_predictions_strat_summary2) %>%
  mutate(rel.slope.sig = case_when(is.na(rel.slope.sig) ~ "Ref.",
                                   T ~ rel.slope.sig),
         # Map DPSO values to chronological study phases
         timepoint = case_when(DPSO == "20" ~ "Cv-3M",
                               DPSO == "100" ~ "3M-6M",
                               DPSO == "200" ~ "6M-18M"),
         # 7. Antibody Half-Life (Years)
         hl.estimate = round((-1/norm.estimate)/365.25,2),
         hl.norm.lower = round((-1/norm.lower)/365.25,2),
         hl.norm.upper = round((-1/norm.upper)/365.25,2)) 


### --- PART 5: Marginal Prediction (MFI curves) ########################################

# previous part summarized the slopes (how fast levels change), this script summarizes the actual levels (how much antibody is present) at days 20, 100, 200, and 550.


ant.combined_predictions_strat_timepoints2 = ant.combined_predictions_strat_timepoints %>%
  # 1. Standardize column names from ggeffects output
  # 'group' corresponds to Antigen, 'x' corresponds to DPSO
  rename(Antigen = group, DPSO = x) %>%
  # 2. Deconstruct Metadata
  # Split the 'Combination' string into individual analysis variables
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         # Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                Virus.resp == "ZV" ~ "ZV_XR",
                                T ~ Virus.resp)) %>%
  # 3. Aggregate Bootstrap Predictions
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(predicted ),
            sd = sd(predicted ),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd,
            norm.upper = mean + qnorm(0.975) * sd
  ) %>%
  # 4. Clean up and Factorize
  select(-mean, -sd, - n) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR"))) 






######################################################################3
## Viral response question  ------------------------------------------ 
######################################################################3

### PART 1: Parallel bootstrapping pipeline  ##########

tic()
n.simulations = 3 # 1000 simulations 

# Set up parallel backend to use multiple processors
cl <- makeCluster(detectCores() - 1) # Leave one core free
registerDoParallel(cl)

# Loop through each combination in parallel
vir.results_list_strat_paral <- foreach(i = secondary_ab, .combine=list, .multicombine=TRUE, .packages=c('lme4', 'boot', 'ggeffects', 'dplyr', "lmeresampler", "lspline", "emmeans"), .errorhandling = 'pass') %:%
  foreach(j = antigen_bead, .combine='c', .errorhandling = 'pass') %:%
  foreach(k = IR, .combine='c', .errorhandling = 'pass')  %dopar% {
    
    subset_data <- subset(aim1.dataset2, `Secondary Ab` == i  & IR == k & Antigen.bead == j) 
    
    model.ref = lmer(MFI ~ Virus.resp*lspline(DPSO, knots = c(100,200)) + Age + Sex  + (1|Code), 
                     data = subset_data) 
    
    
    bootstraps = case_bootstrap(model.ref, fixef, n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    
    bootstraps.at18 = case_bootstrap(model.ref, 
                                     
                                     function(m) {
                                       out <- emtrends(
                                         m,               # pass the refitted model here
                                         specs = ~ Virus.resp | DPSO, 
                                         var = "DPSO", 
                                         at = list(DPSO = c(20,100,200)),
                                         # (other emtrends arguments go here)
                                         lmerTest.limit = 19530, 
                                         pbkrtest.limit = 19530,
                                         data = subset_data
                                       )
                                       as_tibble(summary(out, infer = T, confint = T))
                                     }
                                     
                                     , n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    
    bootstraps.predictions = case_bootstrap(model.ref, function(model) {tibble(ggpredict(model, terms = c("DPSO [20, 100, 200, 550]", "Virus.resp"), ci.lvl = NA))}, n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    bootstraps.confint = bootstraps %>% confint()
    
    combo_name <- paste(i, j, k, sep = "_")
    named_result <- list(combo_name = list(boot.lmer = bootstraps, boot.confint = bootstraps.confint, predictions = bootstraps.predictions, ratesat18 = bootstraps.at18))
    names(named_result) <- combo_name
    named_result
    
  }


stopCluster(cl) # Stop the cluster
toc()

# do not print this
vir.results_list_strat_paral = unlist(vir.results_list_strat_paral, recursive = FALSE)



### --- PART 2: Result Aggregation &  Formatting --- ######################


# Initialize an empty data frame to combine all predictions
vir.combined_predictions_strat <- data.frame()
vir.combined_predictions_strat_observed <- data.frame()
vir.combined_predictions_strat_summary <- data.frame()
vir.combined_predictions_strat_timepoints <- data.frame()
vir.combined_predictions_strat_rates <- data.frame()

# Loop through the results_list to extract and combine predicted.lmer data frames
for (name in names(vir.results_list_strat_paral)) {
  # Extract the predicted.lmer data frame
  predicted_df <- vir.results_list_strat_paral[[name]]$boot.lmer$replicates 
  summary_df <- vir.results_list_strat_paral[[name]]$boot.lmer %>% confint()
  predicted.mg <- vir.results_list_strat_paral[[name]]$predictions$replicates 
  observed_df <- vir.results_list_strat_paral[[name]]$boot.lmer$observed 
  extracted_rates <- vir.results_list_strat_paral[[name]]$ratesat18$replicates
  
  # Add a column to preserve the combination name
  predicted_df <- predicted_df %>%
    mutate(Combination = name)
  # Add a column to preserve the combination name
  summary_df <- summary_df %>%
    mutate(Combination = name)
  
  predicted.mg <- predicted.mg %>%
    mutate(Combination = name)
  
  observed_df <- observed_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "coef") %>%
    rename(obs = 2) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    mutate(Combination = name) 
  
  
  extracted_rates <- extracted_rates %>%
    mutate(Combination = name)
  
  
  # Combine the data frames
  vir.combined_predictions_strat <- bind_rows(vir.combined_predictions_strat, predicted_df)
  
  vir.combined_predictions_strat_observed <- bind_rows(vir.combined_predictions_strat_observed, observed_df)
  
  # Combine the data frames
  vir.combined_predictions_strat_summary <- bind_rows(vir.combined_predictions_strat_summary, summary_df)
  
  vir.combined_predictions_strat_timepoints <- bind_rows(vir.combined_predictions_strat_timepoints, predicted.mg)
  
  vir.combined_predictions_strat_rates <- bind_rows(vir.combined_predictions_strat_rates, extracted_rates)
}

### --- PART 3: Statistical Inference & Slope Categorization ###########################

vir.combined_predictions_strat_summary2 = vir.combined_predictions_strat_summary %>%
  mutate(DPSO = case_when(grepl("))1", term) ~ 20, 
                          grepl("))2", term) ~ 100,
                          grepl("))3", term) ~ 200),
         Virus.resp = case_when(grepl("DV2_XR", term) ~ "DV2_XR",
                                grepl("DV4_XR", term) ~ "DV4_XR",
                                grepl("ZV_XR", term) ~ "ZV_XR"),
         interaction = grepl(":", term)) %>%
  filter(interaction == T,
         type == "norm") %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,3]) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR"))) %>%
  select(-term) %>%
  mutate(rel.slope.sig = case_when(  (lower < 0 & upper < 0) ~ "Sig. Waning",
                                     (lower > 0 & upper > 0) ~ "Sig. Rising",
                                     T ~ "Not Sig.")) %>%
  rename(rel.estimate = estimate, rel.lower = lower, rel.upper = upper)


### --- PART 4: Trend Summarization & Half-Life Calculation #############################

vir.combined_predictions_strat_rates.summary = vir.combined_predictions_strat_rates %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,3],
  ) %>%
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(DPSO.trend),
            sd = sd(DPSO.trend),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd,
            norm.upper = mean + qnorm(0.975) * sd
  ) %>%
  select(-mean, -sd, - n) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR", "DV4_XR","ZV_XR"))) %>%
  group_by(Combination, DPSO ,Isotype, IR, Antigen) %>%
  group_modify(~ ci_overlap_test(.x, filter_var = Virus.resp, reference_label = "HOMOTYPIC")) %>%
  mutate(overlap_with_ref = case_when(Virus.resp == "HOMOTYPIC" ~ "Ref.", T ~ as.character(overlap_with_ref ))) %>%
  mutate(slope.sig = case_when(  (norm.lower < 0 & norm.upper < 0) ~ "Sig. Waning",
                                 (norm.lower > 0 & norm.upper > 0) ~ "Sig. Rising",
                                 T ~ "Not Sig.")) %>%
  left_join(vir.combined_predictions_strat_summary2) %>%
  mutate(rel.slope.sig = case_when(is.na(rel.slope.sig) ~ "Ref.",
                                   T ~ rel.slope.sig),
         timepoint = case_when(DPSO == "20" ~ "Cv-3M",
                               DPSO == "100" ~ "3M-6M",
                               DPSO == "200" ~ "6M-18M"),
         hl.estimate = round((-1/norm.estimate)/365.25,2),
         hl.norm.lower = round((-1/norm.lower)/365.25,2),
         hl.norm.upper = round((-1/norm.upper)/365.25,2)) 


### --- PART 5: Marginal Prediction (MFI curves) ########################################

vir.combined_predictions_strat_timepoints2 = vir.combined_predictions_strat_timepoints %>%
  rename(Virus.resp = group, DPSO = x) %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         IR = str_split(Combination, pattern = "\\_", simplify = T)[,3]) %>%
  # left_join(cutoffs) %>%
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(predicted ),
            sd = sd(predicted ),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd,
            norm.upper = mean + qnorm(0.975) * sd
  ) %>%
  select(-mean, -sd, - n) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR"))) 






######################################################################3
## IR question  ------------------------------------------------------ 
######################################################################3

### PART 1: Parallel bootstrapping pipeline  ##########


tic()
n.simulations = 3 # 1000 simulations

# Set up parallel backend to use multiple processors
cl <- makeCluster(detectCores() - 1) # Leave one core free
registerDoParallel(cl)

# Loop through each combination in parallel
ir.results_list_strat_paral <- foreach(i = secondary_ab, .combine=list, .multicombine=TRUE, .packages=c('lme4', 'boot', 'ggeffects', 'dplyr', "lmeresampler", "lspline", "emmeans"), .errorhandling = 'pass') %:%
  foreach(j = antigen_bead, .combine='c', .errorhandling = 'pass') %:%
  foreach(l = virus_resp, .combine='c', .errorhandling = 'pass') %dopar%{
    
    subset_data <- subset(aim1.dataset2, `Secondary Ab` == i  & Antigen.bead == j & Virus.resp == l) 
    
    model.ref = lmer(MFI ~ IR*lspline(DPSO, knots = c(100,200)) + Age + Sex  + (1|Code), 
                     data = subset_data) 
    
    
    bootstraps = case_bootstrap(model.ref, fixef, n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    bootstraps.at18 = case_bootstrap(model.ref, 
                                     
                                     function(m) {
                                       out <- emtrends(
                                         m,               # pass the refitted model here
                                         specs = ~ IR | DPSO, 
                                         var = "DPSO", 
                                         at = list(DPSO = c(20,100,200)),
                                         # (other emtrends arguments go here)
                                         lmerTest.limit = 19530, 
                                         pbkrtest.limit = 19530,
                                         data = subset_data
                                       )
                                       as_tibble(summary(out, infer = T, confint = T))
                                     }
                                     
                                     , n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    
    bootstraps.predictions = case_bootstrap(model.ref, function(model) {tibble(ggpredict(model, terms = c("DPSO [20, 100, 200, 550]", "IR"), ci.lvl = NA))}, n.simulations, .refit = TRUE, resample = c(TRUE, F), orig_data = subset_data)
    
    bootstraps.confint = bootstraps %>% confint()
    
    combo_name <- paste(i, j, l, sep = "_")
    named_result <- list(combo_name = list(boot.lmer = bootstraps, boot.confint = bootstraps.confint, predictions = bootstraps.predictions, ratesat18 = bootstraps.at18))
    names(named_result) <- combo_name
    named_result
    
  }


stopCluster(cl) # Stop the cluster
toc()

# do not print this
ir.results_list_strat_paral = unlist(ir.results_list_strat_paral, recursive = FALSE)

### --- PART 2: Result Aggregation &  Formatting --- ######################


# Initialize an empty data frame to combine all predictions
ir.combined_predictions_strat <- data.frame()
ir.combined_predictions_strat_observed <- data.frame()
ir.combined_predictions_strat_summary <- data.frame()
ir.combined_predictions_strat_timepoints <- data.frame()
ir.combined_predictions_strat_rates <- data.frame()

# Loop through the results_list to extract and combine predicted.lmer data frames
for (name in names(ir.results_list_strat_paral)) {
  # Extract the predicted.lmer data frame
  predicted_df <- ir.results_list_strat_paral[[name]]$boot.lmer$replicates 
  summary_df <- ir.results_list_strat_paral[[name]]$boot.lmer %>% confint()
  predicted.mg <- ir.results_list_strat_paral[[name]]$predictions$replicates 
  observed_df <- ir.results_list_strat_paral[[name]]$boot.lmer$observed 
  extracted_rates <- ir.results_list_strat_paral[[name]]$ratesat18$replicates
  
  # Add a column to preserve the combination name
  predicted_df <- predicted_df %>%
    mutate(Combination = name)
  # Add a column to preserve the combination name
  summary_df <- summary_df %>%
    mutate(Combination = name)
  
  predicted.mg <- predicted.mg %>%
    mutate(Combination = name)
  
  observed_df <- observed_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "coef") %>%
    rename(obs = 2) %>%
    t() %>%
    janitor::row_to_names(row_number = 1) %>%
    as.data.frame() %>%
    mutate(Combination = name) 
  
  
  extracted_rates <- extracted_rates %>%
    mutate(Combination = name)
  
  
  # Combine the data frames
  ir.combined_predictions_strat <- bind_rows(ir.combined_predictions_strat, predicted_df)
  
  ir.combined_predictions_strat_observed <- bind_rows(ir.combined_predictions_strat_observed, observed_df)
  
  # Combine the data frames
  ir.combined_predictions_strat_summary <- bind_rows(ir.combined_predictions_strat_summary, summary_df)
  
  ir.combined_predictions_strat_timepoints <- bind_rows(ir.combined_predictions_strat_timepoints, predicted.mg)
  
  ir.combined_predictions_strat_rates <- bind_rows(ir.combined_predictions_strat_rates, extracted_rates)
}

### --- PART 3: Statistical Inference & Slope Categorization ###########################


ir.combined_predictions_strat_summary2 = ir.combined_predictions_strat_summary %>%
  mutate(DPSO = case_when(grepl("))1", term) ~ 20, 
                          grepl("))2", term) ~ 100,
                          grepl("))3", term) ~ 200),
         IR = case_when(grepl("IRS", term) ~ "S"),
         interaction = grepl(":", term)) %>%
  filter(interaction == T,
         type == "norm") %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                Virus.resp == "ZV" ~ "ZV_XR",
                                T ~ Virus.resp)) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR", "DV4_XR","ZV_XR"))) %>%
  select(-term) %>%
  mutate(rel.slope.sig = case_when(  (lower < 0 & upper < 0) ~ "Sig. Waning",
                                     (lower > 0 & upper > 0) ~ "Sig. Rising",
                                     T ~ "Not Sig.")) %>%
  rename(rel.estimate = estimate, rel.lower = lower, rel.upper = upper)


### --- PART 4: Trend Summarization & Half-Life Calculation #############################


ir.combined_predictions_strat_rates.summary = ir.combined_predictions_strat_rates %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                Virus.resp == "ZV" ~ "ZV_XR",
                                T ~ Virus.resp)) %>%
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(DPSO.trend),
            sd = sd(DPSO.trend),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd,
            norm.upper = mean + qnorm(0.975) * sd
  ) %>%
  select(-mean, -sd, - n) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR", "DV4_XR","ZV_XR"))) %>%
  group_by(Combination, DPSO ,Isotype, Virus.resp, Antigen) %>%
  group_modify(~ ci_overlap_test(.x, filter_var = IR, reference_label = "P")) %>%
  mutate(overlap_with_ref = case_when(IR == "P" ~ "Ref.", T ~ as.character(overlap_with_ref ))) %>%
  mutate(slope.sig = case_when(  (norm.lower < 0 & norm.upper < 0) ~ "Sig. Waning",
                                 (norm.lower > 0 & norm.upper > 0) ~ "Sig. Rising",
                                 T ~ "Not Sig.")) %>%
  left_join(ir.combined_predictions_strat_summary2) %>%
  mutate(rel.slope.sig = case_when(is.na(rel.slope.sig) ~ "Ref.",
                                   T ~ rel.slope.sig),
         timepoint = case_when(DPSO == "20" ~ "Cv-3M",
                               DPSO == "100" ~ "3M-6M",
                               DPSO == "200" ~ "6M-18M"),
         hl.estimate = round((-1/norm.estimate)/365.25,2),
         hl.norm.lower = round((-1/norm.lower)/365.25,2),
         hl.norm.upper = round((-1/norm.upper)/365.25,2)) 


### --- PART 5: Marginal Prediction (MFI curves) ########################################

ir.combined_predictions_strat_timepoints2 = ir.combined_predictions_strat_timepoints %>%
  rename(IR = group, DPSO = x) %>%
  mutate(Isotype = str_split(Combination, pattern = "\\_", simplify = T)[,1],
         Antigen = str_split(Combination, pattern = "\\_", simplify = T)[,2],
         Virus.resp = str_split(Combination, pattern = "\\_", simplify = T)[,3],
         Virus.resp = case_when(Virus.resp == "DV2" ~ "DV2_XR",
                                Virus.resp == "DV4" ~ "DV4_XR",
                                Virus.resp == "ZV" ~ "ZV_XR",
                                T ~ Virus.resp)) %>%
  group_by(DPSO, Antigen, Combination, Isotype, IR, Virus.resp) %>%
  summarise(mean = mean(predicted ),
            sd = sd(predicted ),
            n = n.simulations,
            norm.estimate = mean,
            norm.lower = mean - qnorm(0.975) * sd,
            norm.upper = mean + qnorm(0.975) * sd
  ) %>%
  select(-mean, -sd, - n) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG", "IgG1", "IgG2", "IgG3", "IgG4", "IgA", "IgM")),
         IR = factor(IR, levels = c("P", "S")),
         Antigen = factor(Antigen, levels = c("E", "EDIII", "NS1")),
         Virus.resp = factor(Virus.resp, levels = c("HOMOTYPIC", "DV2_XR","DV4_XR", "ZV_XR"))) 



aim1.dataset2.forind = aim1.dataset2 %>%
  select(Code, DPSO, IR, Isotype = `Secondary Ab`, Virus.resp, Antigen = Antigen.bead, norm.estimate = MFI, cutoff) %>%
  mutate(norm.estimate = norm.estimate - cutoff)



#######################################################################################3
# Bayesian model --------------------------------------------------------------------
#######################################################################################3

make_spline_basis <- function(dpso_vector, knots = c(100, 200)) {
  
  # Generate the linear spline matrix
  B <- lspline(dpso_vector, knots = knots)
  
  # Convert to tibble and apply  specific naming convention
  B_df <- as_tibble(B)
  colnames(B_df) <- c("DPSO_0_100", "DPSO_100_200", "DPSO_200p")
  
  return(B_df)
}

### --- PART 1: Data preparation --- ########################

iter = 40
warmup = 21
# iter = 1500 for final model
# warmup = 500 for final model

# 1. Create Linear Spline Basis for DPSO

B <- lspline(aim1.dataset2$DPSO, knots = c(100, 200))
colnames(B) <- c("DPSO_0_100", "DPSO_100_200", "DPSO_200p")



dat <- bind_cols(aim1.dataset2, as.data.frame(B))

dat$SecondaryAb <-  dat$`Secondary Ab`

# 'Group4' combines all "antibody compartments" to allow for varying slopes/intercepts
dat$Group4 <- interaction(
  dat$Virus.resp, dat$Antigen.bead, dat$SecondaryAb, dat$IR,
  drop = TRUE
)

# make sure I can identify "make" and pathway is ok for modeling with C language integration
Sys.setenv(PATH = paste(
  "C:/rtools44/usr/bin",
  "C:/rtools44/mingw64/bin",
  Sys.getenv("PATH"),
  sep = ";"
))

Sys.which("make")
check_cmdstan_toolchain(fix = TRUE)

### --- PART 2: The model ########################
# Formula Logic:
# - Population effects: Interaction of spline segments by antibody compartments (IR Antigen.bead Virus.resp SecondaryAb)
# - Random slopes: (Spline terms | Group4) allows kinetics to vary by antibody compartments
# - Random intercepts: (1 | Code) accounts for individual baseline variation
# - Nested intercepts: (1 | Code:Group4) accounts for person-specific antibody compartment baseline variation

set.seed(123)

fit_pragmatic2B <- brm(
  MFI ~ (DPSO_0_100 + DPSO_100_200 + DPSO_200p) * (IR + Antigen.bead + Virus.resp + SecondaryAb) + Age + Sex + (DPSO_0_100 + DPSO_100_200 + DPSO_200p | Group4) + (1 | Code:Group4) + (1 | Code) ,
  data    = dat,
  family  = gaussian(),
  prior   = c(
    prior(normal(0, 1), class = "b"), # Weakly informative fixed effects
    prior(student_t(3, 0, 1), class = "sd", group = "Group4"), # Robust priors for variance
    prior(student_t(3, 0, 1), class = "sd", group = "Code:Group4"),
    prior(student_t(3, 0, 1), class = "sd", group = "Code"),
    prior(student_t(3, 0, 1), class = "sigma")
  ),
  chains  = 4,
  cores   = 4,        
  iter    = iter,
  warmup  = warmup,
  threads = threading(16),
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97, max_treedepth = 18)
)

--HERE
# 4. Model Diagnostics & Performance
R2_IQR = quantile(brms::bayes_R2(fit_pragmatic2B), c(0.25, 0.75))
pp_check(fit_pragmatic2B, type = "scatter_avg")
pp_check(fit_pragmatic2B, type = "dens_overlay", ndraws = 10)
hist(brms::bayes_R2(fit_pragmatic2B))

# plot(fit_pragmatic2B)

# check Rsquares
r2s = brms::bayes_R2(fit_pragmatic2B)

# model summary
broom::tidy(fit_pragmatic2B)

### --- PART 3: Posterior Predictive Slope Extraction ######################
# Compute Posterior Slopes for Specific Time Windows
# This function simulates antibody decay/growth rates between two timepoints (t1, t2)
# by drawing from the posterior distribution of the fitted Bayesian model.

# @param t1 Start day (DPSO)
# @param t2 End day (DPSO)
# @param dat The original training dataset (for metadata lookup)
# @param fit The fitted brms object
# @param codes Character vector of specific participant Codes to represent 'typical' variation in one primary and secondary individual

compute_posterior_slopes <- function(t1, t2,
                                     dat,
                                     fit,
                                     codes = c("1007", "1008")) {
  
  # 1. Map Group4 identifiers back their corresponding categories
  group4_lookup <- dat %>%
    transmute(
      Code = as.factor(Code),
      Group4,
      Virus.resp,
      Antigen.bead,
      SecondaryAb,
      IR
    ) %>%
    distinct()
  
  # 2. Filter for representative participant IDs
  lookup_ids <- group4_lookup %>%
    filter(as.character(Code) %in% codes)
  
  # 3. Construct Prediction Grid (Newdata)
  # Generates 2 rows per group (one for t1, one for t2) with standardized Age/Sex
  
  newdata <- lookup_ids %>%
    distinct(Code, Group4, Virus.resp, Antigen.bead, SecondaryAb, IR) %>%
    tidyr::crossing(DPSO = c(t1, t2)) %>%
    mutate(
      Age = median(dat$Age, na.rm = TRUE),
      Sex = names(sort(table(dat$Sex), decreasing = TRUE))[1]
    ) %>%
    bind_cols(make_spline_basis(.$DPSO)) %>% # Apply the same spline transformation used in fitting
    mutate(.row = row_number())
  
  # 4. Extract Posterior Draws
  # re_formula = NULL ensures individual-level random effects are included
  eta_draws <- posterior_linpred(
    fit,
    newdata    = newdata,
    re_formula = NULL,
    transform  = FALSE
  )
  
  # 5. Reshape Draws to Long Format
  eta_long <- as_tibble(eta_draws) %>%
    mutate(draw = row_number()) %>%
    pivot_longer(-draw, names_to = ".col", values_to = "eta") %>%
    mutate(.row = readr::parse_number(.col)) %>%
    select(draw, .row, eta)
  
  # 6. Calculate Slopes (Rise/Run)
  slopes <- eta_long %>%
    left_join(newdata %>% select(.row, Code, Group4, DPSO), by = ".row") %>%
    group_by(Code, Group4, draw) %>%
    summarise(
      slope = (eta[DPSO == t2] - eta[DPSO == t1]) / (t2 - t1),
      .groups = "drop"
    )
  
  # 7. Summarize the Posterior Distribution
  slopes_summary <- slopes %>%
    group_by(Group4) %>%
    summarise(
      median = median(slope),
      q2.5   = quantile(slope, 0.025),
      q97.5  = quantile(slope, 0.975),
      .groups = "drop"
    ) %>%
    # Split the dot-separated Group4 string back into tidy columns
    separate(
      Group4,
      into   = c("Virus.resp","Antigen.bead","SecondaryAb","IR"),
      sep    = "\\.",
      remove = FALSE
    )
  
  return(slopes_summary)
}


# Calculate kinetics for the 'Early' phase (Day 20 to 90)
slopes_Cv_3M <- compute_posterior_slopes(
  t1  = 20,
  t2  = 90,
  dat = dat,
  fit = fit_pragmatic2B
) %>%
  mutate(DPSO = "Cv-3M")

# Calculate kinetics for the 6-18M period (180–550 days)
slopes_6M_18M <- compute_posterior_slopes(
  t1  = 180,
  t2  = 550,
  dat = dat,
  fit = fit_pragmatic2B
) %>%
  mutate(DPSO = "6M-18M")

# Combine the early-phase (Cv-3M) and late-phase (6M-18M) posterior estimates
slopes.combined = bind_rows(slopes_Cv_3M, slopes_6M_18M) %>%
  relocate(DPSO, 1)

### --- PART 4: Significance & Overlap Testing ---##########################

# Apply the ci_overlap_test iteratively across every experimental antibody compartments

abfactorlevels = c("IgG", "IgA", "IgM", "IgG1", "IgG2", "IgG3", "IgG4")

abfactorlevelsfun = function(x)(x %>% mutate(SecondaryAb  = factor(SecondaryAb, levels = abfactorlevels)))

# 2. Point-in-Interval Overlap Function
# Checks if the reference median falls within the 95% CrI of the comparison group.
ci_overlap_test <- function(data, filter_var, reference_label) {
  
  filter_var <- enquo(filter_var)
  overlap_col <- paste0(quo_name(filter_var), "_overlap_with_ref")
  # Identify the baseline value (median slope of the reference group)
  ref_row <- data %>% filter(!!filter_var == reference_label)
  
  ref <- ref_row$median
  
  data %>%
    mutate(
      !!overlap_col := case_when(
        !!filter_var == reference_label ~ "Ref.",
        q2.5 <= ref & q97.5 >= ref      ~ "Not Sig. Diff",
        TRUE                            ~ "Sig. Diff"
      )
    )
}



# 3. Comprehensive Slope Comparisons
slope_comparisons = slopes.combined %>%
  mutate(slope.sig =
           case_when(
             (q2.5 < 0 & q97.5 < 0) ~ "Sig. Waning",
             (q2.5 > 0 & q97.5 > 0) ~ "Sig. Rising",
             T ~ "Not Sig."),
         
         hl.estimate = round((-1/median)/365.25,2),
         hl.q2.5 = round((-1/q2.5)/365.25,2),
         hl.q97.5 = round((-1/q97.5)/365.25,2)
  ) %>%
  # 4. Iterative Reference Testing (The "Correction" Layers)
  # Compare Infection History: Primary (P) vs Secondary (S)
  group_by(DPSO, SecondaryAb, Virus.resp, Antigen.bead) %>%
  group_modify(~ ci_overlap_test(.x, filter_var = IR, reference_label = "P")) %>%
  ungroup() %>%
  # Compare Antigen Targets: Reference = Envelope (E)
  group_by(DPSO, SecondaryAb, IR, Virus.resp) %>%
  group_modify(~ ci_overlap_test(.x, filter_var = Antigen.bead, reference_label = "E")) %>%
  ungroup() %>%
  # Compare Viral Cross-Reactivity: Reference = HOMOTYPIC
  group_by(DPSO, SecondaryAb, IR, Antigen.bead) %>%
  group_modify(~ ci_overlap_test(.x, filter_var = Virus.resp, reference_label = "HOMOTYPIC")) %>%
  ungroup() %>%
  group_by(DPSO, Virus.resp, IR, Antigen.bead) %>%
  # Compare Antibody Isotypes: Reference = Total IgG
  group_modify(~ ci_overlap_test(.x, filter_var = SecondaryAb, reference_label = "IgG")) %>%
  ungroup() %>%
  arrange(DPSO, Virus.resp, IR, Antigen.bead, SecondaryAb)


### --- PART 5: Longitudinal Trajectory & Cutoff Normalization #######################

# 1. Define the Prediction Grid
# predicted levels at these specific days
dpsogrid2 = c(20,90,180, 540)

# 2. Extract unique experimental strata (Group4)
group4_levels <- dat %>%
  distinct(Virus.resp, Antigen.bead, IR, SecondaryAb) %>%
  mutate(Group4 = interaction(
    Virus.resp, Antigen.bead, SecondaryAb, IR,
    drop = TRUE
  ))


# 3. Create 'Newdata' for Prediction
newdata <- group4_levels %>%
  tidyr::crossing(DPSO = dpsogrid2) %>%
  mutate(
    Age = median(dat$Age, na.rm = TRUE),
    Sex = names(sort(table(dat$Sex), decreasing = TRUE))[1]
  ) %>%
  bind_cols(make_spline_basis(.$DPSO)) %>%
  # Make sure Group4 has same factor levels as in the original data
  mutate(Group4 = factor(Group4, levels = levels(dat$Group4)))

newdata2 <- newdata %>% mutate(.row = row_number())


# 4. Generate Posterior Predictions (Fitted Values)
# re_formula = NULL would include all random effects; 
# Using specific Group4 random effects gives the "typical" trajectory for each assay.

pred_traj <- fitted(
  fit_pragmatic2B,
  newdata = newdata2,
  # when running 1000 simulations run the line below
  re_formula = ~(DPSO_0_100 + DPSO_100_200 + DPSO_200p | Group4),
  
  summary = TRUE
)

pred_tbl <- as_tibble(pred_traj) %>%
  mutate(.row = row_number())

# 5. extract Background Cutoffs
cutoffs = aim1.dataset2 %>%
  select(SecondaryAb = `Secondary Ab`, Virus.resp , Antigen.bead, cutoff) %>%
  group_by(SecondaryAb, Virus.resp, Antigen.bead) %>%
  summarise(cutoff = max(cutoff))

# 6. Final Normalization
pred_traj_df <- left_join(newdata2, pred_tbl, by = ".row") %>%
  left_join(cutoffs) %>%
  mutate(norm.estimate = Estimate - cutoff,
         norm.lower = Q2.5 - cutoff,
         norm.upper = Q97.5 - cutoff)



###  Individual-Level Fold Change & Half-Life Calculation #######################

# 1. Prepare Wide Dataset
# Reshape data so each row represents one unique assay per participant,
# with timepoints and MFI values spread across columns.

aim1.dataset2.hl = aim1.dataset2 %>%
  dplyr::select(Code, Timepoint, x = DPSO, Virus.resp, IR, Antigen = Antigen.bead, Isotype = `Secondary Ab`, MFI) %>%
  mutate(netMFI.Naive = MFI) %>%
  pivot_wider(id_cols = c(Code,  Virus.resp, IR, Antigen, Isotype), names_from = c("Timepoint"), values_from = c("x","netMFI.Naive")) 

# 2. Define Kinetic Function
# Calculate Rate of Change (Slope)
# @param time1,time2 Days post-symptom onset
# @param conc1,conc2 MFI values at those times

# Function to calculate fold change per year
calculate_fold_change_per_year <- function(time1, conc1, time2, conc2) {
  delta_t <- (time2 - time1)  
  fold_change <- (conc2 - conc1) / delta_t
  return(fold_change)
}


# 3. Compute Individual Slopes and Half-Lives
# Iterate row-by-row to calculate kinetics for each visit interval

aim1.dataset2.hl <- aim1.dataset2.hl %>%
  rowwise() %>%
  mutate(
    # A. Calculate the daily rate of change for each interval
    fchange_Cv_3 = calculate_fold_change_per_year(x_Cv, netMFI.Naive_Cv, x_3, netMFI.Naive_3),
    fchange_3_6 = calculate_fold_change_per_year(x_3, netMFI.Naive_3, x_6, netMFI.Naive_6),
    fchange_6_18 = calculate_fold_change_per_year(x_6, netMFI.Naive_6, x_18, netMFI.Naive_18),
    # B. Convert rate of change to Half-Life (t1/2) in Years
    half_life_Cv_3 = (-1/(fchange_Cv_3))/365.25,
    half_life_3_6 = (-1/fchange_3_6)/365.25,
    half_life_6_18 = (-1/fchange_6_18)/365.25)


#############################################################################3
# Bivariate Statistical Comparisons --------------------------------------------------------
#############################################################################3

### --- PART 1: Data preparation --- ########################

# Isolate specific kinetic outcomes for univariate testing

univariate.dataset = aim1.dataset2.hl %>%
  select(Code:Isotype, fchange_Cv_3, fchange_6_18 )

### --- PART 2: Unpaired Group Comparisons --- ########################

# ---------------------------------------------------------3
# 2. Function: Unpaired Group Comparisons (e.g., Primary vs Secondary IR)
# Uses the Wilcoxon Rank-Sum test (Mann-Whitney U) to see if two 
# independent groups have different kinetic rates.
# ---------------------------------------------------------3

run_ir_wilcox <- function(dat, outcome,
                          strat_vars = c("Virus.resp", "Antigen", "Isotype"),
                          group_var  = "IR",
                          p_adj      = "bonferroni") {
  
  
  f <- as.formula(paste0(outcome, " ~ ", group_var))
  
  tests <- dat %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[group_var]])) %>%
    group_by(across(all_of(strat_vars))) %>%
    wilcox_test(f, detailed = TRUE) %>%
    rstatix::adjust_pvalue(method = p_adj) %>%
    ungroup()
  
  eff <- dat %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[group_var]])) %>%
    group_by(across(all_of(strat_vars))) %>%
    wilcox_effsize(f) %>%
    ungroup()
  
  tests %>%
    left_join(
      eff %>% select(all_of(strat_vars), group1, group2, effsize, magnitude),
      by = c(strat_vars, "group1", "group2")
    )
}


### --- PART 3: Paired Group Comparisons --- ########################

run_paired_ref_wilcox <- function(dat, outcome,
                                  id_var    = "Code",
                                  group_var = "Antigen",     # or "Virus.resp"
                                  ref_group,
                                  strat_vars,
                                  p_adj = "bonferroni") {
  
  # Define the statistical formula (e.g., fchange_6_18 ~ IR)
  f <- as.formula(paste0(outcome, " ~ ", group_var))
  
  dat2 <- dat %>%
    select(all_of(c(id_var, strat_vars, group_var, outcome))) %>%
    filter(!is.na(.data[[outcome]]), !is.na(.data[[group_var]]))
  
  # Perform the Wilcoxon test across all specified strata
  tests <- dat2 %>%
    group_by(across(all_of(strat_vars))) %>%
    pairwise_wilcox_test(
      f,
      paired = TRUE,
      detailed = T,
      ref.group = ref_group,
      p.adjust.method = p_adj
    ) %>%
    ungroup()
  
  # Calculate Effect Size (r) to determine the magnitude of the difference
  eff <- dat2 %>%
    group_by(across(all_of(strat_vars))) %>%
    wilcox_effsize(
      f,
      paired = TRUE,
      ref.group = ref_group,
      p.adjust.method = "none"
    ) %>%
    ungroup()
  
  # Join p-values with effect sizes for a comprehensive result table
  tests %>%
    left_join(
      eff %>% select(all_of(strat_vars), group1, group2, effsize, magnitude),
      by = c(strat_vars, "group1", "group2")
    )
}



abfactorlevels = c("IgG", "IgA", "IgM", "IgG1", "IgG2", "IgG3", "IgG4")

abfactorlevelsfun = function(x)(x %>% mutate(SecondaryAb  = factor(Isotype, levels = abfactorlevels)))

### --- PART 3: Comparison by IR --- ########################

# Do Primary (P) vs. Secondary (S) individuals have different waning rates?

IRtestsCV_3M_results <- run_ir_wilcox(
  dat     = univariate.dataset,
  outcome = "fchange_Cv_3"
)

IRtests6_18_results <- run_ir_wilcox(
  dat     = univariate.dataset,
  outcome = "fchange_6_18"
)


### --- PART 4: Comparison by Antigen --- ########################
# Comparison of Antigen  (Reference = E)
# Does EDIII or NS1 wane faster/slower than the Envelope (E) protein?
# Uses a PAIRED test because each 'Code' has measurements for all antigens.
antigenCv_3_results <- run_paired_ref_wilcox(
  dat       = univariate.dataset,
  outcome   = "fchange_Cv_3",
  id_var    = "Code",
  group_var = "Antigen",
  ref_group = "E",
  strat_vars = c("IR","Virus.resp", "Isotype")
)

antigen6_18_results <- run_paired_ref_wilcox(
  dat       = univariate.dataset,
  outcome   = "fchange_6_18",
  id_var    = "Code",
  group_var = "Antigen",
  ref_group = "E",
  strat_vars = c("IR","Virus.resp", "Isotype")
)

### --- PART 4: Comparison by Viral response --- ########################
# Comparison of Viral Response (Reference = HOMOTYPIC)
# Do cross-reactive responses (DV2_XR, ZV_XR, etc.) wane differently than 
# the homotypic response? Also uses a PAIRED test.

virusCv_3_results <- run_paired_ref_wilcox(
  dat       = univariate.dataset,
  outcome   = "fchange_Cv_3",
  id_var    = "Code",
  group_var = "Virus.resp",
  ref_group = "HOMOTYPIC",
  strat_vars = c("IR","Antigen", "Isotype")
)

virus6_18_results <- run_paired_ref_wilcox(
  dat       = univariate.dataset,
  outcome   = "fchange_6_18",
  id_var    = "Code",
  group_var = "Virus.resp",
  ref_group = "HOMOTYPIC",
  strat_vars = c("IR","Antigen", "Isotype")
)
#############################################################################
# MS plots
#############################################################################

## Palettes ------------------------------------------------------------------

ir_palette <- c("P" = "#C08619",
                "S" = "#A9A9A9")

ir_palette2 <- c("P" = "#B07B17",
                 "S" = "#4B4B4B")

ir_palette_DVXR <- c("P"       = "#2DA2D4",
                     "S"       = "#A9A9A9",
                     "P_DV4_XR" = "#227FA6",
                     "S_DV4_XR" = "#A9A9A9",
                     "P_DV2_XR" = "#6CC08C")

ir_palette2_DVXR <- c("P"       = "#227FA6",
                      "S"       = "#4B4B4B",
                      "P_DV4_XR" = "#227FA6",
                      "S_DV4_XR" = "#4B4B4B",
                      "P_DV2_XR" = "#6CC08C")

virus_resp_palette <- c("HOMOTYPIC" = "#C08619",
                        "DV_XR"    = "#2DA2D4",
                        "ZV_XR"    = "#8A8A8A")

virus_resp_palette3 <- c("HOMOTYPIC" = "#C08619",
                         "DV4_XR"   = "#2DA2D4",
                         "DV2_XR"   = "#6CC08C",
                         "ZV_XR"    = "#8A8A8A")

virus_resp_palette2 <- c("HOMOTYPIC" = "black",
                         "DV_XR"    = "black",
                         "ZV_XR"    = "black",
                         "DV2_XR"   = "black",
                         "DV4_XR"   = "black")

antigen_palette <- c("E"    = "#013872",
                     "EDIII" = "#8D9000",
                     "NS1"  = "#C81671")

antigen_palette3 <- c("E"    = "#005686",
                      "EDIII" = "#7C7E1D",
                      "NS1"  = "#C81671")

antigen_palette2 <- c("E"    = "black",
                      "EDIII" = "black",
                      "NS1"  = "black")

isotype_palette <- c("IgG"  = "#3B0192",
                     "IgG1" = "#4D1368",
                     "IgG2" = "#962D63",
                     "IgG3" = "#F54421",
                     "IgG4" = "#FD963C",
                     "IgA"  = "#C6A517",
                     "IgM"  = "#038861")

isotype_palette2 <- c("IgG"  = "black",
                      "IgG1" = "black",
                      "IgG2" = "black",
                      "IgG3" = "black",
                      "IgG4" = "black",
                      "IgA"  = "black",
                      "IgM"  = "black")

ir_palette_corr <- c("P_HOMOTYPIC" = "#B07B17",
                     "S_HOMOTYPIC" = "#4B4B4B",
                     "P_DV4_XR"   = "#227FA6",
                     "S_DV4_XR"   = "#4B4B4B")

hl_palette <- c(rise   = "#C87105",
                wane   = "#52C1A0",
                stable = "#30433B")

## Functions -----------------------------------------------------------------

formatter <- function() {
  function(y) round((-1/y)/365.25, 2)
}

formatter.years <- function() {
  function(y) round((-1/y), 2)
}

## Fig 1 - Homotypic E -------------------------------------------------------

cutoffs2 <- cutoffs %>%
  rename(Isotype = SecondaryAb, Antigen = Antigen.bead)

### 1A - Curves by IR and Isotype --------------------------------------------

ir.combined_predictions_strat_timepoints2 %>%
  left_join(cutoffs2) %>%
  mutate(
    norm.estimate = norm.estimate - cutoff,
    norm.lower    = norm.lower    - cutoff,
    norm.upper    = norm.upper    - cutoff,
    cutoff        = 0
  ) %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
  filter(Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
         Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  ggplot(aes(x = DPSO, y = norm.estimate)) +
  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = IR), alpha = 0.6) +
  geom_line(aes(linetype = IR), size = 1) +
  scale_fill_manual(values = ir_palette) +
  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

### 1B - Half Life Cv-3M -----------------------------------------------------

bracket_data <- ir.combined_predictions_strat_rates.summary %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
  filter(timepoint == "Cv-3M", Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  group_by(Isotype) %>%
  filter(any(overlap_with_ref == FALSE)) %>%
  summarize(
    xmin  = 1 - 0.2,
    xmax  = 1 + 0.2,
    y_top = max(norm.upper[IR %in% c("P","S")]) + 0.01,
    .groups = "drop"
  )

ir.combined_predictions_strat_rates.summary %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
  filter(timepoint == "Cv-3M",
         Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
         Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
             color = IR, alpha = slope.sig, linetype = IR)) +
  geom_linerange(size = 0.3, position = position_dodge(1)) +
  geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
  scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = c(0.4, 0.5)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
  scale_color_manual(values = ir_palette2) +
  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
  labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
  geom_segment(data = bracket_data,
               aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
               inherit.aes = FALSE, size = 0.2) +
  geom_segment(data = bracket_data,
               aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
               inherit.aes = FALSE, size = 0.2) +
  geom_segment(data = bracket_data,
               aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
               inherit.aes = FALSE, size = 0.2) +
  geom_text(data = bracket_data,
            aes(x = (xmin + xmax)/2, y = y_top + 0.005, label = "*"),
            inherit.aes = FALSE, size = 3) +
  theme_classic() +
  theme(panel.grid = element_blank())

### 1C - Half Life 6-18M -----------------------------------------------------

bracket_data <- ir.combined_predictions_strat_rates.summary %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
  filter(timepoint == "6M-18M", Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  group_by(Isotype) %>%
  filter(any(overlap_with_ref == FALSE)) %>%
  summarize(
    xmin  = 1 - 0.2,
    xmax  = 1 + 0.2,
    y_top = max(norm.upper[IR %in% c("P","S")]) + 0.001,
    .groups = "drop"
  )

ir.combined_predictions_strat_rates.summary %>%
  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
  filter(timepoint == "6M-18M",
         Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
         Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
             color = IR, alpha = slope.sig, linetype = IR)) +
  geom_linerange(size = 0.3, position = position_dodge(1)) +
  geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
  scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = c(0.4, 0.5)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
  scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                     guide = guide_legend(override.aes = list(fill = "black"))) +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
  scale_color_manual(values = ir_palette2) +
  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
  labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
  geom_segment(data = bracket_data,
               aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
               inherit.aes = FALSE, size = 0.2) +
  geom_segment(data = bracket_data,
               aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.0005),
               inherit.aes = FALSE, size = 0.2) +
  geom_segment(data = bracket_data,
               aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.0005),
               inherit.aes = FALSE, size = 0.2) +
  geom_text(data = bracket_data,
            aes(x = (xmin + xmax)/2, y = y_top + 0.0005, label = "*"),
            inherit.aes = FALSE, size = 3) +
  theme_classic() +
  theme(panel.grid = element_blank())

### 1D-a - Curve Merge IgG-A-M by IR ----------------------------------------

ir.combined_predictions_strat_timepoints2 %>%
  left_join(cutoffs2) %>%
  mutate(
    norm.estimate = norm.estimate - cutoff,
    norm.lower    = norm.lower    - cutoff,
    norm.upper    = norm.upper    - cutoff,
    cutoff        = 0
  ) %>%
  filter(Isotype %in% c("IgG","IgA","IgM"), Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  ggplot(aes(x = DPSO, y = norm.estimate)) +
  facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
  geom_line(aes(linetype = IR, color = Isotype), size = 1) +
  scale_color_manual(values = isotype_palette2) +
  scale_fill_manual(values = isotype_palette) +
  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))

### 1D-b - Curve merge IgG1-4 by IR -----------------------------------------

ir.combined_predictions_strat_timepoints2 %>%
  left_join(cutoffs2) %>%
  mutate(
    norm.estimate = norm.estimate - cutoff,
    norm.lower    = norm.lower    - cutoff,
    norm.upper    = norm.upper    - cutoff,
    cutoff        = 0
  ) %>%
  filter(Isotype %in% c("IgG1","IgG2","IgG3","IgG4"), Virus.resp == "HOMOTYPIC", Antigen == "E") %>%
  ggplot(aes(x = DPSO, y = norm.estimate)) +
  facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
  geom_line(aes(linetype = IR, color = Isotype), size = 1) +
  scale_color_manual(values = isotype_palette2) +
  scale_fill_manual(values = isotype_palette) +
  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1))
        
### 1E - Seropositivity at 18M -----------------------------------------------
        
      aim1.dataset2 %>%
          mutate(netMFI.Naive2 = MFI - cutoff) %>%
          filter(Virus.resp == "HOMOTYPIC", Timepoint == "18", Antigen.bead == "E") %>%
          mutate(
            seropositivity = case_when(netMFI.Naive2 > 0 ~ "Seropositive",
                                       netMFI.Naive2 <= 0 ~ "Seronegative"),
            `Secondary Ab` = factor(`Secondary Ab`, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))
          ) %>%
          count(seropositivity, IR, `Secondary Ab`, Antigen.bead, Timepoint, Virus.resp) %>%
          ungroup() %>%
          ggplot(aes(y = n, x = Antigen.bead, alpha = seropositivity, fill = `Secondary Ab`)) +
          geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.3, color = "grey60") +
          facet_grid(vars(IR), vars(`Secondary Ab`)) +
          scale_fill_manual(values = isotype_palette) +
          scale_alpha_manual(values = c("Seropositive" = 0.75, "Seronegative" = 0.15)) +
          labs(y = "Seroprevalence at\n18 months post infection") +
          geom_col(position = "fill") +
          theme_classic() +
          theme(axis.text.x  = element_text(angle = 90),
                panel.spacing.x = unit(0.2, "lines"),
                panel.spacing.y = unit(0.2, "lines"))
        
        ## Fig 2 - Homotypic NS1 -----------------------------------------------------
        
        ### 2A - Curves by isotype ---------------------------------------------------
        
        ir.combined_predictions_strat_timepoints2 %>%
          left_join(cutoffs2) %>%
          mutate(
            norm.estimate = norm.estimate - cutoff,
            norm.lower    = norm.lower    - cutoff,
            norm.upper    = norm.upper    - cutoff,
            cutoff        = 0
          ) %>%
          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
          filter(Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                 Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          ggplot(aes(x = DPSO, y = norm.estimate)) +
          facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = IR), alpha = 0.6) +
          geom_line(aes(linetype = IR), size = 1) +
          scale_fill_manual(values = ir_palette) +
          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
          scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
          geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
          theme_classic() +
          theme(axis.text.x = element_text(hjust = 1))
        
        ### 2B - Half Life Cv-3M -----------------------------------------------------
        
        bracket_data <- ir.combined_predictions_strat_rates.summary %>%
          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
          filter(timepoint == "Cv-3M", Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          group_by(Isotype) %>%
          filter(any(overlap_with_ref == FALSE)) %>%
          summarize(
            xmin  = 1 - 0.2,
            xmax  = 1 + 0.2,
            y_top = max(norm.upper[IR %in% c("P","S")]) + 0.01,
            .groups = "drop"
          )
        
        ir.combined_predictions_strat_rates.summary %>%
          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
          filter(timepoint == "Cv-3M",
                 Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                 Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                     color = IR, alpha = slope.sig, linetype = IR)) +
          geom_linerange(size = 0.3, position = position_dodge(1)) +
          geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
          facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
          scale_x_discrete(expand = c(0.4, 0.5)) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                   fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
          scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                             guide = guide_legend(override.aes = list(fill = "black"))) +
          scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
          scale_color_manual(values = ir_palette2) +
          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
          labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
          geom_segment(data = bracket_data,
                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                       inherit.aes = FALSE, size = 0.2) +
          geom_segment(data = bracket_data,
                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
                       inherit.aes = FALSE, size = 0.2) +
          geom_segment(data = bracket_data,
                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
                       inherit.aes = FALSE, size = 0.2) +
          geom_text(data = bracket_data,
                    aes(x = (xmin + xmax)/2, y = y_top + 0.005, label = "*"),
                    inherit.aes = FALSE, size = 3) +
          theme_classic() +
          theme(panel.grid = element_blank())
        
        ### 2C - Half Life 6-18M -----------------------------------------------------
        
        bracket_data <- ir.combined_predictions_strat_rates.summary %>%
          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
          filter(timepoint == "6M-18M", Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          group_by(Isotype) %>%
          filter(any(overlap_with_ref == FALSE)) %>%
          summarize(
            xmin  = 1 - 0.2,
            xmax  = 1 + 0.2,
            y_top = max(norm.upper[IR %in% c("P","S")]) + 0.001,
            .groups = "drop"
          )
        
        ir.combined_predictions_strat_rates.summary %>%
          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
          filter(timepoint == "6M-18M",
                 Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                 Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                     color = IR, alpha = slope.sig, linetype = IR)) +
          geom_linerange(size = 0.3, position = position_dodge(1)) +
          geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
          facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
          scale_x_discrete(expand = c(0.4, 0.5)) +
          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                   fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
          scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                             guide = guide_legend(override.aes = list(fill = "black"))) +
          scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
          scale_color_manual(values = ir_palette2) +
          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
          labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
          geom_segment(data = bracket_data,
                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                       inherit.aes = FALSE, size = 0.2) +
          geom_segment(data = bracket_data,
                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.0005),
                       inherit.aes = FALSE, size = 0.2) +
          geom_segment(data = bracket_data,
                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.0005),
                       inherit.aes = FALSE, size = 0.2) +
          geom_text(data = bracket_data,
                    aes(x = (xmin + xmax)/2, y = y_top + 0.0005, label = "*"),
                    inherit.aes = FALSE, size = 3) +
          theme_classic() +
          theme(panel.grid = element_blank())
        
        ### 2D-a - Curve Merge IgG-A-M by IR ----------------------------------------
        
        ir.combined_predictions_strat_timepoints2 %>%
          left_join(cutoffs) %>%
          mutate(
            norm.estimate = norm.estimate - cutoff,
            norm.lower    = norm.lower    - cutoff,
            norm.upper    = norm.upper    - cutoff,
            cutoff        = 0
          ) %>%
          filter(Isotype %in% c("IgG","IgA","IgM"), Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          ggplot(aes(x = DPSO, y = norm.estimate)) +
          facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
          geom_line(aes(linetype = IR, color = Isotype), size = 1) +
          scale_color_manual(values = isotype_palette2) +
          scale_fill_manual(values = isotype_palette) +
          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
          scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
          geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
          theme_classic() +
          theme(axis.text.x = element_text(hjust = 1))
        
        ### 2D-b - Curve merge IgG1-4 by IR -----------------------------------------
        
        ir.combined_predictions_strat_timepoints2 %>%
          left_join(cutoffs2) %>%
          mutate(
            norm.estimate = norm.estimate - cutoff,
            norm.lower    = norm.lower    - cutoff,
            norm.upper    = norm.upper    - cutoff,
            cutoff        = 0
          ) %>%
          filter(Isotype %in% c("IgG1","IgG2","IgG3","IgG4"), Virus.resp == "HOMOTYPIC", Antigen == "NS1") %>%
          ggplot(aes(x = DPSO, y = norm.estimate)) +
          facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
          geom_line(aes(linetype = IR, color = Isotype), size = 1) +
          scale_color_manual(values = isotype_palette2) +
          scale_fill_manual(values = isotype_palette) +
          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
          scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
          geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
          theme_classic() +
          theme(axis.text.x  = element_text(hjust = 1))
                
                ### 2E - Seropositivity at 18M -----------------------------------------------
                
                aim1.dataset2 %>%
                  mutate(netMFI.Naive2 = MFI - cutoff) %>%
                  filter(Virus.resp == "HOMOTYPIC", Timepoint == "18", Antigen.bead == "NS1") %>%
                  mutate(
                    seropositivity = case_when(netMFI.Naive2 > 0 ~ "Seropositive",
                                               netMFI.Naive2 <= 0 ~ "Seronegative"),
                    `Secondary Ab` = factor(`Secondary Ab`, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))
                  ) %>%
                  count(seropositivity, IR, `Secondary Ab`, Antigen.bead, Timepoint, Virus.resp) %>%
                  ungroup() %>%
                  ggplot(aes(y = n, x = Antigen.bead, alpha = seropositivity, fill = `Secondary Ab`)) +
                  geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.3, color = "grey60") +
                  facet_grid(vars(IR), vars(`Secondary Ab`)) +
                  scale_fill_manual(values = isotype_palette) +
                  scale_alpha_manual(values = c("Seropositive" = 0.75, "Seronegative" = 0.15)) +
                  labs(y = "Seroprevalence at\n18 months post infection") +
                  geom_col(position = "fill") +
                  theme_classic() +
                  theme(axis.text.x  = element_text(angle = 90),
                        panel.spacing.x = unit(0.2, "lines"),
                        panel.spacing.y = unit(0.2, "lines"))
                
                ## Fig 3 - Cross-reactive E --------------------------------------------------
                
                ### 3A - Curves by isotype ---------------------------------------------------
                
                ir.combined_predictions_strat_timepoints2 %>%
                  left_join(cutoffs2) %>%
                  mutate(
                    norm.estimate = norm.estimate - cutoff,
                    norm.lower    = norm.lower    - cutoff,
                    norm.upper    = norm.upper    - cutoff,
                    cutoff        = 0
                  ) %>%
                  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                  filter(Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                         Virus.resp == "DV4_XR", Antigen == "E") %>%
                  ggplot(aes(x = DPSO, y = norm.estimate)) +
                  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
                  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = IR), alpha = 0.6) +
                  geom_line(aes(linetype = IR), size = 1) +
                  scale_fill_manual(values = ir_palette_DVXR) +
                  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
                  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
                  theme_classic() +
                  theme(axis.text.x = element_text(hjust = 1))
                
                ### 3B - Half Life Cv-3M -----------------------------------------------------
                
                bracket_data <- ir.combined_predictions_strat_rates.summary %>%
                  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                  filter(timepoint == "Cv-3M", Virus.resp == "DV4_XR", Antigen == "E") %>%
                  group_by(Isotype) %>%
                  filter(any(overlap_with_ref == FALSE)) %>%
                  summarize(
                    xmin  = 1 - 0.2,
                    xmax  = 1 + 0.2,
                    y_top = max(norm.upper[IR %in% c("P","S")]) + 0.01,
                    .groups = "drop"
                  )
                
                ir.combined_predictions_strat_rates.summary %>%
                  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                  filter(timepoint == "Cv-3M",
                         Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                         Virus.resp == "DV4_XR", Antigen == "E") %>%
                  ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                             color = IR, alpha = slope.sig, linetype = IR)) +
                  geom_linerange(size = 0.3, position = position_dodge(1)) +
                  geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
                  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
                  scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
                  scale_x_discrete(expand = c(0.4, 0.5)) +
                  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                           fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
                  scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                     guide = guide_legend(override.aes = list(fill = "black"))) +
                  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
                  scale_color_manual(values = ir_palette2_DVXR) +
                  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                  labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
                  geom_segment(data = bracket_data,
                               aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_segment(data = bracket_data,
                               aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_segment(data = bracket_data,
                               aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_text(data = bracket_data,
                            aes(x = (xmin + xmax)/2, y = y_top + 0.005, label = "*"),
                            inherit.aes = FALSE, size = 3) +
                  theme_classic() +
                  theme(panel.grid = element_blank())
                
                ### 3C - Half Life 6-18M -----------------------------------------------------
                
                bracket_data <- ir.combined_predictions_strat_rates.summary %>%
                  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                  filter(timepoint == "6M-18M", Virus.resp == "DV4_XR", Antigen == "E") %>%
                  group_by(Isotype) %>%
                  filter(any(overlap_with_ref == FALSE)) %>%
                  summarize(
                    xmin  = 1 - 0.2,
                    xmax  = 1 + 0.2,
                    y_top = max(norm.upper[IR %in% c("P","S")]) + 0.0001,
                    .groups = "drop"
                  )
                
                ir.combined_predictions_strat_rates.summary %>%
                  mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                  filter(timepoint == "6M-18M",
                         Isotype %in% c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
                         Virus.resp == "DV4_XR", Antigen == "E") %>%
                  ggplot(aes(x = timepoint, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                             color = IR, alpha = slope.sig, linetype = IR)) +
                  geom_linerange(size = 0.3, position = position_dodge(1)) +
                  geom_point(size = 1, stroke = 0, aes(fill = IR), position = position_dodge(1)) +
                  facet_wrap(vars(Isotype), ncol = 9, scales = "free_x") +
                  scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.2))) +
                  scale_x_discrete(expand = c(0.4, 0.5)) +
                  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                           fill = "grey", alpha = 0.2, inherit.aes = FALSE) +
                  scale_alpha_manual(values = c("Not Sig." = 0.3, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                     guide = guide_legend(override.aes = list(fill = "black"))) +
                  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = NA)) +
                  scale_color_manual(values = ir_palette2_DVXR) +
                  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                  labs(y = expression(t[1/2] ~ "(years)"), alpha = "Dynamics") +
                  geom_segment(data = bracket_data,
                               aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_segment(data = bracket_data,
                               aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.00005),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_segment(data = bracket_data,
                               aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.00005),
                               inherit.aes = FALSE, size = 0.2) +
                  geom_text(data = bracket_data,
                            aes(x = (xmin + xmax)/2, y = y_top + 0.00005, label = "*"),
                            inherit.aes = FALSE, size = 3) +
                  theme_classic() +
                  theme(panel.grid = element_blank())
                
                ### 3D-a - Curve Merge IgG-A-M by IR ----------------------------------------
                
                ir.combined_predictions_strat_timepoints2 %>%
                  left_join(cutoffs2) %>%
                  mutate(
                    norm.estimate = norm.estimate - cutoff,
                    norm.lower    = norm.lower    - cutoff,
                    norm.upper    = norm.upper    - cutoff,
                    cutoff        = 0
                  ) %>%
                  filter(Isotype %in% c("IgG","IgA","IgM"), Virus.resp == "DV4_XR", Antigen == "E") %>%
                  ggplot(aes(x = DPSO, y = norm.estimate)) +
                  facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
                  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
                  geom_line(aes(linetype = IR, color = Isotype), size = 1) +
                  scale_color_manual(values = isotype_palette2) +
                  scale_fill_manual(values = isotype_palette) +
                  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
                  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
                  theme_classic() +
                  theme(axis.text.x = element_text(hjust = 1))
                
                ### 3D-b - Curve merge IgG1-4 by IR -----------------------------------------
                
                ir.combined_predictions_strat_timepoints2 %>%
                  left_join(cutoffs2) %>%
                  mutate(
                    norm.estimate = norm.estimate - cutoff,
                    norm.lower    = norm.lower    - cutoff,
                    norm.upper    = norm.upper    - cutoff,
                    cutoff        = 0
                  ) %>%
                  filter(Isotype %in% c("IgG1","IgG2","IgG3","IgG4"), Virus.resp == "DV4_XR", Antigen == "E") %>%
                  ggplot(aes(x = DPSO, y = norm.estimate)) +
                  facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
                  geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Isotype), alpha = 0.6) +
                  geom_line(aes(linetype = IR, color = Isotype), size = 1) +
                  scale_color_manual(values = isotype_palette2) +
                  scale_fill_manual(values = isotype_palette) +
                  scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                  labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                  scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                  scale_y_continuous(limits = c(-1.5,8.5), breaks = c(0,2,4,6,8)) +
                  geom_rect(aes(xmin = 0, xmax = 540, ymin = -Inf, ymax = cutoff), fill = "grey", alpha = 0.02) +
                  theme_classic() +
                  theme(axis.text.x  = element_text(hjust = 1))
                        
                        ### 3E - Seropositivity at 18M -----------------------------------------------
                        
                        aim1.dataset2 %>%
                          mutate(netMFI.Naive2 = MFI - cutoff) %>%
                          filter(Virus.resp == "DV4_XR", Timepoint == "18", !Antigen.bead == "NS1") %>%
                          mutate(
                            seropositivity = case_when(netMFI.Naive2 > 0 ~ "Seropositive",
                                                       netMFI.Naive2 <= 0 ~ "Seronegative"),
                            `Secondary Ab` = factor(`Secondary Ab`, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))
                          ) %>%
                          count(seropositivity, IR, `Secondary Ab`, Antigen.bead, Timepoint, Virus.resp) %>%
                          ungroup() %>%
                          ggplot(aes(y = n, x = Antigen.bead, alpha = seropositivity, fill = `Secondary Ab`)) +
                          geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.3, color = "grey60") +
                          facet_grid(vars(IR), vars(`Secondary Ab`)) +
                          scale_fill_manual(values = isotype_palette) +
                          scale_alpha_manual(values = c("Seropositive" = 0.75, "Seronegative" = 0.15)) +
                          labs(y = "Seroprevalence at\n18 months post infection") +
                          geom_col(position = "fill") +
                          theme_classic() +
                          theme(axis.text.x  = element_text(angle = 90),
                                panel.spacing.x = unit(0.2, "lines"),
                                panel.spacing.y = unit(0.2, "lines"))
                        
                        ## Fig 4 - E vs EDIII vs NS1 IgG ---------------------------------------------
                        
                        ### 4A - Curves Homotypic ----------------------------------------------------
                        
                        ant.combined_predictions_strat_timepoints2 %>%
                          left_join(cutoffs2) %>%
                          mutate(
                            norm.estimate = norm.estimate - cutoff,
                            norm.lower    = norm.lower    - cutoff,
                            norm.upper    = norm.upper    - cutoff,
                            cutoff        = 0
                          ) %>%
                          filter(Isotype %in% c("IgG"), Virus.resp == "HOMOTYPIC") %>%
                          ggplot(aes(x = DPSO, y = norm.estimate)) +
                          facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
                          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Antigen), alpha = 0.6) +
                          geom_line(aes(linetype = IR, color = Antigen), size = 1) +
                          scale_color_manual(values = antigen_palette2) +
                          scale_fill_manual(values = antigen_palette) +
                          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 4B - Half Life Cv-3M -----------------------------------------------------
                        
                        df <- ant.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", Virus.resp == "HOMOTYPIC", timepoint == "Cv-3M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          group_by(IR) %>%
                          filter(overlap_with_ref == FALSE) %>%
                          summarize(
                            xmin  = min(as.numeric(Antigen)) - 0.15,
                            xmax  = max(as.numeric(Antigen)) + 0.15,
                            y_top = max(norm.upper) + 0.01,
                            .groups = "drop"
                          )
                        
                        df %>%
                          ggplot(aes(x = Antigen, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Antigen, alpha = slope.sig, linetype = IR, group = Antigen)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = IR), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(IR), ncol = 2) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter()) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(limits = c("E","EDIII","NS1"), values = antigen_palette) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.005),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 4C - Half Life 6M-18M ----------------------------------------------------
                        
                        df <- ant.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", Virus.resp == "HOMOTYPIC", timepoint == "6M-18M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          group_by(IR) %>%
                          filter(overlap_with_ref == FALSE) %>%
                          summarize(
                            xmin  = min(as.numeric(Antigen)) - 0.15,
                            xmax  = max(as.numeric(Antigen)) + 0.15,
                            y_top = max(norm.upper) + 0.002,
                            .groups = "drop"
                          )
                        
                        df %>%
                          ggplot(aes(x = Antigen, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Antigen, alpha = slope.sig, linetype = IR, group = Antigen)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = IR), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(IR), ncol = 2) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.1))) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(limits = c("E","EDIII","NS1"), values = antigen_palette) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.0004),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.0004),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.0004),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 4D - Curves DV XR --------------------------------------------------------
                        
                        ant.combined_predictions_strat_timepoints2 %>%
                          left_join(cutoffs2) %>%
                          mutate(
                            norm.estimate = norm.estimate - cutoff,
                            norm.lower    = norm.lower    - cutoff,
                            norm.upper    = norm.upper    - cutoff,
                            cutoff        = 0
                          ) %>%
                          filter(Isotype %in% c("IgG"), Virus.resp == "DV4_XR") %>%
                          ggplot(aes(x = DPSO, y = norm.estimate)) +
                          facet_wrap(vars(IR), ncol = 9, scales = "free_x") +
                          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Antigen), alpha = 0.6) +
                          geom_line(aes(linetype = IR, color = Antigen), size = 1) +
                          scale_color_manual(values = antigen_palette2) +
                          scale_fill_manual(values = antigen_palette) +
                          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 4E - DV_XR Cv-3M ---------------------------------------------------------
                        
                        df <- ant.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", Virus.resp == "DV4_XR", timepoint == "Cv-3M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          group_by(IR) %>%
                          filter(overlap_with_ref == FALSE) %>%
                          summarize(
                            xmin  = min(as.numeric(Antigen)) - 0.15,
                            xmax  = max(as.numeric(Antigen)) + 0.15,
                            y_top = max(norm.upper) + 0.01,
                            .groups = "drop"
                          )
                        
                        df %>%
                          ggplot(aes(x = Antigen, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Antigen, alpha = slope.sig, linetype = IR, group = Antigen)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = IR), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(IR), ncol = 2) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.2))) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(limits = c("E","EDIII","NS1"), values = antigen_palette) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.005),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 4F - DV XR 6M-18M --------------------------------------------------------
                        
                        df <- ant.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", Virus.resp == "DV4_XR", timepoint == "6M-18M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          group_by(IR) %>%
                          filter(overlap_with_ref == FALSE) %>%
                          summarize(
                            xmin  = min(as.numeric(Antigen)) - 0.15,
                            xmax  = max(as.numeric(Antigen)) + 0.15,
                            y_top = max(norm.upper) + 0.002,
                            .groups = "drop"
                          )
                        
                        df %>%
                          ggplot(aes(x = Antigen, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Antigen, alpha = slope.sig, linetype = IR, group = Antigen)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = IR), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(IR), ncol = 2) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0))) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(limits = c("E","EDIII","NS1"), values = antigen_palette) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.0004),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.0004),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.0004),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ## Fig 5 - Homotypic vs DV XR vs ZV XR --------------------------------------
                        
                        ### 5A - Curves Homotypic vs XR ----------------------------------------------
                        
                        vir.combined_predictions_strat_timepoints2 %>%
                          left_join(cutoffs2) %>%
                          mutate(
                            norm.estimate = norm.estimate - cutoff,
                            norm.lower    = norm.lower    - cutoff,
                            norm.upper    = norm.upper    - cutoff,
                            cutoff        = 0
                          ) %>%
                          filter(Isotype %in% c("IgG"), IR == "P") %>%
                          ggplot(aes(x = DPSO, y = norm.estimate)) +
                          facet_wrap(vars(Antigen), ncol = 9, scales = "free_x") +
                          geom_ribbon(aes(ymin = norm.lower, ymax = norm.upper, fill = Virus.resp), alpha = 0.6) +
                          geom_line(aes(linetype = IR, color = Virus.resp), size = 1) +
                          scale_color_manual(values = virus_resp_palette2) +
                          scale_fill_manual(values = virus_resp_palette3) +
                          scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                          labs(x = "Days post symptoms onset", y = expression(log[2]~"(MFI) - BG")) +
                          scale_x_log10(breaks = c(0,30,90,540), expand = c(0,0)) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 5B - Half Life Cv-3M -----------------------------------------------------
                        
                        df <- vir.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", IR == "P", timepoint == "Cv-3M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          filter(overlap_with_ref == FALSE) %>%
                          mutate(
                            x_numeric = as.numeric(factor(Virus.resp, levels = levels(df$Virus.resp))),
                            xmin  = x_numeric - 0.15,
                            xmax  = x_numeric + 0.15,
                            y_top = norm.upper + 0.01
                          )
                        
                        df %>%
                          ggplot(aes(x = Virus.resp, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Virus.resp, alpha = slope.sig, linetype = IR, group = Virus.resp)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = Virus.resp), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(Antigen), ncol = 3) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.5))) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(values = virus_resp_palette3) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.005),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic() +
                          theme(axis.text.x = element_text(hjust = 1))
                        
                        ### 5C - Half Life 6M-18M ----------------------------------------------------
                        
                        df <- vir.combined_predictions_strat_rates.summary %>%
                          filter(Isotype == "IgG", IR == "P", timepoint == "6M-18M") %>%
                          mutate(Antigen = factor(Antigen, levels = c("E","EDIII","NS1")))
                        
                        bracket_data_ant <- df %>%
                          filter(overlap_with_ref == FALSE) %>%
                          mutate(
                            x_numeric = as.numeric(factor(Virus.resp, levels = levels(df$Virus.resp))),
                            xmin  = x_numeric - 0.15,
                            xmax  = x_numeric + 0.15,
                            y_top = norm.upper + 0.001
                          )
                        
                        df %>%
                          ggplot(aes(x = Virus.resp, y = norm.estimate, ymin = norm.lower, ymax = norm.upper,
                                     color = Virus.resp, alpha = slope.sig, linetype = IR, group = Virus.resp)) +
                          geom_linerange(size = 0.3, position = position_dodge(0.5)) +
                          geom_point(aes(fill = Virus.resp), size = 1, stroke = 1, position = position_dodge(0.5)) +
                          scale_x_discrete(expand = c(0.4, 0)) +
                          facet_wrap(vars(Antigen), ncol = 3) +
                          annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey", alpha = 0.2) +
                          scale_y_continuous(labels = formatter(), expand = expansion(mult = c(0, 0.5))) +
                          scale_linetype_manual(values = c(P = "dashed", S = "solid")) +
                          scale_alpha_manual(values = c("Not Sig." = 0.2, "Sig. Rising" = 1, "Sig. Waning" = 1),
                                             guide = guide_legend(override.aes = list(fill = "black"))) +
                          scale_shape_manual(values = c("Ref." = 1, "FALSE" = 17, "TRUE" = 19)) +
                          scale_color_manual(values = virus_resp_palette3) +
                          labs(y = expression(t[1/2] ~ "(years)"), x = NULL, shape = "CI overlap", alpha = "Dynamics") +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmin, xend = xmin, y = y_top, yend = y_top - 0.0005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_segment(data = bracket_data_ant,
                                       aes(x = xmax, xend = xmax, y = y_top, yend = y_top - 0.0005),
                                       inherit.aes = FALSE, size = 0.2) +
                          geom_text(data = bracket_data_ant,
                                    aes(x = (xmin + xmax)/2, y = y_top + 0.0005),
                                    label = "*", inherit.aes = FALSE, size = 3) +
                          theme_classic()
                        
                        ## Fig 7 - Interindividual variation -----------------------------------------
                        
                        formula <- y ~ x
                        
                        aim1.dataset2.hl_2 <- aim1.dataset2.hl %>%
                          group_by(Code, Virus.resp, IR, Isotype, Antigen) %>%
                          mutate(
                            fchange_6_18_cat = case_when(
                              fchange_6_18 < 0 ~ "wane",
                              fchange_6_18 > 0 ~ "rise"
                            ),
                            fchange_6_18_cat = factor(fchange_6_18_cat, levels = c("rise","wane"))
                          )
                        
                        ### 7A - Violin plot Magnitude ------------------------------------------------
                        
                        aim1.dataset2.hl_2 %>%
                          filter(Antigen == "E", Virus.resp == "DV4_XR", !is.na(netMFI.Naive_18)) %>%
                          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                          ggplot(aes(x = IR, y = netMFI.Naive_18)) +
                          geom_violin(fill = "grey75", color = NA, alpha = 0.5, width = 1,
                                      draw_quantiles = c(0.25, 0.5, 0.75)) +
                          geom_jitter(aes(color = fchange_6_18_cat), width = 0.2, size = 0.2, alpha = 0.8) +
                          facet_grid(cols = vars(Isotype), scales = "fixed") +
                          scale_color_manual(values = hl_palette) +
                          labs(x = "Infection group",
                               y = "(log2) MFI - BG\n18 Months post symptoms onset",
                               color = "Category") +
                          theme_classic() +
                          theme(panel.spacing.x = unit(0.2, "lines"),
                                panel.spacing.y = unit(0.2, "lines"))
                        
                        ### 7B - Violin plot Half-life ------------------------------------------------
                        
                        aim1.dataset2.hl_2 %>%
                          filter(Antigen == "E", Virus.resp == "DV4_XR", !is.na(fchange_6_18_cat)) %>%
                          mutate(Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))) %>%
                          ggplot(aes(x = IR, y = fchange_6_18)) +
                          geom_violin(fill = "grey75", color = NA, alpha = 0.5, width = 1,
                                      draw_quantiles = c(0.25, 0.5, 0.75)) +
                          geom_jitter(aes(color = fchange_6_18_cat), width = 0.15, size = 0.2, alpha = 0.8) +
                          geom_hline(yintercept = 0, color = "grey60", size = 0.2) +
                          scale_y_continuous(labels = formatter.years(), expand = expansion(mult = c(0, 0.1))) +
                          facet_grid(cols = vars(Isotype), scales = "fixed") +
                          scale_color_manual(values = hl_palette) +
                          labs(x = "Infection group", y = "Half-life (years)", color = "Category") +
                          theme_classic() +
                          theme(panel.spacing.x = unit(0.2, "lines"),
                                panel.spacing.y = unit(0.2, "lines"))
                        
                        ### 7C - Barplot % rise vs wane ----------------------------------------------
                        
                        aim1.dataset2.hl_2 %>%
                          filter(Antigen %in% c("E"), Virus.resp %in% c("DV4_XR"), !is.na(fchange_6_18_cat)) %>%
                          mutate(
                            fchange_6_18_cat = factor(fchange_6_18_cat, levels = c("rise","stable","wane")),
                            Isotype = factor(Isotype, levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"))
                          ) %>%
                          ggplot(aes(x = IR, fill = fchange_6_18_cat)) +
                          geom_bar(position = "fill", colour = NA, width = 0.9, alpha = 0.8) +
                          facet_grid(cols = vars(Isotype), scales = "fixed") +
                          scale_fill_manual(values = hl_palette) +
                          scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                             expand = expansion(mult = c(0, 0.02))) +
                          labs(x = "Infection group", y = "Percent of Observations") +
                          theme_classic() +
                          theme(panel.spacing.x = unit(0.2, "lines"),
                                panel.spacing.y = unit(0.2, "lines"))
                        
                        ## Bayesian curves ============================================================
                        
                        y_preset_E   <- list(y_limits = c(-1.5, 8.5), y_breaks = c(0,2,4,6,8))
                        y_preset_NS1 <- list(y_limits = c(-1.5, 8.5), y_breaks = c(0,2,4,6,8))
                        
                        pick_y_preset <- function(antigen) {
                          if (antigen == "E")   return(y_preset_E)
                          if (antigen == "NS1") return(y_preset_NS1)
                          list()
                        }
                        
                        pick_ir_palette <- function(virus_resp) {
                          if (virus_resp == "HOMOTYPIC") return(ir_palette)
                          if (virus_resp == "DV4_XR")   return(ir_palette_DVXR)
                          stop("No palette rule for Virus.resp = ", virus_resp)
                        }
                        

                        plot_pred_traj <- function(
    df,
    secondary_ab = c("IgM","IgA","IgG","IgG1","IgG2","IgG3","IgG4"),
    facet_col = "SecondaryAb",
    x_col     = "DPSO",
    y_col     = "norm.estimate",
    lo_col    = "norm.lower",
    hi_col    = "norm.upper",
    ir_col    = "IR",
    ir_palette,
    facet_ncol = 7,
    x_breaks   = c(30, 90, 540),
    y_limits   = c(-1.5, 8.5),
    y_breaks   = c(0,2,4,6,8),
    top_label  = NULL
                        ) {
                          d <- df %>%
                            filter(.data[[facet_col]] %in% secondary_ab) %>%
                            mutate(
                              "{facet_col}" := factor(.data[[facet_col]], levels = secondary_ab),
                              .x = suppressWarnings(as.numeric(.data[[x_col]]))
                            ) %>%
                            filter(is.finite(.x), .x > 0)
                            
                          
                          x_min <- min(d$.x, na.rm = TRUE)
                          x_max <- max(d$.x, na.rm = TRUE)
                          
                          ggplot(d, aes(x = .x, y = .data[[y_col]])) +
                            facet_wrap(vars(.data[[facet_col]]), ncol = facet_ncol, scales = "free_x") +
                            geom_rect(
                              data = data.frame(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
                              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                              inherit.aes = FALSE,
                              fill = "grey", alpha = 0.2
                            ) +
                            # geom_ribbon(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]],
                            #                 fill = .data[[ir_col]]), alpha = 0.6) +
                            geom_line(aes(group = .data[[ir_col]], linetype = .data[[ir_col]],color = .data[[ir_col]]),
                                      linewidth = 1) +
                            scale_color_manual(values = ir_palette) +
                            scale_linetype_manual(values = c("P" = "dashed", "S" = "solid")) +
                            labs(x = "Days post symptoms onset",
                                 y = expression(log[2]~"(MFI) - BG"),
                                 title = top_label) +
                            scale_x_log10(breaks = x_breaks, expand = c(0, 0)) +
                            scale_y_continuous(limits = y_limits, breaks = y_breaks) +
                            theme_classic() +
                            theme(axis.text.x = element_text(hjust = 1))
                          }
                        
                        make_pred_traj_plot <- function(data, virus_resp, antigen,
                                                        w = 15, h = 4.5, units = "cm",
                                                        cutoff = NULL) {
                          df <- data %>%
                            filter(Virus.resp == virus_resp, Antigen.bead == antigen)
                          
                          if (nrow(df) == 0) {
                            message("SKIP (no data): ", virus_resp, " / ", antigen)
                            return(NULL)
                          }
                          
                          pal <- pick_ir_palette(virus_resp)
                          yp  <- pick_y_preset(antigen)
                          
                          do.call(
                            plot_pred_traj,
                            c(list(df = df, ir_palette = pal,
                                   top_label = paste0(antigen, "  ", virus_resp)), yp)
                          )
                        }
                        
                        Virus.resp.values <- c("HOMOTYPIC", "DV4_XR")
                        antigens          <- c("E", "NS1", "EDIII")
                        
                        pred_plots <- list()
                        for (v in Virus.resp.values) {
                          for (a in antigens) {
                            key <- paste(v, a, sep = "__")
                            pred_plots[[key]] <- make_pred_traj_plot(data = pred_traj_df,
                                                                     virus_resp = v, antigen = a)
                          }
                        }
                        
                        pred_plots[["HOMOTYPIC__E"]]
                        pred_plots[["HOMOTYPIC__NS1"]]
                        pred_plots[["DV4_XR__E"]]
                        pred_plots[["DV4_XR__NS1"]]
                        
                        ## Bayesian slopes half-life visualization ====================================
                        
                        bracket_preset_Cv3M <- list(
                          bracket_halfwidth = 0.2,
                          bracket_y_pad     = 0.01,
                          tick_drop         = 0.005,
                          star_rise         = 0.005
                        )
                        
                        bracket_preset_6M18M <- list(
                          bracket_halfwidth = 0.25,
                          bracket_y_pad     = 0.001,
                          tick_drop         = 0.0005,
                          star_rise         = 0.0005
                        )
                        
                        pick_bracket <- function(dpso) {
                          if (dpso == "Cv-3M")  return(bracket_preset_Cv3M)
                          if (dpso == "6M-18M") return(bracket_preset_6M18M)
                          list()
                        }
                        
                        ylims_by_DPSO <- slope_comparisons %>%
                          filter(Virus.resp %in% c("HOMOTYPIC","DV4_XR"),
                                 Antigen.bead %in% c("E","NS1")) %>%
                          group_by(DPSO) %>%
                          summarise(ymin = min(q2.5, na.rm = TRUE),
                                    ymax_ci = max(q97.5, na.rm = TRUE),
                                    .groups = "drop") %>%
                          rowwise() %>%
                          mutate(
                            preset = list(pick_bracket(DPSO)),
                            pad    = preset$bracket_y_pad %||% 0.001,
                            rise   = preset$star_rise     %||% 0.0005,
                            ymax   = ymax_ci + pad + rise
                          ) %>%
                          select(DPSO, ymin, ymax)
                        
                        pick_palette <- function(virus_resp) {
                          if (virus_resp == "HOMOTYPIC") return(ir_palette2)
                          if (virus_resp == "DV4_XR")   return(ir_palette2_DVXR)
                          stop("No palette rule for Virus.resp = ", virus_resp)
                        }
                        
                        plot_slope <- function(
    df,
    ir_palette,
    x_levels = c("IgG","IgA","IgM","IgG1","IgG2","IgG3","IgG4"),
    bracket_flag_col  = "IR_overlap_with_ref",
    bracket_flag_val  = "Sig. Diff",
    x_col   = "SecondaryAb",
    y_col   = "median",
    lo_col  = "q2.5",
    hi_col  = "q97.5",
    ir_col  = "IR",
    sig_col = "slope.sig",
    dodge_width       = 0.6,
    bracket_halfwidth = 0.2,
    bracket_y_pad     = 0.001,
    tick_drop         = 0.0005,
    star_rise         = 0.0005,
    x_lab        = "Antibody isotype and subtype",
    alpha_levels = c("Not Sig." = 0.3, "Sig. Waning" = 1, "Sig. Rising" = 1)
                        ) {
                          df <- df %>%
                            mutate(
                              "{x_col}" := factor(.data[[x_col]], levels = x_levels),
                              facet_label = paste0(.data[["Antigen.bead"]], "  ", .data[["DPSO"]])
                            )
                          
                          bracket_data <- df %>%
                            group_by(.data[[x_col]]) %>%
                            filter(any(.data[[bracket_flag_col]] == bracket_flag_val, na.rm = TRUE)) %>%
                            summarize(
                              x_num = as.numeric(first(.data[[x_col]])),
                              xmin  = x_num - bracket_halfwidth,
                              xmax  = x_num + bracket_halfwidth,
                              y_top = max(.data[[hi_col]], na.rm = TRUE) + bracket_y_pad,
                              .groups = "drop"
                            )
                          
                          ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]],
                                         color = .data[[ir_col]], alpha = factor(.data[[sig_col]]))) +
                            annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
                                     fill = "grey", alpha = 0.2) +
                            geom_point(size = 2.1, position = position_dodge(width = dodge_width)) +
                            geom_errorbar(aes(ymin = .data[[lo_col]], ymax = .data[[hi_col]]),
                                          width = 0.12, linewidth = 0.35,
                                          position = position_dodge(width = dodge_width)) +
                            geom_segment(data = bracket_data,
                                         aes(x = xmin, xend = xmax, y = y_top, yend = y_top),
                                         inherit.aes = FALSE, linewidth = 0.2) +
                            geom_segment(data = bracket_data,
                                         aes(x = xmin, xend = xmin, y = y_top, yend = y_top - tick_drop),
                                         inherit.aes = FALSE, linewidth = 0.2) +
                            geom_segment(data = bracket_data,
                                         aes(x = xmax, xend = xmax, y = y_top, yend = y_top - tick_drop),
                                         inherit.aes = FALSE, linewidth = 0.2) +
                            geom_text(data = bracket_data,
                                      aes(x = (xmin + xmax)/2, y = y_top + star_rise, label = "*"),
                                      inherit.aes = FALSE, size = 3) +
                            labs(x = x_lab, y = expression(t[1/2] ~ "(years)"), title = NULL) +
                            scale_color_manual(values = ir_palette) +
                            scale_alpha_manual(values = alpha_levels, guide = "none") +
                            facet_wrap(vars(facet_label)) +
                            theme_classic() +
                            theme(
                              strip.background = element_blank(),
                              panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.2),
                              plot.title       = element_blank()
                            )
                        }
                        
                        make_slope_plot <- function(data, virus_resp, antigen, dpso,
                                                    w = 7.5, h = 4.5, units = "cm") {
                          df <- data %>%
                            filter(Virus.resp == virus_resp, Antigen.bead == antigen, DPSO == dpso)
                          
                          if (nrow(df) == 0) {
                            message("SKIP (no data): ", virus_resp, " / ", antigen, " / ", dpso)
                            return(NULL)
                          }
                          
                          pal <- pick_palette(virus_resp)
                          br  <- pick_bracket(dpso)
                          yl  <- ylims_by_DPSO %>% filter(DPSO == dpso)
                          
                          do.call(plot_slope, c(list(df = df, ir_palette = pal), br)) +
                            coord_cartesian(ylim = c(yl$ymin, yl$ymax))
                        }
                        
                        Virus.resp.values <- c("HOMOTYPIC", "DV4_XR")
                        antigens          <- c("E", "NS1")
                        dpsos             <- sort(unique(slope_comparisons$DPSO))
                        
                        plots <- list()
                        for (v in Virus.resp.values) {
                          for (a in antigens) {
                            for (t in dpsos) {
                              key         <- paste(v, a, t, sep = "__")
                              plots[[key]] <- make_slope_plot(data = slope_comparisons,
                                                              virus_resp = v, antigen = a, dpso = t)
                            }
                          }
                        }
                        
                        plots
                        
                        plots[["HOMOTYPIC__E__6M-18M"]]
                        plots[["HOMOTYPIC__E__Cv-3M"]]
                        plots[["HOMOTYPIC__NS1__6M-18M"]]
                        plots[["HOMOTYPIC__NS1__Cv-3M"]]
                        plots[["DV4_XR__E__6M-18M"]]
                        plots[["DV4_XR__E__Cv-3M"]]
                        plots[["DV4_XR__NS1__6M-18M"]]
                        plots[["DV4_XR__NS1__Cv-3M"]]
                        
                        
                        
   
                        
          