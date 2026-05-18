#SARS-CoV-2 Correlates Analysis with Columbia neut data
#Programmer: Jose Victor Zambrana
#Date created: 05/21/2024
#Date modified: 06/26/2025


# Libraries ---------------------------------------------------------------------------------------------------------------------------------------------------------------

library(haven) # read SAS
library(mgcv) # GAM modeling
library(ggplot2) # plotting
library(gridExtra) # combine plots
library(ggpubr) #  plotting
library(tidyverse) # data manipulation
library(rstatix) # univariate statistics with nice output
library(readxl) # read excel
library(dplyr) # data manipulation
library(lubridate) # date manipulation
library(broom) # tidy model output
library(mediation) # for mediation analyiss
library(gtsummary) # gtsummary
library(GGally) # for correlation plot
library(geepack) # for gee analysis


# data reading ---------------------------------------------------------------------------------------------------------------------------------------------------------------

# dataset containing the neuts
datasero <-read_sas("Data/phi_covcop_ho_neut_sero_24sep24.sas7bdat")

# dataset containing the spike
spike.titers <- read_sas("Data/covcop_spiketiters_24sep24.sas7bdat")

# dataset containing infection history
inf.hist.curated <- read_csv("Data/JAN17-2024 Infection history curation.csv")

# dataset containing aggregated weekly infection counts in HICS and sequencing
aggregated.seq = readxl::read_xlsx("Data/Abby/HICS Cases Sequencing Over Time 14Jan25.xlsx", sheet = 2) 

# dataset containing aggregated weekly infection counts in Managua
aggregated.seq.minsa = readxl::read_xlsx("Data/freq_sars_cov_2_linajes.xlsx", sheet = 2) 


# data manipulation ---------------------------------------------------------------------------------------------------------------------------------------------------------------



# pre-processing on spike titers
pre.titers2 = spike.titers %>%
  dplyr::select(codigo, sdt = sampledate, titer_num, assay_dt) %>%
  arrange(codigo,  sdt, assay_dt) %>% # sort by codigo, sampledate and assaydate
  right_join(datasero %>% dplyr::select(codigo, sampledate), relationship = "many-to-many") %>% # many-to-many cause a couple of indviduals were in both waves
  filter(sdt <= sampledate) %>% # filter those from activation or before
  group_by(codigo, sampledate) %>%
  filter(sdt  == max(sdt)) %>%  # select the most immediate sample prior or during activation
  summarise(spike_titer.jv  = median(titer_num ), sampledate_spike = first(sdt))


# manipulating main dataset containing neuts
datasero2 = datasero %>%
  
  # Step 1: Add a 'variant' column based on 'sequence2' (imputed) 
  mutate(
    variant = factor(
      case_when(
        grepl(pattern = "BA\\.1", x = sequence2) ~ "BA.1 wave",  # Assign "BA.1 wave" if 'sequence2' contains "BA.1"
        grepl(pattern = "BA\\.2", x = sequence2) ~ "BA.2 wave"   # Assign "BA.2 wave" if 'sequence2' contains "BA.2"
      )
    )
  ) %>%
  
  # Step 2: Filter out rows where 'variant' is NA, excluding participants with uncertain subvariant information
  filter(!is.na(variant)) %>% 
  
  # Step 3: Merge with the curated infection history dataset
  left_join(inf.hist.curated, by = join_by(codigo, sampledate)) %>%
  
  # Step 4: Create and transform multiple variables
  mutate(
    
    inf.hist.i = case_when(is.na(inf.hist.i) ~ "Uninfected",
                           T ~ inf.hist.i),
    
    inf.hist.v = case_when(is.na(inf.hist.v) ~ "Unvaccinated",
                           T ~ inf.hist.v),
    
    # Define vaccination status ('any.vax') based on 'inf.hist.v'
    any.vax = case_when(
      inf.hist.v == "Unvaccinated" ~ "Unvaccinated",                               # Label as "Unvaccinated"
      inf.hist.v %in% c("Vax", "2Vax", "+2Vax") ~ "Vax"                            # Label as "Vax" for various vaccination statuses
    ),
    any.vax = factor(any.vax, levels = c("Vax", "Unvaccinated")),                  # Convert 'any.vax' to a factor with specified levels
    
    # Define infection status ('any.inf') based on 'inf.hist.i'
    any.inf = case_when(
      inf.hist.i == "Uninfected" ~ "Uninfected",                                   # Label as "Uninfected"
      inf.hist.i %in% c("CoV", "2CoV", "+2CoV") ~ "Infected"                      # Label as "Infected" for various infection statuses
    ),
    any.inf = factor(any.inf, levels = c("Infected", "Uninfected")),                # Convert 'any.inf' to a factor with specified levels
    
    # Calculate the time difference between sample date and infection date
    t.sample.and.infection = as.integer(sampledate - inf_date),                    # Time between sample and infection dates
    
    # Define symptomatic infection ('sympt.inf') based on infection and severity
    sympt.inf = factor(
      case_when(
        pt_inf == 1 & severity >= 1 ~ 1,                                           # Label as 1 if infected and severity is at mild
        TRUE ~ 0                                                                    # Otherwise, label as 0
      )
    ),
    
    # Convert 'sexo' to a factor 
    sexo = as.factor(sexo),
    
    # Convert 'pt_inf' to a factor 
    pt_inf = as.factor(pt_inf)
  ) %>%
  
  # Step 5: Merge with pre-existing titers datasets
  left_join(pre.titers2, by = join_by(codigo, sampledate)) %>% 
  
  # Step 6:  calculate additional time differences
  mutate(
    
    # Calculate time difference between spike sample date and infection date
    t.sample.and.infection.spike = as.integer(sampledate_spike - inf_date)      # Time between spike sample and infection dates
    
  )



# -----------------------------------3
# Transforming datasero2 into Long Format
# -----------------------------------3
datasero2.long <- datasero2 %>%
  
  # Step 1: Select Relevant Columns
  select(
    codigo,          # Participant code or identifier
    hh_code,         # household code
    sexo,            # Sex of the participant
    any.vax,         # Vaccination status
    any.inf,         # Infection status
    age_ingreso,     # Age at admission or entry
    variant,         # Variant type (e.g., BA.1 wave, BA.2 wave)
    pt_inf,          # Participant infection flag
    sympt.inf,       # Symptomatic infection flag
    ID50_num_BA_1,   # Neutralization titer for BA.1 variant
    ID50_num_BA_2,   # Neutralization titer for BA.2 variant
    ID50_num_D614G,  # Neutralization titer for D614G variant
    spike_titer.jv,  # spike titer elisa
    time.since.inf,  # time since last infection
    time.since.vax   # time since vaccination
  ) %>%
  
  # Step 2: Pivot Data from Wide to Long Format
  pivot_longer(
    cols = c(
      ID50_num_BA_1,    # Neutralization titer for BA.1 variant
      ID50_num_BA_2,    # Neutralization titer for BA.2 variant
      ID50_num_D614G,   # Neutralization titer for D614G variant
      spike_titer.jv    # Spike protein binding titer
    ),
    names_to = "Assay",       # New column to store the assay type
    values_to = "Titers"      # New column to store the corresponding titer values
  ) %>%
  
  # Step 3: Apply Log Transformation to Titer Values
  mutate(
    Titers = log(Titers, base = 4)  # Log-transform the 'Titers' column using base 4
  ) %>%
  
  # Step 4: Relabel and Factorize the 'Assay' Column, and Factorize Other Categorical Variables
  mutate(
    # Relabel the 'Assay' names to more descriptive labels with Unicode subscripts
    Assay = case_when(
      Assay == "ID50_num_BA_1"  ~ "BA.1 Neutralization, Log\u2084[ID\u2085\u2080]",
      Assay == "ID50_num_BA_2"  ~ "BA.2 Neutralization, Log\u2084[ID\u2085\u2080]",
      Assay == "ID50_num_D614G" ~ "D614G Neutralization, Log\u2084[ID\u2085\u2080]",
      Assay == "spike_titer.jv" ~ "Spike Binding, Log\u2084[Titer]"
     
    ),
    
    # Convert the 'Assay' column to a factor with specified level order
    Assay = factor(
      Assay, 
      levels = c(
        "BA.1 Neutralization, Log\u2084[ID\u2085\u2080]",
        "BA.2 Neutralization, Log\u2084[ID\u2085\u2080]",
        "D614G Neutralization, Log\u2084[ID\u2085\u2080]",
        "Spike Binding, Log\u2084[Titer]"
      )
    ),
    
    # Re-factorize 'any.vax' to ensure consistent ordering of levels
    any.vax = factor(any.vax, levels = c("Unvaccinated", "Vax")),
    
    # Re-factorize 'any.inf' to ensure consistent ordering of levels
    any.inf = factor(any.inf, levels = c("Uninfected", "Infected"))
  )


# filter NA values
datasero2.long.noNA = datasero2.long %>%
  filter(!is.na(Titers)) # 1 in spike


# manipulating dataset for mediation analysis
datasero3 <- datasero2 %>%
  
   #Create 'neut.homotypic' Based on Variant and Apply Log Transformations
  mutate(
    # Assign homotypic neutralization titer based on the variant
    neut.homotypic = case_when(
      variant == "BA.1 wave" ~ ID50_num_BA_1,
      variant == "BA.2 wave" ~ ID50_num_BA_2
    ),
    
    # Log-transform 
    neut.homotypic = log(neut.homotypic, base = 4),
    neut.ancestral = log(ID50_num_D614G, base = 4),
    spike_titer.jv = log(spike_titer.jv, base = 4)
  )



# data analysis ---------------------------------------------------------------------------------------------------------------------------------------------------------------

# useful functions for pasting 95% CI
paste_ci <- function(estimate, lower, upper) {
  paste(round(estimate, 2), " (", round(lower, 2), ", ", round(upper, 2), ")", sep = "")
}

paste_ci_round <- function(estimate, lower, upper) {
  paste(round(estimate), " (", round(lower), ", ", round(upper), ")", sep = "")
}


###########################################################################3
# Table 1 -------------------------------------------------------------------
###########################################################################3

datasero2 %>%
  group_by(hh_code, variant) %>%
  summarise(count = n()) %>%
  summarise(count = n()) %>%
  arrange(desc(count))


# notes, I built the labels manually in Word
datasero2 %>%
  mutate(hh_seq = as.factor(hh_seq)) %>%
      tbl_summary(
         by = variant,
        include = c(pt_inf, 
                    severity,
                    age_ingreso,
                    sexo,
                    inf.hist.i,
                    time.since.inf,
                    inf.hist.v,
                    time.since.vax,
                    ID50_num_BA_1, 
                    ID50_num_BA_2, 
                    ID50_num_D614G,
                   spike_titer.jv,
                   t.sample.and.infection,
                   t.sample.and.infection.spike
                    ), missing = "no") |>
  add_overall()  %>%
  add_p() %>%
  as_gt() %>%
  gt::gtsave(filename = "Results/Table 1 - Participants characteristics.docx")


# appendix to the table
# get overall number of households and mean number of participants 
datasero2 %>%
  group_by(hh_code) %>%
  summarise(N = n()) %>%
  ungroup() %>%
  summarise(HouseholdN = n(), Participants25 = quantile(N,probs = c(0.25)),
            Participants50 = quantile(N,probs = c(0.50)),
            Participants75 = quantile(N,probs = c(0.75)),
            Participants = paste0(Participants50, " (",Participants25, ", ", Participants75, ")")
  ) %>%
  mutate(estimate = paste_ci_round(Participants50,Participants25,Participants75)) %>%
  select(estimate)


# get number of households and mean number of participants per omicron wave
datasero2 %>%
  group_by(hh_code, variant) %>%
  summarise(N = n()) %>%
  ungroup() %>%
  group_by(variant) %>%
  summarise(HouseholdN = n(), Participants25 = quantile(N,probs = c(0.25)),
            Participants50 = quantile(N,probs = c(0.50)),
            Participants75 = quantile(N,probs = c(0.75)),
            Participants = paste0(Participants50, " (",Participants25, ", ", Participants75, ")")
            ) %>%
  mutate(estimate = paste_ci_round(Participants50,Participants25,Participants75)) %>%
  select(variant, estimate)





###########################################################################3
# Table S1-2 -------------------------------------------------------------------
###########################################################################3


process.univ = function(x){
  x %>%mutate(`Median difference` = paste_ci(estimate, conf.low, conf.high)) %>%
    select(-estimate, -conf.low, -conf.high, -group1, -group2, - alternative, -method, -`.y.`) %>%
    relocate(`Median difference`, .after = Assay) %>%
    gt::gt(rownames_to_stub = TRUE)
  
}



datasero2.long %>%
  group_by(Assay, variant) %>%
  rstatix::wilcox_test(Titers ~ pt_inf, detailed = T) %>%
  rename(Infected = n1, Uninfected = n2) %>%
  process.univ() %>%
  gt::gtsave(filename = "Results/Table S1 - Univariate analysis by infection.docx")



datasero2.long %>%
  group_by(Assay, variant) %>%
  rstatix::wilcox_test(Titers ~ sympt.inf, detailed = T) %>%
  rename(`Symptomatic` = n1, `Not symptomatic` = n2) %>%
  process.univ() %>%
  gt::gtsave(filename =  "Results/Table S2 - Univariate analysis by symptomatic infections.docx")


  
###########################################################################3
# Figure S5 - Correlation  -------------------------------------------------------------------
###########################################################################3

# funny warnings, do not know how to skip it
ggpairs(datasero2 %>%
          select(`BA.1 Neutralization` = "ID50_num_BA_1",
                 `BA.2 Neutralization` = "ID50_num_BA_2",
                 `D614G Neutralization` = "ID50_num_D614G",
                 `Spike Binding` = spike_titer.jv,
                 variant) %>%
          mutate_if(is.numeric,.funs = function(x)(log(x, base = 4))) %>%
          na.omit(),                 # filtering 1 guy without spike titers
        columns = c(1:4),        # Columns
        aes(color = variant, fill = variant,  alpha = 0.8), # Color by group (cat. variable), # Transparency
        upper = list(continuous = wrap("smooth", method = "loess")),
        lower = list(continuous = wrap("cor", method = "pearson"))) +
  theme_classic2(base_size = 20) +
  scale_color_manual(values = c(`BA.1 wave` = "#9c954d", `BA.2 wave` = "#b067a3"))+
  scale_fill_manual(values = c(`BA.1 wave` = "#9c954d", `BA.2 wave` = "#b067a3")) 

ggsave("Results/Figure S5 - Correlation.svg", device = "svg", height = 11, width = 11)  


######################################################################################4
# Figure S6-7  ###################################################################
##################################################################################33

create_violin_plot <- function(data, grouping_var, x_label, output_file) {
  
  # Convert 0/1 to No/Yes in the grouping_var column
  data[[grouping_var]] <- factor(
    data[[grouping_var]],
    levels = c(0, 1),    # The original numeric values
    labels = c("No", "Yes")  # The new labels
  )
  
  # Perform the Wilcoxon tests
  wilcox_results <- data %>%
    group_by(Assay, variant) %>%
    rstatix::wilcox_test(as.formula(paste("Titers ~", grouping_var)), detailed = TRUE) %>%
    mutate(p_label = paste0("p = ", signif(p, digits = 3)))
  
  # Merge p-values back to the dataset for annotation
  data_with_p <- data %>%
    left_join(wilcox_results %>% select(Assay, variant, p_label), by = c("Assay", "variant")) 
  
  # Create the violin plot with facet_wrap
  plot <- ggplot(data_with_p, aes_string(x = grouping_var, y = "Titers")) +
    geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.7) +
    ggbeeswarm::geom_beeswarm() +
    facet_wrap(~ Assay + variant, ncol = 4) +  # Facet by Assay and variant
    theme_classic(base_size = 18) +
    labs(
      x = x_label,
      y = "Titers"
    ) +
    # Annotate p-values inside each facet
    geom_text(aes(label = p_label),
              x = 1.5, y = Inf, vjust = 1.5, size = 5)
  
  # Save the plot
  ggsave(output_file, plot, height = 10, width = 15)
}


create_violin_plot(
  data = datasero2.long.noNA,
  grouping_var = "pt_inf",                # <- Will be converted inside the function
  x_label = "Infected",
  output_file = "Results/Figure S6 Univariate_distributions_by_Infection.svg"
)


create_violin_plot(
  data = datasero2.long.noNA,
  grouping_var = "sympt.inf",             # <- Will be converted inside the function
  x_label = "Symptomatic infection",
  output_file = "Results/Figure S7 Univariate_distributions_by_Symptomatic_infection.svg"
)


###########################################################################3
# Table S3-4 -------------------------------------------------------------------
###########################################################################3

# Define the function
create_multivariate_table <- function(data, outcome_var, output_file) {
  table_results <- data %>%
    group_by(Assay, variant) %>%
    group_modify(
      ~ {
        mod <- glm(
          as.formula(paste(outcome_var, "~ Titers + any.vax + any.inf + age_ingreso")),
          data = .x,
          family = binomial(link = "logit")
        )
        
        mod_tidy <- broom::tidy(mod, exponentiate = TRUE, conf.int = TRUE)
        mod_glance <- broom::glance(mod)
        
        mod_tidy %>%
          mutate(
            df.residual = mod_glance$df.residual,
            df.null = mod_glance$df.null
          )
      },
      .keep = TRUE
    ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = case_when(
        term == "Titers" ~ "Titer",
        term == "age_ingreso" ~ "Age",
        TRUE ~ term
      ),
      pstar = gtools::stars.pval(p.value),
      `P-value` = insight::format_p(p.value, stars_only = TRUE),
      p.value = print(formatC(signif(p.value, digits = 3), digits = 3, format = "fg", flag = "#")),
      protection.estimate = (1 - estimate) * 100,
      protection.conf.low = (1 - conf.low) * 100,
      protection.conf.high = (1 - conf.high) * 100,
      OR = paste_ci(estimate, conf.low, conf.high),
      Protection = paste_ci_round(protection.estimate, protection.conf.high, protection.conf.low),
      Outcome = ifelse(outcome_var == "pt_inf", "Infection", "Symptomatic infection"),
      std.error = round(std.error, digits = 3),
      statistic = round(statistic, digits = 3)
    ) %>%
    select(
      Assay, variant, `Variable` = term, OR, Protection, p.value, pstar, Outcome,
      `P-value`, std.error, statistic, df.residual, df.null,
      protection.estimate, protection.conf.low, protection.conf.high
    )
}

# "Infection" outcome
tables3 = create_multivariate_table(
  data = datasero2.long,
  outcome_var = "pt_inf",
  output_file = "Results/Table_S3_Multivariate_analysis_by_infection.docx"
)

# "Symptomatic infection" outcome
tables4 = create_multivariate_table(
  data = datasero2.long,
  outcome_var = "sympt.inf",
  output_file = "Results/Table_S4_Multivariate_analysis_by_symptomatic_infection.docx"
)


outcomes <- c("pt_inf", "sympt.inf")


# glm_models_raw_named <- list()
# 
# for (outcome in c("pt_inf", "sympt.inf")) {
#   # Prepare grouped data
#   df_list <- datasero2.long %>%
#     group_by(Assay, variant) %>%
#     group_split()
#   
#   # Get names in format Assay_Variant
#   names_list <- datasero2.long %>%
#     group_by(Assay, variant) %>%
#     group_keys() %>%
#     transmute(name = paste(Assay, variant, sep = " | ")) %>%
#     pull(name)
#   
#   # Run model and gtsummary::tbl_regression with headers
#   glm_models_raw_named[[outcome]] <- Map(function(df, name_label) {
#     model <- glm(
#       as.formula(paste0(outcome, " ~ Titers + any.vax + any.inf + age_ingreso")),
#       data = df,
#       family = binomial(link = "logit")
#     )
#     
#     gtsummary::tbl_regression(model, exponentiate = TRUE,
#                               label = list(
#                                 any.vax ~ "Any prior vax",
#                                 any.inf ~ "Any prior infection",
#                                 age_ingreso ~ "Age",
#                                 Titers ~ "Titer"),
#                               pvalue_fun = ~style_sigfig(., digits = 3)) %>%
#       modify_column_unhide(column = std.error) %>%
#       modify_column_unhide(column = statistic) %>%
#       modify_spanning_header(everything() ~ paste0("Outcome: ",  ifelse(outcome == "pt_inf", "Infection", "Symptomatic infection"), " | ", name_label))
#   },
#   df = df_list,
#   name_label = names_list
#   )
#   
#   names(glm_models_raw_named[[outcome]]) <- names_list
# }


########################################################################################################################################3
# Tables S5-8  Sensitivity analyses  Robust CIs-------------------------------------------------------------------------------------------------------------------------------------------------------------- 
#######################################################################################################################################3

# Define the function
create_multivariate_table_gee <- function(data, outcome_var, output_file) {
  
  # Perform multivariate logistic regression with GEE and extract model + residuals
  table_results <- data %>%
    group_by(Assay, variant) %>%
    mutate(
      hh_code = as.factor(hh_code),
      pt_inf = as.numeric(as.character(pt_inf)),
      sympt.inf = as.numeric(as.character(sympt.inf))
    ) %>%
    arrange(hh_code) %>%
    group_modify(
      ~ {
        model <- geepack::geeglm(
          as.formula(paste(outcome_var, "~ Titers + any.vax + any.inf + age_ingreso")),
          data = .x,
          id = .x$hh_code,
          family = binomial(link = "logit"),
          corstr = "exchangeable"
        )
        
        mod_tidy <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
        mod_glance <- broom::glance(model)
        
        mod_tidy %>%
          mutate(
            df.residual = mod_glance$df.residual,
            alpha = mod_glance$alpha,
            n.clusters = mod_glance$n.clusters
          )
      },
      .keep = TRUE
    ) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "Titers" ~ "Titer",
      term == "age_ingreso" ~ "Age",
      TRUE ~ term
    ),
    pstar = gtools::stars.pval(p.value),
    `P-value` = insight::format_p(p.value, stars_only = TRUE),
    p.value = print(formatC(signif(p.value, digits = 3), digits = 3, format = "fg", flag = "#")),
    protection.estimate = (1 - estimate) * 100,
    protection.conf.low = (1 - conf.low) * 100,
    protection.conf.high = (1 - conf.high) * 100,
    OR = paste_ci(estimate, conf.low, conf.high),
    Protection = paste_ci_round(protection.estimate, protection.conf.high, protection.conf.low),
    Outcome = ifelse(outcome_var == "pt_inf", "Infection", "Symptomatic infection"),
    std.error = round(std.error, digits = 3),
    statistic = round(statistic, digits = 3),
    alpha = round(alpha, digits = 3)
  ) %>%
  select(
    Assay, variant, `Variable` = term, OR, Protection, p.value, pstar, Outcome,
    `P-value`, std.error, statistic, df.residual, alpha, n.clusters,
    protection.estimate, protection.conf.low, protection.conf.high
  )
}




# "Infection" outcome
tables5 = create_multivariate_table_gee(
  data = datasero2.long.noNA,
  outcome_var = "pt_inf",
  output_file = "Results/Table_S5_Multivariate_analysis_by_infection_robust.docx"
)


# "symptomatic" outcome
tabless6 = create_multivariate_table_gee(
  data = datasero2.long.noNA,
  outcome_var = "sympt.inf",
  output_file = "Results/Table_S6_Multivariate_analysis_by_symp_inf_robust.docx"
)


# gee_models_raw_named <- list()
# 
# for (outcome in c("pt_inf", "sympt.inf")) {
#   # Prepare data
#   df_list <- datasero2.long.noNA %>%
#     mutate(
#       hh_code = as.factor(hh_code),
#       pt_inf = as.numeric(as.character(pt_inf)),
#       sympt.inf = as.numeric(as.character(sympt.inf))
#     ) %>%
#     group_by(Assay, variant)
#     group_split()  %>%
#     arrange(hh_code) 
#   
#   # Get group names
#   names_list <- datasero2.long.noNA %>%
#     group_by(Assay, variant) %>%
#     group_keys() %>%
#     transmute(name = paste(Assay, variant, sep = " | ")) %>%
#     pull(name)
#   
#   # Use Map to pair data frames with names
#   gee_models_raw_named[[outcome]] <- Map(function(df, name_label) {
#     model <- geeglm(
#       as.formula(paste0(outcome, " ~ Titers + any.vax + any.inf + age_ingreso")),
#       data = df,
#       id = df$hh_code,
#       family = binomial(link = "logit"),
#       corstr = "exchangeable"
#     )
#     
#     gtsummary::tbl_regression(model, exponentiate = TRUE,
#                               label = list(
#                                 any.vax ~ "Any prior vax",
#                                 any.inf ~ "Any prior infection",
#                                 age_ingreso ~ "Age",
#                                 Titers ~ "Titer"),
#                               pvalue_fun = ~style_sigfig(., digits = 3)
#                               ) %>%
#       modify_column_unhide(column = std.error) %>%
#       modify_column_unhide(column = statistic) %>%
#       modify_spanning_header(everything() ~ paste0("Outcome: ",  ifelse(outcome == "pt_inf", "Infection", "Symptomatic infection"), " | ", name_label))
#   },
#   df = df_list,
#   name_label = names_list
#   )
#   
#   names(gee_models_raw_named[[outcome]]) <- names_list
# }






create_multivariate_table_gee2 <- function(data, outcome_var, output_file) {
  
  # Perform multivariate logistic regression with GEE and format results
  table_results <- data %>%
    group_by(Assay, variant) %>%
  mutate(hh_code = as.factor(hh_code),
         pt_inf = as.numeric(as.character(pt_inf)),
         sympt.inf = as.numeric(as.character(sympt.inf)),
         time.since.exp = pmin(time.since.inf, time.since.vax, na.rm = T)) %>% # get minimum time since exposure between vaccination and infection
    filter(!is.na(time.since.exp)) %>% 
    arrange(hh_code) %>%
    group_modify(
      ~ {
        model <- geepack::geeglm(
          as.formula(paste(outcome_var, "~ Titers + any.vax + any.inf + age_ingreso + time.since.exp")),
          data = .x,
          id = .x$hh_code,
          family = binomial(link = "logit"),
          corstr = "exchangeable"
        )
        
        mod_tidy <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
        mod_glance <- broom::glance(model)
        
        mod_tidy %>%
          mutate(
            df.residual = mod_glance$df.residual,
            alpha = mod_glance$alpha,
            n.clusters = mod_glance$n.clusters
          )
      },
      .keep = TRUE
    ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      term = case_when(
        term == "Titers" ~ "Titer",
        term == "age_ingreso" ~ "Age",
        TRUE ~ term
      ),
      pstar = gtools::stars.pval(p.value),
      `P-value` = insight::format_p(p.value, stars_only = TRUE),
      p.value = print(formatC(signif(p.value, digits = 3), digits = 3, format = "fg", flag = "#")),
      protection.estimate = (1 - estimate) * 100,
      protection.conf.low = (1 - conf.low) * 100,
      protection.conf.high = (1 - conf.high) * 100,
      OR = paste_ci(estimate, conf.low, conf.high),
      Protection = paste_ci_round(protection.estimate, protection.conf.high, protection.conf.low),
      Outcome = ifelse(outcome_var == "pt_inf", "Infection", "Symptomatic infection"),
      std.error = round(std.error, digits = 3),
      statistic = round(statistic, digits = 3),
      alpha = round(alpha, digits = 3)
    ) %>%
    select(
      Assay, variant, `Variable` = term, OR, Protection, p.value, pstar, Outcome,
      `P-value`, std.error, statistic, df.residual, alpha, n.clusters,
      protection.estimate, protection.conf.low, protection.conf.high
    )
}







tables7 = create_multivariate_table_gee2(
  data = datasero2.long.noNA,
  outcome_var = "pt_inf",
  output_file = "Results/Table_S7_Multivariate_analysis_by_infection_robust_time_since_exposure.docx"
)


# "symptomatic" outcome
tables8= create_multivariate_table_gee2(
  data = datasero2.long.noNA,
  outcome_var = "sympt.inf",
  output_file = "Results/Table_S8_Multivariate_analysis_by_symp_inf_robust_robust_time_since_exposure.docx"
)


# gee_models_with_time_named <- list()
# 
# for (outcome in c("pt_inf", "sympt.inf")) {
#   # Prepare data
#   df_filtered <- datasero2.long.noNA %>%
#     mutate(
#       hh_code = as.factor(hh_code),
#       pt_inf = as.numeric(as.character(pt_inf)),
#       sympt.inf = as.numeric(as.character(sympt.inf)),
#       time.since.exp = pmin(time.since.inf, time.since.vax, na.rm = TRUE)
#     ) %>%
#     filter(!is.na(time.since.exp)) %>%
#     arrange(hh_code) 
#   
#   # Group data
#   df_list <- df_filtered %>%
#     group_by(Assay, variant) %>%
#     group_split() %>%
#     arrange(hh_code) 
#   
#   # Create names like "sVNT_BA.1"
#   names_list <- df_filtered %>%
#     group_by(Assay, variant) %>%
#     group_keys() %>%
#     transmute(name = paste(Assay, variant, sep = " | ")) %>%
#     pull(name)
#   
#   # Map through each group
#   gee_models_with_time_named[[outcome]] <- Map(function(df, name_label) {
#     model <- geeglm(
#       as.formula(paste0(outcome, " ~ Titers + any.vax + any.inf + age_ingreso + time.since.exp")),
#       data = df,
#       id = df$hh_code,
#       family = binomial(link = "logit"),
#       corstr = "exchangeable"
#     )
#     
#     gtsummary::tbl_regression(model, exponentiate = TRUE,
#                               label = list(
#                                 any.vax ~ "Any prior vax",
#                                 any.inf ~ "Any prior infection",
#                                 age_ingreso ~ "Age",
#                                 Titers ~ "Titer",
#                                 time.since.exp ~ "Time since last exposure"),
#                               pvalue_fun = ~style_sigfig(., digits = 3)) %>%
#       modify_column_unhide(column = std.error) %>%
#       modify_column_unhide(column = statistic) %>%
#       modify_spanning_header(everything() ~ paste0("Outcome: ", ifelse(outcome == "pt_inf", "Infection", "Symptomatic infection"), " | ", name_label))
#   },
#   df = df_list,
#   name_label = names_list
#   )
#   
#   names(gee_models_with_time_named[[outcome]]) <- names_list
# }
# 



###########################################################################3
# Figure S8  -------------------------------------------------------------------
###########################################################################3
 


# overall analysis
median.assay <- datasero2.long.noNA %>%
  group_by(Assay) %>%
  summarize(median = median(Titers, na.rm = T),
            lower = quantile(Titers, probs = c(0.25), na.rm = T),
            upper = quantile(Titers, probs = c(0.75), na.rm = T),
            estimate = paste_ci_round(4^median, 4^lower, 4^upper))



# Create the violin plot with facet_wrap
figS8A = ggplot(datasero2.long.noNA %>% left_join(median.assay), aes(x = factor(1, label = ""), y = Titers)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.7) +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(. ~ Assay, ncol = 2) +  # Facet by Assay and variant
  theme_classic(base_size = 18) +
  labs(
    x = "",
    y = "Titers",
    title = "A") +
  # Annotate p-values inside each facet
  geom_text(aes(label = estimate), x = 1, y = 12, vjust = 1, size = 6) +
  scale_y_continuous(limits = c(-1,13)) 

  
figS8B =  bind_rows(tables3, tables4) %>%
   filter(Variable == "Titer") %>%
   mutate(Assay = fct_rev(Assay)) %>%
   ggplot(aes(x = Assay, y = protection.estimate, ymin = protection.conf.low, ymax = protection.conf.high)) +
   geom_point(aes(shape = `P-value`), size = 4) +  # Add points for estimates
   geom_segment(aes(x = Assay, xend = Assay, y = protection.conf.low, yend =  protection.conf.high), size =2) +  # Add CI segments
   geom_hline(yintercept = 0, linetype = "dashed", color = "blue2") +  # Add dashed line at y = 0
   facet_wrap(variant ~ Outcome, ncol = 2) +  # Facet by variant
   coord_flip() +  # Flip coordinates for horizontal orientation
   labs(title = "B", x = "Assay", y = "Protection (%)") +
   theme_classic2(base_size = 18) +
   geom_text(aes(label = Protection), nudge_x = 0.4, nudge_y = 0.4)


 ggsave(grid.arrange(figS8A, figS8B, ncol = 2, widths = c(1, 2)),filename = "Results/Figure S6 - Linear protection.svg", width = 22, height = 10)
 
 
 ###########################################################################3
 # Table S9-10 -------------------------------------------------------------------
 ###########################################################################3
 
 # Define the function
 create_predictions <- function(data, outcome_var, method, newdata_params, k = 3, output = "predictions") {
   if (method == "GLM") {
     predictions <- data %>%
       group_by(Assay, variant) %>%
       group_modify(
         ~ augment(
           glm(
             as.formula(paste(outcome_var, "~ Titers + age_ingreso + any.vax + any.inf")),
             data = .,
             family = binomial(link = "logit")
           ),
           type.predict = "response",
           se_fit = TRUE,
           newdata = tibble(
             Titers = seq(
               from = 0,
               to = max(.$Titers, na.rm = TRUE),
               length.out = newdata_params$length_out
             ),
             any.inf = factor(
               newdata_params$any_inf,
               levels = c("Uninfected", "Infected")
             ),
             any.vax = factor(
               newdata_params$any_vax,
               levels = c("Unvaccinated", "Vax")
             ),
             age_ingreso = newdata_params$age_ingreso
           )
         )
       ) %>%
       mutate(
         upper = (((1 - .fitted) + (2 * .se.fit)) * 100),
         lower = (((1 - .fitted) - (2 * .se.fit)) * 100),
         .fitted = ((1 - .fitted) * 100)
       )
   } else if (method == "GAM") {
     predictions <- data %>%
       group_by(Assay, variant) %>%
       group_modify(
         ~ augment(
           gam(
             as.formula(paste(outcome_var, "~ s(Titers, k =", k, ", bs = 'cr') + age_ingreso + any.vax + any.inf")),
             data = .,
             family = binomial(link = "logit")
           ),
           type.predict = "response",
           se_fit = TRUE,
           newdata = tibble(
             Titers = seq(
               from = 0,
               to = max(.$Titers, na.rm = TRUE),
               length.out = newdata_params$length_out
             ),
             any.inf = factor(
               newdata_params$any_inf,
               levels = c("Uninfected", "Infected")
             ),
             any.vax = factor(
               newdata_params$any_vax,
               levels = c("Unvaccinated", "Vax")
             ),
             age_ingreso = newdata_params$age_ingreso
           )
         )
       ) %>%
       mutate(
         upper = (((1 - .fitted) + (2 * .se.fit)) * 100),
         lower = (((1 - .fitted) - (2 * .se.fit)) * 100),
         .fitted = ((1 - .fitted) * 100)
       )
   } else {
     stop("Invalid method. Choose either 'GLM' or 'GAM'.")
   }
   
   return(predictions)
 }
 
 # Define parameters for predictions
 newdata_params <- list(
   length_out = 100000,
   any_inf = "Infected",
   any_vax = "Vax",
   age_ingreso = 27 # average age in the dataset which is adult age
 )
 
 # Generate predictions
 glm_predictions_infection <- create_predictions(
   data = datasero2.long,
   outcome_var = "pt_inf",
   method = "GLM",
   newdata_params = newdata_params
 )
 
 gam_predictions_infection <- create_predictions(
   data = datasero2.long,
   outcome_var = "pt_inf",
   method = "GAM",
   newdata_params = newdata_params
 )
 
 glm_predictions_symptomatic <- create_predictions(
   data = datasero2.long,
   outcome_var = "sympt.inf",
   method = "GLM",
   newdata_params = newdata_params
 )
 
 gam_predictions_symptomatic <- create_predictions(
   data = datasero2.long,
   outcome_var = "sympt.inf",
   method = "GAM",
   newdata_params = newdata_params
 )
 
 # Combine all predictions
 bind.allpredictions.plot <- bind_rows(
   glm_predictions_infection %>%
     mutate(method = "GLM", outcome = "Infection"),
   gam_predictions_infection %>%
     mutate(method = "GAM", outcome = "Infection"),
   glm_predictions_symptomatic %>%
     mutate(method = "GLM", outcome = "Symptomatic infection"),
   gam_predictions_symptomatic %>%
     mutate(method = "GAM", outcome = "Symptomatic infection")
 )
 
 # Compute splines for .fitted, lower, and upper. Make it smoother for plotting
 compute_splines <- function(data, value_column) {
   data %>%
     group_by(method, outcome, Assay, variant) %>%
     group_modify(
       ~ as.data.frame(spline(.$Titers, .[[value_column]]))
     ) %>%
     rename(Titers = x, !!value_column := y)
 }
 
 # Compute splines for each column
 splines_list <- list(
   fitted = compute_splines(bind.allpredictions.plot, ".fitted"),
   lower = compute_splines(bind.allpredictions.plot, "lower"),
   upper = compute_splines(bind.allpredictions.plot, "upper")
 )
 
 # Combine the results into a single dataset
 splines <- reduce(splines_list, full_join, by = c("method", "outcome", "Assay", "variant", "Titers"))
 
 
 
 # Define a reusable function to compute thresholds
 compute_thresholds <- function(predictions, method, outcome) {
   predictions %>%
     mutate(.fitted = round(.fitted)) %>%
     filter(.fitted == 50 | .fitted == 80) %>%
     group_by(.fitted, Assay, variant) %>%
     filter(row_number() == 1) %>%
     mutate(
       Titers = 4^Titers,
       method = method,
       outcome = outcome
     )
 }
 
 # Compute thresholds for all prediction sets
 thresholds <- bind_rows(
   compute_thresholds(glm_predictions_infection, "GLM", "Infection"),
   compute_thresholds(gam_predictions_infection, "GAM", "Infection"),
   compute_thresholds(glm_predictions_symptomatic, "GLM", "Symptomatic infection"),
   compute_thresholds(gam_predictions_symptomatic, "GAM", "Symptomatic infection")
 ) %>%
   mutate(
     protection95 = paste_ci_round(.fitted, lower, upper),
     threshold = case_when(
       .fitted == 50 ~ "50%",
       .fitted == 80 ~ "80%"
     ),
     Titers = format(round(Titers), scientific = FALSE)
   ) %>%
   select(Assay, variant, threshold, Titers, protection95, method, outcome)
 
 
 

   thresholds %>%
     filter(
            method == "GAM") %>%
     pivot_wider(id_cols = c(Assay, variant), names_from = c(threshold, outcome), values_from = c(Titers)) %>%
   arrange(variant, Assay) %>%
   ungroup() %>%
   gt::gt(rownames_to_stub = TRUE) %>%
   gt::gtsave(filename =  "Results/Table S9 -Thresholds GAM.docx")
 
 

   thresholds %>%
     filter(
            method == "GLM") %>%
     pivot_wider(id_cols = c(Assay, variant), names_from = c(threshold, outcome), values_from = c(Titers)) %>%
   arrange(variant, Assay) %>%
   ungroup() %>%
   gt::gt(rownames_to_stub = TRUE) %>%
   gt::gtsave(filename =  "Results/Table S10 -Thresholds GLM.docx")
 

 
###########################################################################3
# Figure 1 -------------------------------------------------------------------
###########################################################################3
 
   
 mut.assay =   function(x){x %>% mutate(Assay = case_when(Assay == "BA.1 Neutralization, Log₄[ID₅₀]"~"BA.1 Neutralization",
                                              Assay == "BA.2 Neutralization, Log₄[ID₅₀]"~"BA.2 Neutralization",
                                              Assay == "D614G Neutralization, Log₄[ID₅₀]"~"D614G Neutralization",
                                              Assay == "Spike Binding, Log₄[Titer]"~"Spike Binding",
                                              Assay == "NP Binding, Log₄[Titer]"~"NP Binding"),
                            Assay = factor(Assay, levels = c("BA.1 Neutralization",
                                                             "BA.2 Neutralization",
                                                             "D614G Neutralization",
                                                             "Spike Binding")))}
                                              
   mut.pt = function(x){x %>% mutate(Outcome = case_when(pt_inf == 0 ~ "Uninfected",
                                                         pt_inf == 1 ~ "Infected"),
                                     Outcome = factor(Outcome, levels = c("Uninfected", "Infected")))}
   
   mut.sympt = function(x){x %>% mutate(Outcome = case_when(sympt.inf == 1 ~ "Symptomatic",
                                                            sympt.inf == 0 ~ "Uninfected or asymptomatic"),
                                     Outcome = factor(Outcome, levels = c("Uninfected or asymptomatic", "Symptomatic")))}
   
   splines.plot2 = splines %>% mut.assay() %>% filter(method == "GAM", !Assay %in% c("Spike Binding"))
   datasero2.long.plot2 = datasero2.long %>% mut.assay() %>% filter(!Assay %in% c("Spike Binding"))  
   
                             
   Fig1A = ggplot(splines.plot2 %>% filter(outcome == "Infection"), aes(y = .fitted, x = Titers)) + 
     geom_line() + 
     facet_wrap(variant ~ Assay, ncol = 6, drop = T) + 
     labs(x = "Titers", y = "Protection (%)") +
     geom_hline(yintercept = 80, color = "#004080") + # Dark Navy Blue for 80%
     geom_hline(yintercept = 50, color = "#1F77B4") + # Medium Blue for 50%
     geom_ribbon(aes(ymin = lower, ymax = upper, alpha = 0.5), fill = "lightgray") + # Light Orange ribbon for confidence intervals
     theme_classic2(base_size = 26) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     geom_density(data = datasero2.long.plot2 %>% mut.pt(), aes(fill = Outcome, x = Titers, y = ..density.. * 200), 
                  color = NA, alpha = 0.45, kernel = "optcosine") +
     scale_x_continuous(breaks = seq(0, 10, by = 2),
                        limits = c(0, 10),
                        labels = round(4^(seq(0, 10, by = 2)))) +
     geom_text(data = thresholds %>% mut.assay() %>%
                 filter(outcome == "Infection", method == "GAM", .fitted == 50, !Assay %in% c("Spike Binding")), 
               aes(x = 1.8, y = 54, label = Titers), color = "#1F77B4", size = 7) +
     geom_text(data = thresholds %>% mut.assay() %>%
                 filter(outcome == "Infection", method == "GAM", .fitted == 80, !Assay %in% c("Spike Binding")), 
               aes(x = 1.8, y = 84, label = Titers), color = "#004080", size = 7) + 
     guides(alpha = "none") +
     scale_fill_manual(values = c(Infected = "#D35400", Uninfected = "#229954")) + # Darker colors for density curves
     theme(legend.position = "bottom") + labs(title = "A")
   
   Fig1B = ggplot(splines.plot2 %>% filter(outcome == "Symptomatic infection"), aes(y = .fitted, x = Titers)) + 
     geom_line() + 
     facet_wrap(variant ~ Assay, ncol = 6, drop = T) + 
     labs(x = "Titers", y = "Protection (%)") +
     geom_hline(yintercept = 80, color = "#004080") + # Dark Navy Blue for 80%
     geom_hline(yintercept = 50, color = "#1F77B4") + # Medium Blue for 50%
     geom_ribbon(aes(ymin = lower, ymax = upper, alpha = 0.5), fill = "lightgray") + # Light Red ribbon for confidence intervals
     theme_classic2(base_size = 26) +
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
     geom_density(data = datasero2.long.plot2 %>% mut.sympt(), aes(fill = Outcome, x = Titers, y = ..density.. * 200), 
                  color = NA, alpha = 0.45, kernel = "optcosine") +
     scale_x_continuous(breaks = seq(0, 10, by = 2),
                        limits = c(0, 10),
                        labels = round(4^(seq(0, 10, by = 2)))) +
     geom_text(data = thresholds %>% mut.assay() %>%  
                 filter(outcome == "Symptomatic infection", method == "GAM", .fitted == 50, !Assay %in% c("Spike Binding", "NP Binding")), 
               aes(x = 1.8, y = 54, label = Titers), color = "#1F77B4", size = 7) +
     geom_text(data = thresholds %>% mut.assay() %>%
                 filter(outcome == "Symptomatic infection", method == "GAM", .fitted == 80, !Assay %in% c("Spike Binding", "NP Binding")), 
               aes(x = 1.8, y = 84, label = Titers), color = "#004080", size = 7) + 
     guides(alpha = "none") +
     scale_fill_manual(values = c(Symptomatic = "#C0392B", `Uninfected or asymptomatic` = "#2874A6"))+ # Darker colors for density curves
     theme(legend.position = "bottom") + labs(title = "B")
    
    
    ggsave(plot =  grid.arrange(Fig1A, Fig1B), "Results/Figure 1 - Protection curves.pdf", device = "pdf", width = 22, height = 18)
   
    

###########################################################################3
# Figure S9 - Protection curve for other assays -------------------------------------------------------------------
###########################################################################3
    
    
    splines.plot3 = splines %>% mut.assay() %>% filter(method == "GAM", Assay %in% c("Spike Binding"))
    datasero2.long.plot3 = datasero2.long %>% mut.assay() %>% filter(Assay %in% c("Spike Binding"))  
    
    FigS9A = ggplot(splines.plot3 %>% filter(outcome == "Infection"), aes(y = .fitted, x = Titers)) + 
      geom_line() + 
      facet_wrap(variant ~ Assay, ncol = 6, drop = T) + 
      labs(x = "Titers", y = "Protection (%)") +
      geom_hline(yintercept = 80, color = "#004080") + # Dark Navy Blue for 80%
      geom_hline(yintercept = 50, color = "#1F77B4") + # Medium Blue for 50%
      geom_ribbon(aes(ymin = lower, ymax = upper, alpha = 0.5), fill = "lightgray") + # Light Orange ribbon for confidence intervals
      theme_classic2(base_size = 26) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      geom_density(data = datasero2.long.plot3 %>% mut.pt() %>% filter(!is.na(Titers)), aes(fill = Outcome, x = Titers, y = ..density.. * 200), 
                   color = NA, alpha = 0.45, kernel = "optcosine") +
      scale_x_continuous(breaks = seq(0, 10, by = 2),
                         limits = c(0, 10),
                         labels = round(4^(seq(0, 10, by = 2)))) +
      geom_text(data = thresholds %>% mut.assay() %>%
                  filter(outcome == "Infection", method == "GAM", .fitted == 50, Assay %in% c("Spike Binding")), 
                aes(x = 1.8, y = 54, label = Titers), color = "#1F77B4", size = 7) +
      geom_text(data = thresholds %>% mut.assay() %>%
                  filter(outcome == "Infection", method == "GAM", .fitted == 80, Assay %in% c("Spike Binding")), 
                aes(x = 1.8, y = 84, label = Titers), color = "#004080", size = 7) + 
      guides(alpha = "none") +
      scale_fill_manual(values = c(Infected = "#D35400", Uninfected = "#229954")) + # Darker colors for density curves
      theme(legend.position = "bottom")
    
    FigS9B = ggplot(splines.plot3 %>% filter(outcome == "Symptomatic infection"), aes(y = .fitted, x = Titers)) + 
      geom_line() + 
      facet_wrap(variant ~ Assay, ncol = 6, drop = T) + 
      labs(x = "Titers", y = "Protection (%)") +
      geom_hline(yintercept = 80, color = "#004080") + # Dark Navy Blue for 80%
      geom_hline(yintercept = 50, color = "#1F77B4") + # Medium Blue for 50%
      geom_ribbon(aes(ymin = lower, ymax = upper, alpha = 0.5), fill = "lightgray") + # Light Red ribbon for confidence intervals
      theme_classic2(base_size = 26) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
      geom_density(data = datasero2.long.plot3 %>% mut.sympt() %>% filter(!is.na(Titers)), aes(fill = Outcome, x = Titers, y = ..density.. * 200), # filtering out the guy without titer value 
                   color = NA, alpha = 0.45, kernel = "optcosine") +
      scale_x_continuous(breaks = seq(0, 10, by = 2),
                         limits = c(0, 10),
                         labels = round(4^(seq(0, 10, by = 2)))) +
      geom_text(data = thresholds %>% mut.assay() %>%  
                  filter(outcome == "Symptomatic infection", method == "GAM", .fitted == 50, Assay %in% c("Spike Binding")), 
                aes(x = 1.8, y = 54, label = Titers), color = "#1F77B4", size = 7) +
      geom_text(data = thresholds %>% mut.assay() %>%
                  filter(outcome == "Symptomatic infection", method == "GAM", .fitted == 80, Assay %in% c("Spike Binding")), 
                aes(x = 1.8, y = 84, label = Titers), color = "#004080", size = 7) + 
      guides(alpha = "none") +
      scale_fill_manual(values = c(Symptomatic = "#C0392B",  `Uninfected or asymptomatic` = "#2874A6"))+ # Darker colors for density curves
      theme(legend.position = "bottom")
    
    
    ggsave(plot =  grid.arrange(FigS9A, FigS9B), "Results/Figure S9 - Protection curves for other assays.pdf", device = "pdf", width = 10, height = 18)
    
    
   
   
   
###########################################################################3
# Table S7  -------------------------------------------------------------------
###########################################################################3
    
    extract_mediation_summary <- function (x) { 
      
      clp <- 100 * x$conf.level
      isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) || 
                       (inherits(x$model.y, "glm") && x$model.y$family$family == 
                          "gaussian" && x$model.y$family$link == "identity") || 
                       (inherits(x$model.y, "survreg") && x$model.y$dist == 
                          "gaussian"))
      
      printone <- !x$INT && isLinear.y
      
      if (printone) {
        
        smat <- c(x$d1, x$d1.ci, x$d1.p)
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        
        rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
        
      } else {
        smat <- c(x$d0, x$d0.ci, x$d0.p)
        smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
        smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
        smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
        smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
        smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
        smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
        smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
        smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
        smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
        
        rownames(smat) <- c("ACME (control)", "ACME (treated)", 
                            "ADE (control)", "ADE (treated)", "Total Effect", 
                            "Prop. Mediated (control)", "Prop. Mediated (treated)", 
                            "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
        
      }
      
      colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""), 
                          paste(clp, "% CI Upper", sep = ""), "p-value")
      smat
      
    }
    
    
    

 # Define a function for mediation analysis with dynamic mediator
 run_mediation_analysis <- function(mediator_column) {
  
   
   datasero3$var = pull(datasero3[,mediator_column])
   
   # infected -> NAb -> infection
   n.sims = 1000
   
   m.inf_nab_inf <- lm(var ~ age_ingreso + any.inf + any.vax + variant, data= datasero3)
   
   c.glm.inf_nab_inf <- glm(pt_inf ~ var  + age_ingreso + any.inf + any.vax + variant, data = datasero3 %>% mutate(pt_inf = as.numeric(as.character(pt_inf))), family = poisson(link = "log"))
   
   set.seed(10000)
   mediate_inf_nab_inf <- mediate(m.inf_nab_inf, c.glm.inf_nab_inf, sims=n.sims, treat="any.inf",
                                  mediator="var", boot=T, boot.ci.type = "bca", control.value = "Uninfected", treat.value = "Infected")
   
   # control.value = "Uninfected", treat.value = "Infected"
   
   
   # infected -> NAb -> symptomatic infection
   n.sims = 1000
   
   m.inf_nab_sympt_inf <- lm(var ~ age_ingreso + any.inf + any.vax + variant, data= datasero3)
   
   c.glm.inf_nab_sympt_inf <- glm(sympt.inf ~ var  + age_ingreso + any.inf + any.vax + variant, data = datasero3 %>% mutate(sympt.inf = as.numeric(as.character(sympt.inf))), family = poisson(link = "log"))
   
   set.seed(10000)
   mediate_inf_nab_sympt_inf <- mediate(m.inf_nab_sympt_inf, c.glm.inf_nab_sympt_inf, sims=n.sims, treat="any.inf",
                                        mediator="var", boot=T, boot.ci.type = "bca", control.value = "Uninfected", treat.value = "Infected")
   
   
   # vax -> NAb -> symptomatic infection
   n.sims = 1000
   
   m.vax_nab_sympt_inf <- lm(var ~ age_ingreso + any.inf + any.vax + variant, data= datasero3)
   
   c.glm.vax_nab_sympt_inf <- glm(sympt.inf ~ var + age_ingreso + any.inf + any.vax + variant, data = datasero3 %>% mutate(sympt.inf = as.numeric(as.character(sympt.inf))), family = poisson(link = "log"))
   
   set.seed(10000)
   mediate_vax_nab_sympt_inf <- mediate(m.vax_nab_sympt_inf, c.glm.vax_nab_sympt_inf, sims=n.sims, treat="any.vax",
                                        mediator="var", boot=T, boot.ci.type = "bca", control.value = "Unvaccinated", treat.value = "Vax")
   
   
   
   
   # vax -> NAb -> infection
   n.sims = 1000
   
   m.vax_nab_inf <- lm(var ~ age_ingreso + any.inf + any.vax + variant, data= datasero3)
   
   c.glm.vax_nab_inf <- glm(pt_inf ~ var  + age_ingreso + any.inf + any.vax + variant, data = datasero3 %>% mutate(pt_inf = as.numeric(as.character(pt_inf))), family = poisson(link = "log"))
   
   set.seed(10000)
   mediate_vax_nab_inf <- mediate(m.vax_nab_inf, c.glm.vax_nab_inf, sims=n.sims, treat="any.vax",
                                  mediator="var", boot=T, boot.ci.type = "bca", control.value = "Unvaccinated", treat.value = "Vax")
   
   
   
   
   # Summarize results
   results <- bind_rows(
     extract_mediation_summary(summary(mediate_inf_nab_inf)) %>% as.data.frame() %>%
       rownames_to_column(var = "Term") %>%
       mutate(Exposure = "Prior infection", Mediator = mediator_column, Outcome = "Infection"),
     extract_mediation_summary(summary(mediate_inf_nab_sympt_inf)) %>% as.data.frame() %>%
       rownames_to_column(var = "Term") %>%
       mutate(Exposure = "Prior infection", Mediator = mediator_column, Outcome = "Symptomatic infection"),
     extract_mediation_summary(summary(mediate_vax_nab_inf)) %>% as.data.frame() %>%
       rownames_to_column(var = "Term") %>%
       mutate(Exposure = "Prior Vax", Mediator = mediator_column, Outcome = "Infection"),
     extract_mediation_summary(summary(mediate_vax_nab_sympt_inf)) %>% as.data.frame() %>%
       rownames_to_column(var = "Term") %>%
       mutate(Exposure = "Prior Vax", Mediator = mediator_column, Outcome = "Symptomatic infection")
   ) %>%
     filter(grepl("av|Total", x = Term)) %>%
     filter(!grepl("Prop", x = Term)) %>%
     mutate(
       Protection = round((1 - exp(Estimate)) * 100),
       `Protection Lower` = round((1 - exp(`95% CI Lower`)) * 100),
       `Protection Upper` = round((1 - exp(`95% CI Upper`)) * 100),
       Estimate = exp(Estimate),
       `95% CI Lower` = exp(`95% CI Lower`),
       `95% CI Upper` = exp(`95% CI Upper`)
     ) %>%
     mutate(
       Protection2 = paste_ci_round(Protection,  `Protection Upper`, `Protection Lower`),
       Estimate2 = paste_ci(Estimate, `95% CI Lower`, `95% CI Upper`)
     ) %>%
     group_by(Exposure, Mediator, Outcome) %>%
     mutate(Mediation.Proportion = round((first(Estimate)/(abs(nth(Estimate,2) + (abs(nth(Estimate,3))))))*100),
            Mediation.lower = round((first(`95% CI Lower`)/(abs(nth(`95% CI Lower`,2) + (abs(nth(`95% CI Lower`,3))))))*100),
            Mediation.upper = round((first(`95% CI Upper`)/(abs(nth(`95% CI Upper`,2) + (abs(nth(`95% CI Upper`,3))))))*100)) %>%
     mutate(
       Mediation.Proportion = paste_ci_round(Mediation.Proportion, Mediation.lower, Mediation.upper)
     ) %>%
     select(-`95% CI Lower`, -`95% CI Upper`, -Estimate, -Mediation.lower, -Mediation.upper) %>%
     ungroup()
   
   # Save table
   return(results)
     
 }
 
 
 

 neut.mediation = run_mediation_analysis(mediator_column = "neut.homotypic")
 ancestral.mediation = run_mediation_analysis(mediator_column = "neut.ancestral")
 spike.mediation = run_mediation_analysis(mediator_column = "spike_titer.jv")


bind_rows(neut.mediation, ancestral.mediation, spike.mediation) %>%
  ungroup() %>%
  gt::gt(rownames_to_stub = TRUE) %>%
  gt::gtsave(filename =  "Results/Table S7 - Mediation results.docx")



########################################################################3
# Figure 2 ###################################
#########################################################################4


combined.mediation = bind_rows(neut.mediation, ancestral.mediation, spike.mediation) %>% 
  filter(Exposure == "Prior infection") %>%
  mutate(Term = case_when(Term == "Total Effect" ~ "Total\n Effect",
                          Term == "ADE (average)" ~ "Direct\n Effect",
                          Term == "ACME (average)" ~ "Mediating\n Effect",
                          T ~ Term),
         Term = fct_rev(factor(Term, levels = c("Mediating\n Effect",
                                                "Direct\n Effect",
                                                "Total\n Effect"))),
         `P-value` = insight::format_p(`p-value`,
                                   stars_only = TRUE),
         
         Mediator = factor(Mediator, levels = c("neut.homotypic", "neut.ancestral", "spike_titer.jv"), labels = c("Homotypic Neutralization", "D614G Neutralization", "Spike Binding")),
         Outcome = factor(Outcome, levels = c("Infection", "Symptomatic infection"))
         
         ) 


fig2 = ggplot(combined.mediation, aes(x = Term, y = Protection, ymin = `Protection Lower`, ymax = `Protection Upper`)) +
  geom_point(aes(shape = `P-value`),size = 4) +  # Add points for estimates
  geom_segment(aes(x = Term, xend = Term, y = `Protection Lower`, yend = `Protection Upper`), size =2) +  # Add CI segments
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue2") +  # Add dashed line at y = 0
  facet_wrap(Mediator ~ Outcome, ncol = 2) +  # Facet by variant
  coord_flip() +  # Flip coordinates for horizontal orientation
  labs(title = "", x = "Effect", y = "Protection (%)") +
  theme_classic2(base_size = 15) +
  geom_text(aes(label = Protection2), nudge_x = 0.4, nudge_y = 0.4)
ggsave(filename = "Results/Figure 2 - Mediation.pdf", width = 10, height = 10)
ggsave(fig2, filename = "Results/Figure 2 - Mediation.tif", device = "tif", width = 10, height = 10)
ggsave(fig2, filename = "Results/Figure 2 - Mediation.eps", device = "eps", width = 10, height = 10)



###########################################################################4
# Figure S4 #########################################################
###########################################################################4

# this goes with warnings on scale problems

aggregated.seq.long = aggregated.seq %>%
  dplyr::select(-tot_cases) %>%
  mutate(date = as.Date(paste(year, week, 7, sep="-"), "%Y-%U-%u")) %>% # just for plotting
  relocate(date, .after = week_year) %>%
  pivot_longer(cols = clade_Not_sequenced:clade_BA_1_20, names_to = "Clade") %>%
  filter(Clade != "clade_Not_sequenced") %>%
  mutate(Clade = gsub("clade_", "", Clade),  # Remove "clade_"
         Clade = gsub("_", ".", Clade),
         Clade = factor(case_when(grepl(x = Clade, pattern = "BA\\.1") ~ "BA.1 (Omicron)",
                                  grepl(x = Clade, pattern = "BA\\.2") ~ "BA.2 (Omicron)")))  %>%# Replace "_" with ".")
  filter(date >  ymd("2022-01-01") &  date <  ymd("2022-07-01"), value > 0) %>%
  rename(Subvariant = Clade) 

FigS4A = ggplot(aggregated.seq %>%
             mutate(date = as.Date(paste(year, week, 7, sep="-"), "%Y-%U-%u")) %>%
             filter(date >=  ymd("2022-01-01") &  date <=  ymd("2022-07-01")), 
           aes(date, tot_cases, fill = "NA")) + geom_col() + theme_classic2(base_size = 18) + labs(y = "# of Infections", title = "A", x = "Date") +
  scale_x_date(date_labels = "%b" ,date_breaks = "1 month", limits = c(ymd("2021-12-31"), ymd("2022-07-15")))  +  
  theme(legend.position="bottom") +
  scale_y_continuous(limits = c(0,60))


FigS4B = ggplot(aggregated.seq.long, aes(x = date, y = value, fill = Subvariant)) + geom_col() + theme_classic2(base_size = 18) +
  labs(y = "# of sequences", title = "B", x = "Date") +
  scale_x_date(date_labels = "%b" ,date_breaks = "1 month", limits = c(ymd("2021-12-31"), ymd("2022-07-15"))) +
  theme(legend.position="bottom") +
  scale_fill_manual(values = c("BA.1 (Omicron)" = "#9c954d", "BA.2 (Omicron)" = "#b067a3"))

FigS4C = aggregated.seq.minsa %>%
  mutate(week = as.Date(week)) %>%
  rename(Subvariant = Lineage) %>%
  filter(Subvariant != "Unassigned", week >=  ymd("2022-01-01") &  week <=  ymd("2022-07-01")) %>%
  ggplot(aes(x = week, y = count, fill = Subvariant)) + geom_col() +
  scale_x_date(date_labels = "%b" ,date_breaks = "1 month", limits = c(ymd("2021-12-31"), ymd("2022-07-15"))) +
  theme_classic2(base_size = 18) +
  labs(y = "# of sequences", title = "C",  x = "Date") +  theme(legend.position="bottom") +
  scale_fill_manual(values = c("BA.1 (Omicron)" = "#9c954d", "BA.2 (Omicron)" = "#b067a3"))


FigS4D = ggplot(datasero2 %>%
             filter(!is.na(inf_date)) %>%
             mutate(inf_date = week(inf_date),
                    sequence2 = factor(case_when(grepl(x = sequence2, pattern = "BA\\.1") ~ "BA.1 (Omicron)",
                                                 grepl(x = sequence2, pattern = "BA\\.2") ~ "BA.2 (Omicron)")),
                    inf_date = as.Date(paste(2022, inf_date, 7, sep="-"), "%Y-%U-%u")) %>%
             count(Date = inf_date, Subvariant = sequence2), 
           aes(x = Date, y = n, fill = Subvariant) ) + geom_col() +
  scale_x_date(date_labels = "%b" ,date_breaks = "1 month", limits = c(ymd("2021-12-31"), ymd("2022-07-15"))) +
  theme_classic2(base_size = 18) +
  labs(y = "# of sequences", title = "D",  x = "Date") +  theme(legend.position="bottom") +
  scale_fill_manual(values = c("BA.1 (Omicron)" = "#9c954d", "BA.2 (Omicron)" = "#b067a3")) +
  scale_y_continuous(limits = c(0,60))




ggsave(gridExtra::grid.arrange(FigS4A,FigS4B,FigS4C, FigS4D, ncol = 1), device = "pdf", filename = "Figure S4 - Sequencing data.pdf", width = 8, height = 11)



