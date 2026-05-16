#!/usr/bin/env Rscript
# NCBI Analysis: File Extensions and SBML Models
# Generates: (1) File extensions plot, (2) PMC6620235 SBML organisms CSV

library(tidyverse)

project_dir <- path.expand("~/Desktop/study/M1/S2/Projet\ tutoré\ Master\ 1\ BBS/ptut-final/ptut")
output_dir <- file.path(project_dir, "ncbi_analysis_outputs")
dir.create(output_dir, showWarnings = FALSE)

# ============================================================================
# PART 1: FILE EXTENSIONS DISTRIBUTION
# ============================================================================

message("Reading disease summary files...")

df_all <- tibble()
diseases <- c("chikungunya_files", "covid_files", "dengue_files", "hiv_files",
              "influenza_files", "mpox_files", "tuberculosis_files", "west_nile_files")

for (disease_folder in diseases) {
  csv_path <- file.path(project_dir, "downloaded_ncbi_models 2", disease_folder, "disease_summary.csv")
  if (file.exists(csv_path)) {
    df_disease <- read.csv(csv_path, check.names = FALSE) %>%
      as_tibble() %>%
      mutate(disease = str_remove(disease_folder, "_files$"))
    df_all <- bind_rows(df_all, df_disease)
  }
}

# Count extensions
extension_stats <- tibble()

data_cols <- colnames(df_all) %>% str_subset("^data_") %>% str_remove("^data_")
for (ext in data_cols) {
  col_name <- paste0("data_", ext)
  if (col_name %in% colnames(df_all)) {
    values <- df_all[[col_name]]
    values[is.na(values)] <- 0
    total <- sum(as.numeric(values), na.rm = TRUE)
    if (total > 0) {
      extension_stats <- bind_rows(extension_stats,
        tibble(extension = ext, category = "DATA", count = total))
    }
  }
}

metadata_cols <- colnames(df_all) %>% str_subset("^metadata_") %>% str_remove("^metadata_")
for (ext in metadata_cols) {
  col_name <- paste0("metadata_", ext)
  if (col_name %in% colnames(df_all)) {
    values <- df_all[[col_name]]
    values[is.na(values)] <- 0
    total <- sum(as.numeric(values), na.rm = TRUE)
    if (total > 0) {
      extension_stats <- bind_rows(extension_stats,
        tibble(extension = ext, category = "METADATA", count = total))
    }
  }
}

# Create plot
p <- extension_stats %>%
  arrange(desc(count)) %>%
  slice(1:15) %>%
  mutate(label = paste0(str_remove(extension, "_"), " (", category, ")")) %>%
  ggplot(aes(x = reorder(label, count), y = count, fill = category)) +
  geom_col() +
  geom_text(aes(label = count), hjust = -0.2, size = 3) +
  scale_fill_manual(values = c("DATA" = "#3498DB", "METADATA" = "#2ECC71")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  coord_flip() +
  labs(title = "File Extensions in NCBI Downloads",
       subtitle = "Top 15 file types - mostly metadata (PDFs, images), not computational models",
       x = NULL, y = "Count", fill = "Type") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13))

ggsave(file.path(output_dir, "01_file_extensions.png"), p, width = 10, height = 7, dpi = 300, bg = "white")
message("✓ Saved: 01_file_extensions.png")

# ============================================================================
# PART 2: PMC6620235 SBML ORGANISMS
# ============================================================================

manifest_path <- file.path(project_dir, "data_manifest.csv")

if (file.exists(manifest_path)) {
  manifest <- read.csv(manifest_path) %>% as_tibble()

  pmc_data <- manifest %>%
    filter(pmcid == "PMC6620235", category == "sbml_xml") %>%
    mutate(organism = str_extract(filename, "(?<=sbml3_|sbm3_)[^.]+"),
           organism = str_replace_all(organism, "_", " ")) %>%
    select(organism) %>%
    distinct() %>%
    arrange(organism) %>%
    mutate(row_num = row_number(),
           type = case_when(
             str_detect(organism, "tuberculosis|Mtuberculosis") ~ "Tuberculosis",
             str_detect(organism, "influenzae|influenza") ~ "Influenza",
             TRUE ~ "Other"))

  if (nrow(pmc_data) > 0) {
    # Save as CSV
    output_csv <- pmc_data %>% select(row_num, organism, type)
    colnames(output_csv) <- c("No.", "Organism", "Type")
    write.csv(output_csv, file.path(output_dir, "02_pmc6620235_sbml_organisms.csv"), row.names = FALSE)

    message(sprintf("✓ PMC6620235 contains %d unique SBML organisms", nrow(pmc_data)))
    message(sprintf("  - Tuberculosis: %d", sum(pmc_data$type == "Tuberculosis")))
    message(sprintf("  - Influenza: %d", sum(pmc_data$type == "Influenza")))
    message("✓ Saved: 02_pmc6620235_sbml_organisms.csv")
  }
} else {
  message("✗ data_manifest.csv not found")
}

message("\nDone. Output files in: ncbi_analysis_outputs/")
