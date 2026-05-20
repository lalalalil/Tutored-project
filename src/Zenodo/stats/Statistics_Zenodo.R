library(tidyverse)
library(jsonlite)

# SETTINGS & PATHS 
# Define where to save plots and which files we care about
output_dir <- "outputs_zenodo"
dir.create(output_dir, showWarnings = FALSE)

EXECUTABLE_EXTENSIONS <- c(
  "py", "r", "m", "ipynb", "sbml", "sedml",
  "omex", "cps", "cellml", "jl", "cpp", "c"
)

# DATA EXTRACTION 
extract_zenodo_data <- function() {
  # Search for all metadata files in the current directory
  json_files <- list.files(pattern = "metadata_.*\\.json$", recursive = TRUE, full.names = TRUE)
  
  if (length(json_files) == 0) {
    stop("No JSON files found. Please check your working directory.")
  }
  
  message("Processing ", length(json_files), " metadata files...")
  
  all_data <- map_df(json_files, function(file_path) {
    # Load JSON content
    meta <- fromJSON(file_path)
    
    # Identify the virus/disease name from folder structure
    path_parts <- strsplit(file_path, "[/\\\\]")[[1]]
    zenodo_idx <- which(path_parts == "Zenodo")
    disease_name <- if(length(zenodo_idx) > 0) path_parts[zenodo_idx - 1] else "Unknown"
    
    # Look into the 'data' sibling folder to count actual files
    doi_folder <- dirname(dirname(file_path))
    data_path <- file.path(doi_folder, "data")
    
    files_on_disk <- if(dir.exists(data_path)) list.files(data_path) else character(0)
    exts_on_disk <- tolower(tools::file_ext(files_on_disk))
    
    # Count occurrences for each extension in our target list
    counts <- map(EXECUTABLE_EXTENSIONS, ~ sum(exts_on_disk == .x))
    names(counts) <- EXECUTABLE_EXTENSIONS
    
    # Merge everything into a single row
    bind_cols(
      tibble(
        disease_category = disease_name,
        publication_year = as.numeric(substr(meta$metadata$publication_date, 1, 4)),
        journal = if (!is.null(meta$metadata$journal_title)) meta$metadata$journal_title else "N/A"
      ),
      as_tibble(counts)
    )
  })
  
  return(all_data)
}

# Execute extraction
df <- extract_zenodo_data()

# PLOT
# Global theme for all charts
theme_bio <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Plot 1: Timeline of publications
p1 <- df %>%
  count(publication_year) %>%
  filter(!is.na(publication_year)) %>%
  ggplot(aes(x = publication_year, y = n)) +
  geom_line(color = "#2E86C1", linewidth = 1) +
  geom_point(color = "#1B4F72") +
  labs(title = "Publications Timeline", x = "Year", y = "Count") +
  theme_bio

# Plot 2: Distribution of file types
p2 <- df %>%
  summarise(across(all_of(EXECUTABLE_EXTENSIONS), sum)) %>%
  pivot_longer(everything(), names_to = "extension", values_to = "total") %>%
  filter(total > 0) %>% 
  ggplot(aes(x = reorder(extension, -total), y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = total), vjust = -0.5, size = 3) +
  labs(title = "Files by Extension", x = "Type", y = "Total") +
  theme_bio

# Plot 3: Distribution per Virus/Disease category
p3 <- df %>%
  count(disease_category) %>%
  filter(disease_category != "Unknown") %>% 
  ggplot(aes(x = reorder(disease_category, -n), y = n)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() + # Horizontal for better readability
  labs(title = "Models per Category", x = "Category", y = "Count") +
  theme_bio

# EXPORT 
ggsave(file.path(output_dir, "01_timeline.png"), p1, width = 8, height = 5)
ggsave(file.path(output_dir, "02_extensions.png"), p2, width = 8, height = 5)
ggsave(file.path(output_dir, "03_categories.png"), p3, width = 8, height = 5)
