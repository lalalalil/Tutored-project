# 1. Load necessary libraries
library(tidyverse)
library(janitor)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

# 2. Load the CSV file
# Ensure the file name matches your exported CSV
df_approach <- read.csv("modelling_approaches_summary.csv") %>% 
  as_tibble()

# 3. Clean-up: Remove generic or unspecified categories for better visualization
df_clean <- df_approach %>%
  filter(!modelling_approach %in% c("Not specified", "Other", "Unknown", "N/A"))

# --- ANALYSIS 1: GLOBAL DISTRIBUTION ---
stats_global <- df_clean %>%
  count(modelling_approach) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  arrange(desc(n))

print("Global Distribution of Modelling Approaches:")
print(stats_global)

# Global bar chart
plot_global <- ggplot(stats_global, aes(x = reorder(modelling_approach, n), y = n)) +
  geom_col(fill = "#5DADE2") +
  coord_flip() +
  labs(
    title = "Overall Distribution of Modelling Approaches",
    subtitle = "Aggregated data from all infectious disease categories",
    x = "Modelling Approach", 
    y = "Number of Models"
  ) +
  theme_minimal()

# --- ANALYSIS 2: DISTRIBUTION BY DISEASE CATEGORY ---
stats_by_disease <- df_clean %>%
  group_by(disease_category, modelling_approach) %>%
  summarise(count = n(), .groups = "drop")

# Stacked bar chart (Normalized to 100% to compare proportions)
plot_disease <- ggplot(stats_by_disease, aes(x = disease_category, y = count, fill = modelling_approach)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  labs(
    title = "Modelling Approaches Proportion by Disease",
    subtitle = "Relative comparison of mathematical methods used",
    x = "Disease Category", 
    y = "Percentage of Models (%)",
    fill = "Approach Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom", 
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 2)) # Organizes legend in 2 columns for readability

# --- DISPLAY OUTPUTS ---
print(plot_global)
print(plot_disease)

# Optional: Save the plots
# ggsave("global_approaches.png", plot_global, width = 8, height = 6)
# ggsave("approaches_by_disease.png", plot_disease, width = 10, height = 8)