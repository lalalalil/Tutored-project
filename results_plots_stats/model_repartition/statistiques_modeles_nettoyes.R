library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)

# 2. Lecture du fichier
# On suppose que le fichier est dans le répertoire de travail
df <- read.csv("statistiques_modeles_nettoyes.csv", check.names = FALSE)

# Remplacer les NA par 0 (si jamais il en reste)
df[is.na(df)] <- 0

# 3. Préparation des données pour l'analyse
# On pivote les données pour avoir un format long (plus facile pour les graphiques)
df_long <- df %>%
  pivot_longer(
    cols = -c(category, model_id),
    names_to = "extension",
    values_to = "count"
  ) %>%
  filter(count > 0) # On ne garde que ce qui existe réellement

# --- ANALYSES ---

# A. Nombre total de fichiers par extension (Top global)
stats_extensions <- df_long %>%
  group_by(extension) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total))

print("Top des extensions les plus utilisées :")
print(head(stats_extensions))

# B. Répartition des formats par Maladie (Category)
stats_par_maladie <- df_long %>%
  group_by(category, extension) %>%
  summarise(total = sum(count)) %>%
  ungroup()

# --- VISUALISATIONS ---

# 1. Graphique des extensions les plus fréquentes
p1 <- ggplot(stats_extensions, aes(x = reorder(extension, -total), y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Volume total de fichiers par type d'extension",
       x = "Extensions",
       y = "Nombre de fichiers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Heatmap de la présence des formats par catégorie
# Idéal pour voir si le COVID utilise plus de .omex que la Malaria par exemple
p2 <- ggplot(stats_par_maladie, aes(x = category, y = extension, fill = total)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Densité des formats de fichiers par maladie",
       x = "Maladies",
       y = "Extensions",
       fill = "Quantité") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. Nombre de modèles par catégorie
p3 <- df %>%
  count(category) %>%
  ggplot(aes(x = reorder(category, -n), y = n)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Nombre de modèles uniques par catégorie",
       x = "Catégorie",
       y = "Nombre de modèles")

# Affichage des graphiques
print(p1)
print(p2)
print(p3)

# 4. Export d'un résumé rapide
# Nombre moyen de fichiers par modèle dans chaque catégorie
resume_moyen <- df_long %>%
  group_by(category, model_id) %>%
  summarise(total_files = sum(count)) %>%
  group_by(category) %>%
  summarise(moyenne_fichiers_par_modele = mean(total_files))

write.csv(resume_moyen, "resume_moyen_par_categorie.csv", row.names = FALSE)


library(tidyverse)

# 1. Chargement et nettoyage (assure-toi que df est bien chargé avant)
df <- read.csv("statistiques_modeles_nettoyes.csv", check.names = FALSE)
df[is.na(df)] <- 0

# 2. Transformation en format long
df_long <- df %>%
  pivot_longer(
    cols = -c(category, model_id),
    names_to = "extension",
    values_to = "count"
  ) %>%
  filter(count > 0)

# 3. Calcul des totaux par extension pour l'affichage
stats_extensions <- df_long %>%
  group_by(extension) %>%
  summarise(total = sum(count)) %>%
  arrange(desc(total))

# 4. Graphique avec effectifs affichés
ggplot(stats_extensions, aes(x = reorder(extension, -total), y = total)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  # Ajout des effectifs au-dessus des barres
  geom_text(aes(label = total), vjust = -0.5, size = 3.5) + 
  theme_minimal() +
  labs(
    title = "Nombre total de fichiers par type d'extension",
    x = "Extensions de fichiers",
    y = "Nombre total"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



