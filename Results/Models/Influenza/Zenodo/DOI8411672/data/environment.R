set.seed(18)
# Paths for data
sequences_path <- "data/rawData/sorted-smith-et-al_68-02_protein-seqs.fasta"
phylo_path <- "data/rawData/smith-et-al_68-02.rooted.tree"
titer_path <- "data/rawData/smith-2004-formatted-and-sorted-data.csv"

# Load data from paths into variables
sequences <- read.alignment(file = sequences_path, format = "fasta")
tree <- read.nexus(file = phylo_path)
titer_data <- read.titerTable(titer_path)

# Specifications for antigenic cartography analysis
antigenic_dimensions <- 3
antigenic_optimizations <- 1000

# Color function for plots
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# Setting theme_bw() for all plots
ggplot2::theme_set(ggplot2::theme_bw())

# The hope is that this will not be necessary but for now here it is...
# Smith et al. 2004 already gives us accurate cluster information. We are 
# basing A LOT of off the clusters that they have published. The only
# way to recreate their results is to do k-means clustering. This can 
# be done but the accuracy (with respect to what Smith et al. published) 
# of the algorithm depends on the antigenic maps and knowing exactly how many 
# clusters there are.

source('scripts/R/antigenic_mapping.R')

# Read in HI data using Racmacs to get the colors
map_file_path <- system.file("extdata/h3map2004.ace", package = "Racmacs")
smith_map <- read.acmap(map_file_path)

for (i in 1:length(agFill(h3_map))) {
  tmp <- agFill(smith_map)[which(agNames(h3_map)[i] == agNames(smith_map))]
  agFill(h3_map)[i] <- tmp
}

clusters <- c("HK68", "EN72", "VI75",
              "TX77", "BK79", "SI87",
              "BE89", "BE92", "WU95",
              "SY97", "FU02")

clusters_for_color <- c()
for (i in agFill(h3_map)) {
  clusters_for_color <- c(clusters_for_color, 
                          clusters[which(unique(agFill(h3_map)) == i)])
}

sequences$cluster <- clusters_for_color
