# List all the packages that will be used
packages <- c("phytools",
              "phangorn",
              "ape",
              "ggplot2",
              "Racmacs",
              "ggacmap",
              "ggpubr",
              "plotly",
              "factoextra",
              "cowplot",
              "seqinr",
              "plyr",
              "reshape2",
              "RColorBrewer",
              "Peptides",
              "cdata", 
              "GGally",
              "xtable",
              "gridExtra",
              "tidyverse",
              "ggraph",
              "igraph",
              "tidygraph",
              "neuralnet",
              "caret",
              "adephylo")
# List any BioConductor packages that may be needed
bioconductor.packages <- c()

packages.to.install <- packages[!(packages %in% installed.packages()[,'Package'])]
bioC.to.install <- bioconductor.packages[!(bioconductor.packages %in% installed.packages()[,'Package'])]

if (length(packages.to.install)) {
  install.packages(packages.to.install, 
                   repos = "https://cloud.r-project.org", 
                   dependencies = TRUE)
}

# Check to see if there are any packages to be installed using BioConductor
if (length(bioC.to.install)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(bioC.to.install)
}

# This will load all the packages listed above...
invisible(lapply(packages, library, character.only = TRUE))
invisible(lapply(bioconductor.packages, library, character.only = TRUE))

## installation
# for the up to date version
devtools::install_github("yulab-smu/ggacmap")
library(ggacmap)
