###################################################################################################################
# Create heat maps for linear separable Genes based on Day 0 CPM
#
# John Hanley
# October 23, 2019
# 
###################################################################################################################

# Get the current working directory
CurDir <- getwd()

# Read in the linear Separable CPMs
Lcpm <- read.csv(paste(CurDir, "/Day_0_CPM_Limits_EitherRash_MaxP.csv", sep = ""),
                 stringsAsFactors = FALSE)
median_CPM <- read.csv(paste(CurDir, "/Median_CPM_Sep20_2019.csv", sep = ""), 
                       row.names="gene_id", stringsAsFactors = FALSE)
# select only the day 0 CPMs
Day0Medcpm <- median_CPM[, seq(1, 33, 3)]

# For the selected GeneIDs, extract the gene names
Gene_Name_IDs_Map <- read.csv(paste(CurDir, "/GeneNames_GeneIDs_GeneBiotype_Description_1_to_1_Oct10_2019.csv", sep = ""),
                              stringsAsFactors = FALSE)

GeneNameV <- Lcpm$Gene_Name
# Create a Mask for the genes that were selected
KeepMask <- rep(FALSE, nrow(Day0Medcpm))

for (i in 1:nrow(Lcpm)){
  # Extract the ith Gene_ID and insert the Gene Name if it exists
  ithID <- Lcpm$Gene_Name[i]
  Mask <- Gene_Name_IDs_Map$Gene_IDs == ithID
  
  if (!is.na(Gene_Name_IDs_Map$Gene_Names[Mask])){
    # Then insert the Gene Name
    GeneNameV[i] <- Gene_Name_IDs_Map$Gene_Names[Mask]
  }
  # Find where the ith ID matches the Day0Medcpm
  ID <- which(ithID == row.names(Day0Medcpm))
  # Now set the ID to TRUE in the KeepMask
  KeepMask[ID] <- TRUE
}

SelDay0cpm <- Day0Medcpm[KeepMask, ]
# Set the rownames and column names
row.names(SelDay0cpm) <- GeneNameV
colnames(SelDay0cpm) <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K')
#colnames(SelDay0cpm) <- c('Rash', 'Rash', 'Rash', 'No Rash', 'Rash', 'Rash', 'Rash', 'Rash', 'No Rash', 'Rash', 'Rash')


###############################################################################################################
# Now Create a matrix with the correlation of genes for each subject

# Convert the Thresh genes to a matrix
HeatMat <- as.matrix(SelDay0cpm)

# Take the log2 of the HeatMat
Log2HeatMat <- log2(HeatMat)
# Replace the -Inf with -4
Mask <- Log2HeatMat == -Inf
Log2HeatMat[Mask] <- -4


library(heatmaply)
library(RColorBrewer)
curcolors <-  colorRampPalette(brewer.pal(9,"YlOrRd"))
#Keytit <- expression('Log'['2']*'(Median CPM)')
Keytit <- 'Log2(Median CPM)'
heatmaply(Log2HeatMat, main = 'Linear Separable Genes (Day 0 CPM) Rash vs No Rash',
          colors = curcolors, fontsize_row = 4, key.title = Keytit, limits = c(-4, 12),
          file = "Z://JohnHanleyWork/DengueProject/Fall_2019/Figures/HeatMap_Day0cpm_LinearSeparableRashGenes.html")


################################################################################################################
# Now remove any MSTRGs and remove ORs

# Extract the row names
CurGenes <- row.names(Log2HeatMat)

# Now remove any row name that has a MSTRG or ENS
KeepMask <- rep(TRUE, length(CurGenes))
RemoveIDs <- grep('MSTRG.', CurGenes)
KeepMask[RemoveIDs] <- FALSE 
# Now remove any MSTRG Genes
NewHeatMat <- Log2HeatMat[KeepMask,]

# Extract the row names
CurGenes <- row.names(NewHeatMat)
# Now find the genes with Or
OrIDs <- grep('Or', CurGenes)

for (i in 1:length(OrIDs)){
  # Remove any Genes after the or
  CurGenes[OrIDs[i]] <- strsplit(CurGenes[OrIDs[i]], " Or")[[1]][1]
}

# Now rename the row names
row.names(NewHeatMat) <- CurGenes

# To Make the plot more enjoyable round to the thousandths place
NewHeatMat <- round(NewHeatMat, 3)

heatmaply(NewHeatMat, main = 'Linear Separable Genes (Day 0 CPM) Rash vs No Rash',
          colors = curcolors, fontsize_row = 4, key.title = Keytit, limits = c(-4, 12),
          file = "Z://JohnHanleyWork/DengueProject/Fall_2019/Figures/HeatMap_Day0cpm_LinearSeparableRashGenes_Dec19_2019.html")
