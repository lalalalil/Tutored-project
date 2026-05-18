##################################################################################################################
# Function to perform logical conjunction analysis (lca)
# 
# Created by: John Hanley
# December 7, 2020
# Last updated: December 7, 2020
#
# Inputs:
# logFC = a dataframe of the log fold change for each subject and gene. The columns are the subjects and the rows
#         are the genes.
# pVal = a dataframe of the p-value for each subject and gene. The columns are the subjects and the rows are
#        the genes.
# FDR = a dataframe of the false discovery rate (FDR) for each subject and gene. The columns are the subjects 
#       and the rows are the genes.
# thresholds = a dataframe that contains the thresholds that are being used determine gene regulation across all
#              subjects. The variables for the threshold should be entered $logFC, $pVal, and $FDR with one 
#              value for each variable.
#
# Outputs: 
# RegGenes = a dataframe where the first two variables are logical vectors of whether the genes underwent 
#            Up regulation (Up_Regulated) or Down Regulation (Down_Regulated) for every subject. The next three
#            variables record the maximum p-value (Max_p_value), the maximum FDR (Max_FDR) and the minimum
#            absolute log fold change (Min_Abs_LogFC) for each gene that underwent gene regulation.
#
##################################################################################################################


lca <- function(logFC, pVal, FDR, thresholds) {
  # Find the number of subjects and the number of genes
  NumSubj <- ncol(logFC)
  NumGenes <- nrow(logFC)
  
  # set up a data frame for the final results
  RegGenes <- matrix(NA, nrow = NumGenes, ncol = 5)
  RegGenes <- as.data.frame(RegGenes)
  # Set the row names and column names
  row.names(RegGenes) <- row.names(logFC)
  colnames(RegGenes) <- c('Up_Regulated', 'Down_Regulated', 'Max_p_value', 'Max_FDR', 'Min_Abs_LogFC')
  
  # Determine the genes that up regulated based on the thresholds
  RegGenes$Up_Regulated <- rowSums(logFC > thresholds$logFC & pVal < thresholds$pVal & FDR < thresholds$FDR) == NumSubj 
  
  # Determine the genes that down regulated based on the thresholds
  RegGenes$Down_Regulated <- rowSums(logFC < thresholds$logFC & pVal < thresholds$pVal & FDR < thresholds$FDR) == NumSubj 
  
  # Create a mask for any of the genes that underwent regulation
  Mask <- RegGenes$Up_Regulated | RegGenes$Down_Regulated
  
  # For the regulated genes determine the maximum p-value
  RegGenes$Max_p_value[Mask] <- apply(pVal[Mask,], 1, max)
  
  # For the regulated genes determine the maximum FDR
  RegGenes$Max_FDR[Mask] <- apply(FDR[Mask,], 1, max)
  
  # For the up regulated genes determine the absolute minimum logFC 
  RegGenes$Min_Abs_LogFC[RegGenes$Up_Regulated] <- apply(logFC[RegGenes$Up_Regulated,], 1, min)
  
  # For the down regulated genes determine the absolute minimum logFC 
  RegGenes$Min_Abs_LogFC[RegGenes$Down_Regulated] <- apply(logFC[RegGenes$Down_Regulated,], 1, max)
  
  return(RegGenes)
}