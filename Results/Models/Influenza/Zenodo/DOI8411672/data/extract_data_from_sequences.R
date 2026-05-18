#' Extract phylogenetic and physiochemical data from sequences
#'
#' Function saves a data frame of phylogenetic and physiochemical measurements
#' associated with sequences and their phylogeny to the path given by the
#' user. This data is what is used in antigenic_coordinate_prediction.py
#'
#' @param seqs Seqinr object with sequence data
#' @param tr Tree object associated with sequence data
#' @param file File path where data should be saved
#' 
#' @return None
#'
library('seqinr')
library('adephylo')
library('Peptides')
library('readr')

extract_data_from_sequences <- function(seqs = sequences, tr = tree) {
    df <- data.frame(stringsAsFactors = FALSE)
    
    # Use the published map if it is already loaded in the environment
    if (is.null(smith_map)) {
        ag_map <- make.acmap (table = titer_data, 
                              number_of_dimensions = 2,
                              number_of_optimizations = 1000)
    } else {
        ag_map <- smith_map
    }
    
    # Get the phylogenetic branch length from each sample to the preceding bifurcation
    branch_lengths <- (
        setNames(tr$edge.length[sapply(1:length(tr$tip.label),
                                       function(x, y) which(y == x),
                                       y = tr$edge[, 2])], tr$tip.label)
    )
    
    for (i in seq_len(seqs$nb)) {
        
        # Get physiochemical measurements
        aa_stats <- AAstat(toupper(s2c(seqs$seq[[i]])), plot = FALSE)

        if (seqs$nam[i] %in% ag_map$agNames) {
            # map_index - index of sample in the data that holds the antigenic map
            map_index <- which(rownames(agCoords(ag_map)) == seqs$nam[i])
            # x- and y- coordinates of the current sequence
            current_x <- agCoords(ag_map)[map_index, 1]
            current_y <- agCoords(ag_map)[map_index, 2]
        } else {
            current_x <- 0
            current_y <- 0
        }
        # Hold current data in temp dataframe
        tmp_df <- data.frame("seq_id" = seqs$nam[i],
                             "antigenic_x" = current_x,
                             "antigenic_y" = current_y,
                             #"cluster" = seqs$cluster[i],
                             "distRoot" = distRoot(tr, tips = which(tr$tip.label == seqs$nam[i])),
                             "branch_length" = branch_lengths[seqs$nam[i]],
                             "hydrophobicity" = hydrophobicity(seqs$seq[i]),
                             "charge" = charge(seqs$seq[i]),
                             "boman" = boman(seqs$seq[i]),
                             "instability" = instaIndex(seqs$seq[i]),
                             "isoelectric_point" = aa_stats$Pi,
                             stringsAsFactors = FALSE)
        # Bind tmp_df to all other data
        df <- rbind(df, tmp_df)
    }
    # save data to file
    write_csv(x = df, path = "data/processedData/sequence_physioproperties2.csv")
}
# Run function when file is sourced
extract_data_from_sequences()

