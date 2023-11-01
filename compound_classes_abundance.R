#!/usr/bin/Rscript

# This script finds compound classes whose abundance is different between 2 groups
# N.B. --group_2 is used as baseline
# N.B. The peaks file should be a tab-delimited text file where columns are compounds and rows
# individual samples. The first two columns must be sample names and groups, in that order.
# N.B. The compound classes file should be a tab-delimited text file where the first column is the
# compound name/symbol (matching the one in the peaks file) and the second is the compound class.

### ---------------------------------------- ###

parseArgs <- function() {
  
  # Read command line arguments
  args <- commandArgs()
  
  # Peaks file
  peaks_file <- args[match("--peaks_file", args) + 1]
  
  # Compound classes file
  compound_classes_file <- args[match("--compound_classes_file", args) + 1]
  
  # Group 1
  group_1 <- args[match("--group_1", args) + 1]
  
  # Group 2
  group_2 <- args[match("--group_2", args) + 1]
  
  # Paired samples parameter (assumes files are ordered)
  if('--paired' %in% args) {
    
    paired_toggle <- TRUE
    
  } else {

    paired_toggle <- FALSE

  }
  
  return(c(peaks_file, compound_classes_file, group_1, group_2, paired_toggle))
  
}

### ---------------------------------------- ###

fillZeros <- function(raw) {
  
  # Convert NAs to 0
  raw[is.na(raw)] <- 0
  
  # Remove columns with no data
  raw <- raw[, c(TRUE, TRUE, colSums(raw[,3:ncol(raw)]) > 0)]
  
  # Convert 0 values to 1/5 of the min of that compound
  for(i in 3:NCOL(raw)) {
    
    if(0 %in% raw[,i]) {
      
      raw[raw[,i] == 0, i] <- min(raw[raw[,i] > 0, i]) / 5
      
    }
    
  }
  
  return(raw)
  
}

### ---------------------------------------- ###

correctLibrarySize <- function(raw) {
  
  total_peaks <- rowSums(raw[, 3:ncol(raw)])
  total_peaks_median <- median(total_peaks)
  norm_factors <- total_peaks_median / total_peaks
  
  raw[, 3:ncol(raw)] <- raw[, 3:ncol(raw)] * norm_factors
  
  return(raw)
  
}

### ---------------------------------------- ###

getStatistics <- function(data, classes, group_1, group_2, paired_toggle, output_name) {
  
  # Extract groups indexes
  group_1_ids <- data[,2] == group_1
  group_2_ids <- data[,2] == group_2
  
  # Removing classes with only one compound species
  good_classes <- NULL
  for(cl in classes[! duplicated(classes[,2]), 2]) {
    
    compounds <- length(classes[classes[,2] == cl, 1])
    
    if(compounds > 1) {
      
      good_classes <- c(good_classes, cl)
      
    }
    
  }
  
  # Init results table
  stats <- data.frame(CompoundClass = good_classes, Group_1_Abundance = rep(0, length(good_classes)), Group_2_Abundance = rep(0, length(good_classes)),
                      Log2FC = rep(0, length(good_classes)), pval = rep(0, length(good_classes)), FDR = rep(0, length(good_classes)))
  
  # Compare groups
  for(cl in good_classes) {
    
    if(sum(colnames(data) %in% classes[classes[,2] == cl, 1]) == 1) {
      
      group_1_abundance <- data[group_1_ids, colnames(data) %in% classes[classes[,2] == cl, 1]]
      
    } else {
      
      group_1_abundance <- rowSums(data[group_1_ids, colnames(data) %in% classes[classes[,2] == cl, 1]])
      
    }
    group_1_average <- mean(group_1_abundance)
    
    if(sum(colnames(data) %in% classes[classes[,2] == cl, 1]) == 1) {
      
      group_2_abundance <- data[group_2_ids, colnames(data) %in% classes[classes[,2] == cl, 1]]
      
    } else {
      
      group_2_abundance <- rowSums(data[group_2_ids, colnames(data) %in% classes[classes[,2] == cl, 1]])
      
    }
    group_2_average <- mean(group_2_abundance)
    
    log2fc <- log(group_1_average / group_2_average, 2)
    pval <- t.test(group_1_abundance, group_2_abundance, alternative = "two.sided", var.equal = T, paired = paired_toggle)$p.value
    stats[stats$CompoundClass == cl, 2:5] <- c(group_1_average, group_2_average, log2fc, pval)
    
  }
  
  # Remove classes not present in either group
  stats <- subset(stats, (Group_1_Abundance > 0) | (Group_2_Abundance > 0))
  
  # FDR correction
  fdr <- p.adjust(stats$pval, method="fdr")
  stats$FDR <- fdr
  
  # Exporting results
  write.table(stats, output_name, sep = "\t", row.names = F)
  
  return(stats)
  
}

### ---------------------------------------- ###

plotStats <- function(stats, output_name) {
  
  # Replace Inf and -Inf values
  stats[stats$Log2FC == Inf, 'Log2FC'] <- max(stats[stats$Log2FC != Inf, 'Log2FC']) + 1
  stats[stats$Log2FC == -Inf, 'Log2FC'] <- min(stats[stats$Log2FC != -Inf, 'Log2FC']) - 1
  
  rel_abundance <- stats$Group_1_Abundance / sum(stats$Group_1_Abundance)
  rel_abundance <- rel_abundance * 10 / max(rel_abundance) # Scaling to max 10 size
  
  min_y <- floor(min(stats$Log2FC))
  max_y <- floor(max(stats$Log2FC)) + 2
  
  stats_plot <- ggplot(stats, aes(x = CompoundClass, y = Log2FC, fill = FDR))+
    geom_hline(yintercept = 0, linewidth = 0.5, color = "gray")+
    scale_fill_gradient(low = "red", high = "white", limits = c(min(stats$FDR), 0.05))+
    ylim(min_y, max_y)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
          axis.title.y = element_text(face = "bold"),
          axis.title.x = element_text(face = "bold"))+
    geom_point(x = stats$CompoundClass, y = max_y - 0.75,
               shape = 21, size = rel_abundance, stroke = 1, color = "black", fill = "green")+
    geom_col(color = "black", size = 1, width = 0.8)
  ggsave(output_name, stats_plot, dpi = 300, width = NROW(stats) * 0.5 + 2, height = 7)
  
}

### ------------------MAIN------------------ ###

library(ggplot2)

# Import parameters
parameters <- parseArgs()

comparison_name <- paste(parameters[3], parameters[4], sep="_vs_")
comparison_name <- gsub(' ', '_', comparison_name)

# Import peaks counts
peaks <- read.delim(as.character(parameters[1]), header=TRUE, sep="\t", check.names=FALSE)

# Change spaces in group names to underscores
peaks[,2] <- gsub(' ', '_', peaks[,2])

# Select groups of interest
if((! parameters[3] %in% peaks$Group) | (! parameters[4] %in% peaks$Group)) {
  
  print("ERROR: wrong categories")
  quit(save="no", status=1)
  
} else {
  
  peaks <- peaks[(peaks$Group == parameters[3]) | (peaks$Group == parameters[4]), ]
  
}

# Preprocessing
peaks <- fillZeros(peaks) # Dealing with 0 values
peaks <- correctLibrarySize(peaks) # Correct for total peaks in samples

# Import compound classes
class_type <- gsub(".tsv", "", basename(parameters[2]))
compound_classes <- read.delim(parameters[2], stringsAsFactors = F, check.names = F, sep = "\t")
compound_classes <- compound_classes[order(compound_classes[,2]),]

# Compare compound classes
output_name <- paste(comparison_name, '_compound_classes_comparison_', class_type, '.tsv', sep='')
stats <- getStatistics(peaks, compound_classes, parameters[3], parameters[4], parameters[5], output_name)

# Plot classes stats
output_name <- paste(comparison_name, '_compound_classes_comparison_', class_type, '.png', sep='')
plotStats(stats, output_name)
