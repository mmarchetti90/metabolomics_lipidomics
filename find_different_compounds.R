#!/usr/bin/Rscript

# This script finds compounds (lipids or metabolites) whose abundance is different between 2 groups
# N.B. --group_2 is used as baseline
# N.B. The peaks file should be a tab-delimited text file where columns are compounds and rows
# individual samples. The first two columns must be sample names and groups, in that order.

# Comparable to MetaboAnalyst, but in plain code

### ---------------------------------------- ###

parseArgs <- function() {

	# Read command line arguments
	args <- commandArgs()

	# Peaks file
	peaks_file <- args[match("--peaks_file", args) + 1]
	
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
	
	# Number of compounds to plot in heatmap
	if('--heatmap_compounds' %in% args) {
	  
	  heatmap_compounds <- as.numeric(args[match("--heatmap_compounds", args) + 1])
	  
	} else {
	  
	  heatmap_compounds <- 25
	  
	}

	return(c(peaks_file, group_1, group_2, paired_toggle, heatmap_compounds))

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
  
  # Corrects for differing sum of peaks across samples

  total_peaks <- rowSums(raw[, 3:ncol(raw)])
  total_peaks_median <- median(total_peaks)
  norm_factors <- total_peaks_median / total_peaks
  
  raw[, 3:ncol(raw)] <- raw[, 3:ncol(raw)] * norm_factors
  
  return(raw)
  
}

### ---------------------------------------- ###

normalizeData <- function(raw) {
  
  raw[,3:NCOL(raw)] <- log(raw[,3:NCOL(raw)], 10) # Log-transform
  
  for(i in 3:NCOL(raw)) { # Pareto-scaling (mean centering and dividing by sqrt(sd))
    
    raw[,i] <- (raw[,i] - mean(raw[,i])) / sqrt(sd(raw[,i]))
    
  }
  
  return(raw)
  
}

### ---------------------------------------- ###

pcaTransform <- function(data) {
  
  # PCA
  pca <- prcomp(t(data[, 3:ncol(data)]))
  
  # Screeplot
  sdev <- pca$sdev
  sdev <- sdev * 100 / sum(sdev)
  sdev <- data.frame("StDev" = sdev,
                     "PC" = seq(1, length(sdev)),
                     check.names=FALSE)
  
  ggplot(data=sdev, aes(x=PC, y=StDev)) +
    geom_bar(stat="identity") +
    labs(x="PC", y="StDev %",) +
    xlim(0, nrow(sdev) + 0.5)
  ggsave("screeplot.png")
  
  # Scatterplot
  pca_coords <- as.data.frame(pca$rotation)
  pca_coords["Group"] <- norm[,2]
  
  ggplot(data=pca_coords, aes(x=PC1, y=PC2, fill=Group)) +
    geom_point(shape=21, size=5, stroke=1, colour="black")
  ggsave("pca_plot.png")
  
}

### ---------------------------------------- ###

analyzeData <- function(raw, norm, params) {
  
  # Init results table
  results <- data.frame(Compound = colnames(peaks[,3:NCOL(peaks)]), Log2FC = rep(0, NCOL(peaks) - 2),
                        pval = rep(1, NCOL(peaks) - 2), stringsAsFactors = F)
  
  # Extract groups indexes
  group_1_ids <- norm[,2] == params[2]
  group_2_ids <- norm[,2] == params[3]
  
  # Analysis
  for(p in 3:NCOL(norm)) {
    
    g1_norm <- norm[group_1_ids, p]
    g2_norm <- norm[group_2_ids, p]
    g1_raw <- raw[group_1_ids, p]
    g2_raw <- raw[group_2_ids, p]
    log2fc <- log(mean(g1_raw) / mean(g2_raw), 2) # MetaboAnalyst calculates fold change on data before log-transform and scaling (but after sun normalization)
    pval <- t.test(g2_norm, g1_norm, alternative = "two.sided", var.equal = T, paired = params[4])$p.value
    results[p - 2, 2] <- log2fc
    results[p - 2, 3] <- pval
    
  }
  
  # FDR correction
  fdr <- p.adjust(results$pval, method="fdr") 
  results <- cbind(results, fdr)
  colnames(results)[NCOL(results)] <- "FDR"
  results <- results[order(results$FDR, decreasing = F),] # Order by FDR
  
  return(results)
  
}

### ---------------------------------------- ###

plotHeatmap <- function(N, comparison_name, params, norm, stats) {
  
  # Subset dataset
  top <- t(norm[,3:NCOL(norm)])
  colnames(top) <- seq(1, NCOL(top))
  top <- top[rownames(top) %in% stats[1:N,]$Compound,]
  
  # Prepare annotation colors
  annotation_col <- data.frame(SampleType = norm[,2], row.names = colnames(top))
  annotation_colors <- list("SampleType" = c(Group1 = "red", Group2 = "green"))
  names(annotation_colors$SampleType) = params[2 : 3]
  
  # Plot
  pheatmap(top,
           cluster_cols = F,
           scale = "row",
           border_color = NA,
           show_colnames = F,
           annotation_col = annotation_col,
           color = colorRampPalette(c("blue", "white", "red"))(50),
           annotation_colors = annotation_colors,
           filename = paste(comparison_name, "_top", N, "_compounds.png", sep = ""))
  
}

### ------------------MAIN------------------ ###

library(ggplot2)
library(pheatmap)

# Import parameters
parameters <- parseArgs()
comparison_name <- paste(parameters[2], parameters[3], sep="_vs_")
comparison_name <- gsub(' ', '_', comparison_name)

# Import peaks counts
peaks <- read.delim(as.character(parameters[1]), header=TRUE, sep="\t", check.names=FALSE)

# Change spaces in group names to underscores
peaks[,2] <- gsub(' ', '_', peaks[,2])

# Select groups of interest
if((! parameters[2] %in% peaks[,2]) | (! parameters[3] %in% peaks[,2])) {
  
  print("ERROR: wrong categories")
  quit(save="no", status=1)
  
} else {
  
  peaks <- peaks[(peaks[,2] == parameters[2]) | (peaks[,2] == parameters[3]), ]
  
}

# Preprocessing
peaks <- fillZeros(peaks) # Dealing with 0 values
peaks <- correctLibrarySize(peaks) # Correct for total peaks in samples
norm <- normalizeData(peaks) # Normalize data

# PCA
pcaTransform(norm)

# Find significantly different compounds
stats <- analyzeData(peaks, norm, parameters) # pairwise t-test
write.table(stats, paste(comparison_name, '_significant_compounds.tsv', sep=''), row.names=F, sep="\t")

# Plot data
plotHeatmap(as.numeric(parameters[5]), comparison_name, parameters, norm, stats)
