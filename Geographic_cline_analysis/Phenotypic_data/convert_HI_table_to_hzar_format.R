# Convert Hybrid Index table from gghybrid into a hzar readable format

#install.packages("dplyr")

# Load required packages
library(dplyr)

#Set Working Directory
setwd("/path/to/hybrid_index_cline/allelic_format_test")

# Read in your data table, a hybrid index table that is the output file from R package gghybrid
HI_table <- read.delim("./HI_table_Long.txt", stringsAsFactors = FALSE)
#HI_table <- read.delim("./HI_table_Brum.txt", sep = " ", stringsAsFactors = FALSE)

# Remove the 9th contemporary pop Miramar, if desired
HI_table <- subset(HI_table, hi.POPID != "095MR")

# Get list of unique populatin IDs from the input Hybrid Index Table
popids <- rev(unique(HI_table$hi.POPID))

# Initialize empty vectors to put in data for final dataframe
n_counts <- c()
mean_HI <- c()
distance <- rep(0, length(popids))

# Loop through file to calcule mean and count number of individuals by population
for (pop in popids){
  pop.data <- subset(HI_table,hi.POPID == pop)
  n <- nrow(pop.data) *2
  avg_HI <- mean(pop.data$hi.h_posterior_mode)
  n_counts <- c(n_counts,n)
  mean_HI <- c(mean_HI, avg_HI)
}

# Get the inverse of the hybrid index as your "minor allele"
inverse_HI <- 1 - mean_HI

# Make final dataframe with desired vectors
hzar_input <- data.frame(popids, distance, mean_HI, inverse_HI, n_counts)

# Rename columns
colnames(hzar_input) <- c("Population","Distance","69_1.A","69_1.B","69_1.N")

# Save as a csv file
write.csv(hzar_input, "./HI_Long_p8.hzar.csv", row.names = FALSE)
