## Rscript for stats on geographic cline centers and moving genetic clines estimated from HZAR

# Install required packages
# install.packages("ggplot2")

## Load required libraries
library(ggplot2)
library(dplyr)
library(data.table)

## Set working directory
setwd("path/to/hzar/Calculate_parental_allele_freq_diffs")

## Set variables for different runs/filter schemes of data
#max_center_ci <- 11.56
max_center <- 92.5
min_center <- 0
river_position <- 75.125
min_parent_allele_diff <- 0.1
## Note that the average distance between populations in the Manacus dataset is 11.56, and the median is 10.53


## Read in the filtered data tables of the cline parameters for both temporal datasets for comparison
# Import the contemporary dataset (Long)
Long_cline_parameters_table <- read.delim("./long.p8.out.snp_cline_parameters_with_parental_allele_freq_diff.tsv", stringsAsFactors = FALSE)

# Add labels for ease of future graphing
Long_cline_parameters_table$experiment <- "Long"
Long_cline_parameters_table$date_range <- "Contemporary"

# Calulate differenct in credibility interval and add to dataframe
Long_cline_parameters_table$center_ci <- Long_cline_parameters_table$center_high - Long_cline_parameters_table$center_low

# Optional, filter to only keep snps with a CI less than the maximum CI specified above (max_center_ci)
#Long_cline_parameters_table <- subset(Long_cline_parameters_table,center_ci<=max_center_ci)

# Optional, filter to only keep snps within the sampled area of min_center and max_center values specified above
Long_cline_parameters_table <- subset(Long_cline_parameters_table,center>=min_center & center<=max_center)

# Optional, filter by parental species allele frequency difference (note you have to run parental_allele_frequency_diff.py script to add this column to your cline paramter file)
long_diagnosic_cline_table <- subset(Long_cline_parameters_table, parent_allele_freq_diff>=min_parent_allele_diff)

# Import the historic dataset (Brum)
Brum_cline_parameters_table <- read.delim("./brum.p8.out.snp_cline_parameters_with_parental_allele_freq_diff.tsv", stringsAsFactors = FALSE)

# Add labels for ease of future graphing
Brum_cline_parameters_table$experiment <- "Brum"
Brum_cline_parameters_table$date_range <- "Historic"

# Calulate differenct in credibility interval and add to dataframe
Brum_cline_parameters_table$center_ci <- Brum_cline_parameters_table$center_high - Brum_cline_parameters_table$center_low

# Optional, filter to only keep snps with a CI less than the maximum CI specified above (max_center_ci)
#Brum_cline_parameters_table <- subset(Brum_cline_parameters_table,center_ci<=max_center_ci)

# Optional, filter to only keep snps within the sampled area of min_center and max_center values specified above
Brum_cline_parameters_table <- subset(Brum_cline_parameters_table,center>=min_center & center<=max_center)

# Optional, filter by parental species allele frequency difference (note you have to run parental_allele_frequency_diff.py script to add this column to your cline parameter file)
brum_diagnosic_cline_table <- subset(Brum_cline_parameters_table, parent_allele_freq_diff>=min_parent_allele_diff)


## Merge the two data tables for ONLY clines with max CI in BOTH datasets for 1 to 1 comparison
# all_cline_parameters_filtered <- merge(Long_cline_parameters_table, Brum_cline_parameters_table, by.x = "snp_id", by.y = "snp_id")

# If running for diagnostic parental alleles, use this one to merge the data tables with clines in both datasets
all_cline_parameters_filtered <- merge(long_diagnosic_cline_table, brum_diagnosic_cline_table, by.x = "snp_id", by.y = "snp_id")


# Find which clines have overlapping CIs (adds TRUE or FALSE to column in dataframe)

all_cline_parameters_filtered$overlaps <- rep(FALSE,nrow(all_cline_parameters_filtered))

for(i in 1:nrow(all_cline_parameters_filtered)){
  snp_row <- all_cline_parameters_filtered[i,]
  low.x <- snp_row$center_low.x
  high.x <- snp_row$center_high.x
  low.y <- snp_row$center_low.y
  high.y <- snp_row$center_high.y
  # mark if overlapping - assume it's not overlapping, mark as true when found
  if ((low.x >= low.y & low.x <= high.y)         # The bottom part of x is overlapped within y
      | (high.x >= low.y & high.x <= high.y)     # The high part of x is overlapped with y
      | (low.y >= low.x & low.y <= high.x)       # The bottom part of y is overlapped within x
      | (high.y >= low.x & high.y <= high.x)){   # The high part of y is overlapped with y
    all_cline_parameters_filtered[i,"overlaps"] <- TRUE
  }
}

#Add the direction of the movement to main all clines dataframe
all_cline_parameters_filtered <- all_cline_parameters_filtered %>%
  #filter(best_model.x == best_model.y) %>%
  mutate(direction = case_when(
    center.x > center.y ~ "toward_candei",
    center.x < center.y ~ "toward_vitellinus",
    center.x == center.y ~ "none"
  ), magnitude = center.x - center.y)

# Get overlaped data
overlap_tally <- all_cline_parameters_filtered %>%
                select(snp_id,overlaps) %>%
                count(overlaps)

displaced_clines_percentage <- overlap_tally[overlap_tally$overlaps == FALSE,"n"]/sum(overlap_tally$n)

# Look at moving clines
chr_displaced_clines <- all_cline_parameters_filtered %>%
                        select(snp_id,chromosome.x,overlaps) %>%
                        filter(overlaps == FALSE) %>%
                        count(chromosome.x) %>%
                        arrange(desc(n))

# Create a data frame for just the cline center data
moving_cline_parameters <- all_cline_parameters_filtered %>% filter(overlaps == FALSE)
# Export table of "moving" clines
write.csv(moving_cline_parameters, file = "./moving_cline_parameters.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Create a data frame for cline center data
date_range  <- c(all_cline_parameters_filtered$date_range.x,all_cline_parameters_filtered$date_range.y)
cline_center <- c(all_cline_parameters_filtered$center.x,all_cline_parameters_filtered$center.y)
cline_data=data.frame(date_range, cline_center)

# Calculate summary stats for both datasets and make a table
summary_cline_table <- all_cline_parameters_filtered %>%
  select(center.x,center.y) %>%
  na.omit() %>%
  summarise(brum_mean = mean(center.y),
            brum_median = median(center.y),
            brum_sd = sd(center.y),
            long_mean = mean(center.x),
            long_median = median(center.x),
            long_sd = sd(center.x))

# How many clines moved past the river?
past_river_clines <- cline_data %>%
                      filter(cline_center >= river_position) %>%
                      count(date_range)

# How many clines between sites 8 and 9?
genomic_center_clines <- cline_data %>%
                          filter(cline_center >= 20.75) %>%
                          filter(cline_center <= 29.5) %>%
                          count(date_range)
genomic_center_clines$percent <- genomic_center_clines$n/nrow(all_cline_parameters_filtered)

# Direction of the moving clines
direction_moving_clines <- all_cline_parameters_filtered %>%
                          select(snp_id,center.x, center.y, overlaps) %>%
                          filter(overlaps == FALSE) %>%
                          mutate(direction = case_when(
                            center.x > center.y ~ "toward_candei",
                            center.x < center.y ~ "toward_vitellinus",
                            center.x == center.y ~ "none"
                          )) %>%
                          count(direction)

# Export table of scaffolds with "moving" clines
write.csv(chr_displaced_clines, file = "./chr_displaced_clines.csv", quote = FALSE, row.names = FALSE, col.names = TRUE)

