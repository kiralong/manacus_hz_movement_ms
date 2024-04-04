#!/usr/bin/env Rscript

# Calculate interspecific heterozygosity from Manacus hybrids
setwd('/path/to/working/dir')

# Load input packages
library(genetics)
library(introgress)
library(ggplot2)

# ===============
# Genotypes table
# ===============
geno_f <- './genotypes_tbl.csv'
geno_df <- read.delim(geno_f, sep=',', header=F)

# ==========
# Loci table
# ==========
loci_f <- './loci_tbl.csv'
loci_df <- read.delim(loci_f, sep=',')

# ================
# Output file name
# ================
out_f <- './intersp_het.csv' 

# ======================================================
# Split the base genotypes df across hybrids and parents
# ======================================================

parent1_id <- 'P1'
parent2_id <- 'P2'

# Get parent 1 genotype data
parent1_df <- geno_df[-c(1,2),c(geno_df[1,] == parent1_id)]
# parent1_df <- geno_df[,c(geno_df[1,] == parent1_id)]
parent1_sams <- as.character(geno_df[2,c(geno_df[1,] == parent1_id)])
  
# Get parent 2 genotype data
parent2_df <- geno_df[-c(1,2),c(geno_df[1,] == parent2_id)]
parent2_sams <- as.character(geno_df[2,c(geno_df[1,] == parent2_id)])

# Get hybrids genotype data
hybrids_df <- geno_df[,c(geno_df[1,] != parent1_id)]
hybrids_df <- hybrids_df[,c(hybrids_df[1,] != parent2_id)]
hybrids_sams <- as.character(hybrids_df[2,])

# ==================================
# Run initial formatting of the data
# ==================================

fmt_genos <- prepare.data(
  admix.gen = hybrids_df,
  parental1 = parent1_df,
  parental2 = parent2_df,
  loci.data = loci_df,
  pop.id = TRUE,
  ind.id = TRUE,
  fixed = FALSE
)

# ======================================
# Calculate interspecific heterozygosity
# ======================================

interspec_het <- calc.intersp.het(fmt_genos)

# =============================
# Format the final output table
# =============================

# Merge all the samples
all_samples <- c(parent1_sams, hybrids_sams, parent2_sams)

# Get a vector of population ids
all_pops <- sapply(strsplit(all_samples, '_'),
                   function(all_samples){all_samples[2]})

# Merge all the insterspecies het values, parents are 0.0 
all_hets <- c(rep(0.0, length(parent1_sams)),
              interspec_het,
              rep(0.0, length(parent2_sams)))

# Export into a CSV
out_df <- data.frame(all_samples, all_pops, all_hets)
colnames(out_df) <- c('SampleID', 'PopID', 'InterSpHet')
out_df <- out_df[order(out_df$PopID),]
out_df$InterSpHet <- out_df$InterSpHet
write.csv(out_df, out_f, quote=F, eol='\n', row.names=F)

# Make plot to quickly check the results
out_df$Index <- 1:nrow(out_df)
plt <- ggplot(data=out_df, aes(x=Index, y=InterSpHet, color=PopID)) +
  geom_point(alpha=0.6) +
  theme_bw()
plt

ggsave('./hybrid_sims.intersp_het.pdf', plot=plt, device='pdf',
       height=4, width=6)
