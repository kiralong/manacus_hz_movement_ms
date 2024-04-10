# Script for Manacus RAD PCA

setwd("/path/to/working/dir/PCA/maf01_run")

# Load required libraries
library(adegenet)
library(RColorBrewer)
library(ggplot2)

# load GENEPOP object NOTE THAT THE GENEPOP FILE MUST BE RENAMED TO .GEN 
mana = import2genind('./populations.snps.10k.gen')

# Load sample master list with population and dataset information to assign symbols for graphing
sample_list <- read.delim('./pop_labels.tsv')
sample_list$pca.labels <- paste(sample_list$pop,sample_list$experiment, sep = '_')
sample_list$symbol <- 16
sample_list$symbol[sample_list$experiment == "Brum"] <- 17

# Get population list
pop(mana) = sample_list$pca.labels

# scale genind object
mana_s = scaleGen(mana, NA.method="mean")

#Get random distinct colors in a pallete from Rcolorbrewer (only need to do this to generate palette first time)
#palette3_info <- brewer.pal.info[brewer.pal.info$category == "div", ]   # Extract color info
#palette3_all <- unlist(mapply(brewer.pal,                               # Create vector with all colors
#                              palette3_info$maxcolors,
#                              rownames(palette3_info)))
#palette3_all                                                            # Print all colors

#set.seed(9845739)                                                       # Set random seed
#palette3 <- sample(palette3_all, 9)                                     # Sample colors
#palette3                                                                # Print hex color codes
# Set color palette of choice
palette3 <- c("#762A83","#313695","#7F3B08","#1A9850","#F1B6DA","#D9EF8B","#A50026","#3288BD","#003C30")

# Do the PCA
mana_pca = dudi.pca(mana_s, scannf=FALSE)

# Obtain Eigenvalues and PC components %
eigTotal = sum(mana_pca$eig)
pc1 <- (mana_pca$eig[1]/eigTotal) * 100
pc2 <- (mana_pca$eig[2]/eigTotal) * 100
pc3 <- (mana_pca$eig[3]/eigTotal) * 100
pc4 <- (mana_pca$eig[4]/eigTotal) * 100
pc5 <- (mana_pca$eig[5]/eigTotal) * 100


# Capturing PC1 and PC2 from adegenet into sample list for plotting with ggplot
sample_list$pc1 <- mana_pca$li$Axis1
sample_list$pc2 <- mana_pca$li$Axis2

# Make PCA scatter plot in ggplot2
pca_fig <- ggplot(sample_list, aes(x=pc1, y=pc2, color=pop, shape=experiment)) +
  
  # Add origin lines
  geom_hline(yintercept =0, linetype="dashed", color = "gray") +
  geom_vline(xintercept =0, linetype="dashed", color = "gray") +
  
  # Add points
  geom_point(
    size=4,
    alpha=0.65) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_manual(values = c(palette3)) +

# Add titles
  labs(title='Manacus RAD PCA',
     x=paste('PC1 ',format(round(pc1,2)),'%',sep = ''),
     y=paste('PC2 ',format(round(pc2,2)),'%',sep = '')) +
  
# Set Themes
  theme_light(base_size=16) +
  theme(panel.grid = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  
# Centralize title
  theme(plot.title=element_text(hjust=0.5)) +

# Legend
  labs(col="Population ID", shape="Dataset")
  
pca_fig

# Save plot
ggsave('./manacus_rad_pca_maf01.pdf', plot=pca_fig, width=8, height=6)


### Plot PC1 versus Hybrid Index value

HI_table_long <- read.delim("/path/to/Long_p9_r80_mac3_10k/HI_table.txt")
HI_table_brum <- read.csv("/path/to/Brum_p9_r80_10K/HI_table.txt", sep="")
HI_table_all <- rbind(HI_table_long,HI_table_brum)
hi_table_simplified <- data.frame(HI_table_all$hi.INDLABEL,HI_table_all$hi.h_posterior_mode)

master_table <- merge(sample_list,hi_table_simplified, by.x = "sample_id", by.y = "HI_table_all.hi.INDLABEL")

pcxhc_fig <- ggplot(master_table, aes(x=pc1, y=HI_table_all.hi.h_posterior_mode, color=pop, shape=experiment)) +
  
  # Add points
  geom_point(
    size=4,
    alpha=0.65) +
  scale_shape_manual(values = c(17,16)) +
  scale_color_manual(values = c(palette3)) +
  
  # Add titles
  labs(title='Manacus RAD PC1 vs HI',
       x=paste('PC1 ',format(round(pc1,2)),'%',sep = ''),
       y="Hybrid Index") +
  
  # Set Themes
  theme_light(base_size=16) +
  
  # Centralize title
  theme(plot.title=element_text(hjust=0.5)) +
  
  # Legend
  labs(col="Population ID", shape="Dataset")

pcxhc_fig

# Save plot
ggsave('./manacus_rad_pc1_vs_HI_maf01.pdf', plot=pcxhc_fig, width=8, height=6)

### Check correlation
correlation_test <- lm(HI_table_all.hi.h_posterior_mode~pc1,data = master_table)
summary(correlation_test)
