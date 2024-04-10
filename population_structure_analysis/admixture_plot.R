# Script to graph a custom admixture plot

# Load required libraries
library(reshape)
library(ggplot2)
library(forcats)
library(ggthemes)
library(dplyr)
library(magrittr)

setwd("~/path/to/work/dir/Admixture/Long_10k")

#
# For K=2
#
# Load the input file
admix <- read.delim('./admixture_K2.tsv', header = FALSE, stringsAsFactors = FALSE)
admix <- admix[,c(1,2,4,3)]
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

admix.plot <-
  ggplot(madmix, aes(x = sample, y = value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")

admix.plot

num_bars <- 152 + 9
bars_per_inch <- 328/12
wid <- num_bars/(bars_per_inch)

ggsave('./admix_plot_long_k2_10k.pdf', admix.plot,height=3,width=wid)



#
# For K=3
#

# Load the input file
admix <- read.delim('./admixture_K3.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
madmix <-
  madmix %>%
  group_by(variable) %>%
  mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=3", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot

ggsave('./admix_plot_k3_long_10k_reordered.pdf', admix.plot,height=3,width=12)

#  Resize  the plot to compare more easily with the Long plots
num_bars <- 152 + 9
bars_per_inch <- 328/12
wid <- num_bars/(bars_per_inch)

ggsave('./admix_plot_k3_long_10k_reordered.pdf', admix.plot,height=3,width=wid)

#
# For K=4
#

# Load the input file
admix <- read.delim('./admixture_K4.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3','k4')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
# madmix <-
#   madmix %>% 
#   group_by(variable) %>% 
#   mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=4", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot


ggsave('./admix_plot_k4_long_no9.5_10k.pdf', admix.plot,height=3,width=12)

#
# For K=5
#

# Load the input file
admix <- read.delim('./admixture_K5.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3','k4','k5')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
madmix <-
  madmix %>% 
  group_by(variable) %>% 
  mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=5", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot


ggsave('./admix_plot_k5_brum.pdf', admix.plot,height=3,width=12)

#
# For K=6
#

# Load the input file
admix <- read.delim('./admixture_K6.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3','k4','k5','k6')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
madmix <-
  madmix %>% 
  group_by(variable) %>% 
  mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=6", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot


ggsave('./admix_plot_k6_brum.pdf', admix.plot,height=3,width=12)

#
# For K=7
#

# Load the input file
admix <- read.delim('./admixture_K7.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3','k4','k5','k6','k7')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
madmix <-
  madmix %>% 
  group_by(variable) %>% 
  mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=7", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot


ggsave('./admix_plot_k7_brum.pdf', admix.plot,height=3,width=12)

#
# For K=8
#

# Load the input file
admix <- read.delim('./admixture_K8.tsv', header = FALSE, stringsAsFactors = FALSE)
# Add thee column names
colnames(admix) <- c('pop','sample','k1','k2','k3','k4','k5','k6','k7','k8')
# Add the population IDs
admix$loc <- sapply(strsplit(admix$sample, '_'), function(admix){admix[1]})
# Sort by population and K val
admix <- admix[
  with(admix, order(pop,k2)),
]
# Melt so single column for Ks
madmix <- melt(admix, id=c('sample','pop', 'loc'))

# Reorder by k-proportions
madmix <-
  madmix %>% 
  group_by(variable) %>% 
  mutate(sample = fct_reorder(sample, desc(value)))

admix.plot <-
  ggplot(madmix, aes(factor(sample), value, fill = factor(variable))) +
  geom_col() +
  facet_grid(~fct_inorder(loc), scales = "free", switch = "x", space = "free") +
  theme_minimal() + labs(x = "Individuals", title = "K=8", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(
    panel.spacing.x = unit(0.001, "lines"),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank()
  ) +
  scale_fill_gdocs(guide = "none")
admix.plot


ggsave('./admix_plot_k8_brum.pdf', admix.plot,height=3,width=12)
