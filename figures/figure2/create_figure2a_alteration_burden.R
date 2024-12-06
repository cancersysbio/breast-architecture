### CREATE FIGURE 2A #############################################################################
# Alteration burden in primary vs metastatic samples by subtype

### PREAMBLE #####################################################################################
library(tidyverse)
library(ggpubr)

# Set the main path for repo
main_repo_path <- ""
if ((!exists("main_repo_path")) | main_repo_path == "") {
  stop("Error: Path for main repo not set. Please set main_repo_path <- '/path/to/repo/breast-architecture' and try again.")
}

### MAIN ##########################################################################################
source("plot_configs.R")

# Read data -- generated in alteration_burden.R
alterations_data_path <- file.path(main_repo_path, 'data','alteration_burden.txt')
if(!file.exists(alterations_data_path)) {
  stop("Error: Alteration burden data not found, run alteration_burden.R analysis script. ")
} 

alterations <- read.delim(alterations_data_path)

primary_megatable <- read.delim(file.path(main_repo_path, 'data','primary_megatable.txt'))
metastatic_megatable <- read.delim(file.path(main_repo_path, 'data','metastatic_megatable.txt'))

samples <- rbind(
  primary_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Primary"), 
  metastatic_megatable %>% select(Sample, group, genome_doubled) %>% mutate(Stage = "Metastatic")
)

# Plot proportions by groups --------
total_samples <- samples %>% select(Stage, group) %>% table() %>% data.frame() %>% dplyr::rename(total_samples = 3)

alterations$Stage <- factor(alterations$Stage, levels = c("Primary", "Metastatic"))

source_data <- alterations %>% filter(!is.na(group)) %>% 
  group_by(Stage, group, Alteration) %>% summarize(n = n()) %>% group_by(Stage, group) %>%
  mutate(total = sum(n), Proportion = n/total) %>% left_join(total_samples) %>% 
  mutate(n_by_sample = n/total_samples)

a <- ggbarplot(source_data, "Stage", "n_by_sample", 
               color = "Alteration", fill = "Alteration", palette = alteration_colors,
               ylab = "Alterations per sample", xlab = "", facet.by = "group", ncol = 1)  +
  theme_LM +
  guides(fill = guide_legend(ncol = 3)) +
  rotate_x_text(35) 

# Save plots 
ggsave(file.path(main_repo_path, "plots", "Figure2a_2.svg"), 
       plot = a, height = 14, width = 6)

pdf(file.path(main_repo_path, "plots", "Figure2a_2.pdf"), height = 10, width = 4)
a
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(Stage, Subtype = group, Alteration, Alterations_per_Sample = n_by_sample) %>% 
  write.table(file.path(main_repo_path, "data", "Figure2a_alt_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")
