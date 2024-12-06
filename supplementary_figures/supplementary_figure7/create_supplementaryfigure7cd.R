### CREATE SUPPLEMENTARY FIGURE 7C AND 7D #############################################################################
# Immune phenotype distribution by subtype in TCGA
# Immune phenotype distribution for high risk IC subtypes in TCGA 

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

# Read data -----
data <- read.delim(file.path(main_repo_path, 'data','rna_megatable.txt'), check.names = F)
data <- data %>% mutate(Stage = factor(Stage, levels = c("Primary", "Metastatic")), 
                        Immune.phenotype = factor(Immune.phenotype))

source_data_7c <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary")), 
                               !is.na(Immune.phenotype), !is.na(group)) %>% 
  group_by(Immune.phenotype, group) %>% summarize(n = n()) %>% group_by(group) %>% 
  mutate(total = sum(n), Proportion = n/total) 

source_data_7d <- data %>% filter(((Cohort == "TCGA")&(Stage == "Primary")), 
                               !is.na(Immune.phenotype), group == "ER+ High") %>% 
  mutate(eniclust = gsub("ic", "IC", eniclust)) %>% 
  group_by(Immune.phenotype, eniclust) %>% summarize(n = n()) %>% group_by(eniclust) %>% 
  mutate(total = sum(n), Proportion = n/total) 

## Save plot -----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure7C.pdf"), height = 4, width = 6)
ggbarplot(source_data_7c, "group", "Proportion", color = "Immune.phenotype", 
          fill = "Immune.phenotype", palette = immunophenotypes_colors, 
          ylab = "Proportion of samples by immune phenotype", xlab = "",
          title = "Primary (TCGA)") +
  geom_text(aes(label = total, y = 1.03, size = 13), check_overlap = T) +
  guides(size = "none", 
         fill = guide_legend(title = "Immune\nphenotype"), 
         color = guide_legend(title = "Immune\nphenotype")) +
  labs(fill = "Immune phenotype", color = "Immune phenotype") +
  theme_LM  
dev.off()

pdf(file.path(main_repo_path, "plots", "SupplementaryFigure7D.pdf"), height = 4, width = 6)
ggbarplot(source_data_7d, "eniclust", "Proportion", color = "Immune.phenotype", 
          fill = "Immune.phenotype", palette = immunophenotypes_colors, 
          ylab = "Proportion of samples by immune phenotype", xlab = "",
          title = "Primary (TCGA)") +
  geom_text(aes(label = total, y = 1.03), check_overlap = T) +
  labs(fill = "Immune phenotype", color = "Immune phenotype") +
  theme_LM
dev.off()

### SAVE ##########################################################################################
source_data_7c %>% 
  select(Immune_Phenotype = Immune.phenotype, Subtype = group, Total_Samples = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure7c_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")

source_data_7d %>% 
  select(Immune_Phenotype = Immune.phenotype, ENiClust = eniclust, Total_Samples = total, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure7d_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")

