### CREATE SUPPLEMENTARY FIGURE 7E #############################################################################
# Comparison of CIBERSORTx estimates for TNBC

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

source_data <- data %>% filter(group %in% c("IC10", "IC4ER-"), Cohort %in% c("TCGA", "METABRIC")) %>% 
  pivot_longer(cols = c("CAFs", "T-cells"), names_to = "Cell type", values_to = "Proportion") %>%  
  mutate(Cohort = factor(Cohort, levels = c("TCGA", "METABRIC")))

## Save plot -----
pdf(file.path(main_repo_path, "plots", "SupplementaryFigure7E.pdf"), height = 5, width = 5)
ggboxplot(source_data, "group", "Proportion", color = "group", palette = group_colors, 
          facet.by = c("Cohort", "Cell type"), xlab = "", main = "Estimated proportions from CIBERSORTx") +
  theme_LM +
  rotate_x_text(35) +
  rremove("legend")
dev.off()


### SAVE ##########################################################################################
source_data %>% 
  select(Sample, Subtype = group, Cohort, `Cell type`, Proportion) %>% 
  write.table(file.path(main_repo_path, "data", "SupplementaryFigure7e_sourcetable.txt"), 
              quote = F, row.names = F, sep = "\t")


